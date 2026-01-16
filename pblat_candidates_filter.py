#!/usr/bin/env python
"""
BLAT filter for RNA editing candidates.
Python replacement for pblat_candidates_ln.pl.
This script ensures complete functional equivalence with the original Perl script.
"""

import sys
import os
import argparse
import subprocess
from collections import defaultdict


def run_command(command, cwd=None):
    """
    Runs a command using subprocess and handles errors.
    """
    try:
        result = subprocess.run(
            command,
            shell=True,
            check=True,
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        return True
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"Error executing command: {e.cmd}\n")
        sys.stderr.write(f"Stderr: {e.stderr}\n")
        return False


def prepare_blat_input(infile, bam_path, fa_path, min_base_qual, qual_offset, threads, samtools_exe):
    """
    Reads a variant file, extracts reads from a BAM file supporting the variant,
    and writes them to a FASTA file for PBLAT.
    This function replicates the exact logic from the Perl script.
    """
    print("preparing blat input...")
    
    with open(infile, 'r') as sites_file, open(fa_path, 'w') as fa_file:
        for line in sites_file:
            line = line.strip()
            if not line:
                continue
            
            fields = line.split()
            chrom = fields[0]
            position = int(fields[1])
            edit_nuc = fields[4]
            
            # Create temporary file for samtools output
            temp_name = f"{fa_path}_tmp_{chrom}_{position}"
            bam_position = f"{chrom}:{position}-{position}"
            
            # Run samtools view to extract reads at this position
            cmd = f"{samtools_exe} view -@ {threads} {bam_path} {bam_position} > {temp_name}"
            if not run_command(cmd):
                sys.stderr.write(f"Warning: Failed to extract reads for {chrom}:{position}\n")
                continue
            
            mismatch_read_count = 0
            
            # Process each read from samtools output
            try:
                with open(temp_name, 'r') as tmp_file:
                    for bam_line in tmp_file:
                        bam_line = bam_line.strip()
                        if not bam_line:
                            continue
                        
                        bam_fields = bam_line.split('\t')
                        if len(bam_fields) < 11:
                            continue
                        
                        alignment = int(bam_fields[1])
                        read_start = int(bam_fields[3])
                        cigar = bam_fields[5]
                        sequence = bam_fields[9]
                        qualities = bam_fields[10]
                        
                        sequence_bases = list(sequence)
                        qual_scores = list(qualities)
                        
                        current_pos = read_start
                        read_pos = 1
                        base_readpos = 0
                        
                        # Parse CIGAR string - exact replication of Perl logic
                        import re
                        cigarnums = re.split(r'[MIDNSHP]', cigar)
                        cigarletters = re.split(r'[0-9]+', cigar)
                        cigarletters = cigarletters[1:]  # Remove first empty element
                        
                        for i in range(len(cigarnums)):
                            # Optimization: stop if we've passed the position
                            if current_pos > position:
                                break
                            
                            cigar_num = int(cigarnums[i]) if cigarnums[i] else 0
                            cigar_letter = cigarletters[i] if i < len(cigarletters) else ''
                            
                            # Handle soft clipping (S) or hard clipping (I)
                            if cigar_letter in ['S', 'I']:
                                read_pos = read_pos + cigar_num
                            # Handle deletion (D) or N (skip from reference)
                            elif cigar_letter in ['D', 'N']:
                                current_pos = current_pos + cigar_num
                            # Handle match/mismatch (M)
                            elif cigar_letter == 'M':
                                for j in range(cigar_num):
                                    if current_pos == position:
                                        # Check if this is the edited base with sufficient quality
                                        if (read_pos - 1 < len(sequence_bases) and 
                                            read_pos - 1 < len(qual_scores) and
                                            sequence_bases[read_pos - 1] == edit_nuc and
                                            ord(qual_scores[read_pos - 1]) >= min_base_qual + qual_offset):
                                            base_readpos = 1
                                    
                                    current_pos += 1
                                    read_pos += 1
                        
                        # If this read supports the variant, write to FASTA
                        if base_readpos:
                            fa_file.write(f">{chrom}-{position}-{mismatch_read_count}\n{sequence}\n")
                            mismatch_read_count += 1
                            
            except IOError as e:
                sys.stderr.write(f"Warning: Error processing temp file {temp_name}: {e}\n")
            
            # Clean up temporary file
            try:
                os.remove(temp_name)
            except OSError:
                pass
    
    return True


def run_pblat(refgenome, fa_path, psl_path, threads, pblat_exe):
    """
    Runs PBLAT on the generated FASTA file.
    """
    print(f"blatting reads (threads: {threads})...")
    pblat_params = (
        f"-threads={threads} -stepSize=5 -repMatch=2253 -minScore=20 "
        f"-minIdentity=0 -noHead {refgenome} {fa_path} {psl_path}"
    )
    command = f"{pblat_exe} {pblat_params}"
    return run_command(command)


def process_pblat_results(psl_path, score_limit):
    """
    Processes the PSL output to identify uniquely mapped reads.
    Returns two dictionaries: one for valid sites and one for discarded reads.
    """
    print("processing reads...")
    
    # PSL format: https://genome.ucsc.edu/FAQ/FAQformat.html#format12
    # Fields (0-indexed):
    # 0: matches (score)
    # 9: Q name (query name)
    # 13: T name (target name/chromosome)
    # 17: block count
    # 18: block sizes
    # 20: block starts
    
    psl_hash = defaultdict(list)
    
    with open(psl_path, 'r') as psl_file:
        for line in psl_file:
            line = line.strip()
            if not line:
                continue
            
            psl_fields = line.split('\t')
            if len(psl_fields) < 21:
                continue
            
            name = psl_fields[9]
            # Create blatscore string: matches@Tname@blockcount@blocksizes@blockstarts
            blatscore = f"{psl_fields[0]}@{psl_fields[13]}@{psl_fields[17]}@{psl_fields[18]}@{psl_fields[20]}"
            
            if psl_hash[name]:
                psl_hash[name] = f"{psl_hash[name]}-{blatscore}"
            else:
                psl_hash[name] = blatscore
    
    site_hash = defaultdict(int)
    discard_hash = defaultdict(int)
    
    for psl_key in psl_hash.keys():
        # Parse the key: chr-position-mismatchcount
        split_key = psl_key.split('-')
        site = f"{split_key[0]}_{split_key[1]}"
        
        # Parse all alignments for this read
        psl_lines = psl_hash[psl_key].split('-')
        
        largest_score = 0
        largest_score_line = psl_lines[0]
        score_array = []
        
        for score_line in psl_lines:
            scores_array = score_line.split('@')
            line_score = int(scores_array[0])
            score_array.append(line_score)
            
            if line_score > largest_score:
                largest_score_line = score_line
                largest_score = line_score
        
        # Sort scores in descending order
        score_array.sort(reverse=True)
        score_array[1] = 0 if len(score_array) < 2 else score_array[1]
        
        # Parse the best alignment
        split_largest_line = largest_score_line.split('@')
        overlap_found = 0
        
        # Check if best alignment is on correct chromosome and unique enough
        if (split_largest_line[1] == split_key[0] and 
            score_array[1] < (score_array[0] * score_limit)):
            
            num_blocks = int(split_largest_line[2])
            block_sizes = split_largest_line[3].split(',')
            block_starts = split_largest_line[4].split(',')
            
            # Check if variant position is within any aligned block
            for i in range(num_blocks):
                start_pos = int(block_starts[i]) + 1  # Convert from 0-based to 1-based
                end_pos = int(block_starts[i]) + int(block_sizes[i])
                
                if int(split_key[1]) >= start_pos and int(split_key[1]) <= end_pos:
                    overlap_found = 1
                    break
            
            if overlap_found:
                if site_hash[site]:
                    site_hash[site] += 1
                else:
                    site_hash[site] = 1
        
        if not overlap_found:
            if discard_hash[site]:
                discard_hash[site] += 1
            else:
                discard_hash[site] = 1
    
    return site_hash, discard_hash


def finalize_output(infile, outfile, site_hash, discard_hash, min_mismatch):
    """
    Writes the final filtered output.
    """
    with open(infile, 'r') as sites_file, \
         open(outfile, 'w') as out_file, \
         open(f"{outfile}_failed", 'w') as failed_file:
        
        for line in sites_file:
            line = line.strip()
            if not line:
                continue
            
            fields = line.split()
            name = f"{fields[0]}_{fields[1]}"
            
            if name in site_hash:
                new_alter = site_hash[name]
                
                if name in discard_hash:
                    discard_num = discard_hash[name]
                else:
                    discard_num = 0
                
                # Parse original coverage and alter count
                cov_str = fields[2]
                cov_parts = cov_str.split(',')
                cov = int(cov_parts[0])
                old_alter = int(cov_parts[1])
                
                # Recalculate coverage and editing frequency
                new_cov = cov - (old_alter - new_alter)
                new_edit_freq = f"{new_alter/new_cov:.3f}" if new_cov > 0 else "0.000"
                
                if new_alter >= min_mismatch and new_alter > discard_num:
                    out_file.write(f"{fields[0]}\t{fields[1]}\t{new_cov},{new_alter}\t{fields[3]}\t{fields[4]}\t{new_edit_freq}\n")
                else:
                    reason = f"failed_freq #mismatches: {new_alter} #minimumMismatchesNecessary: {min_mismatch} #discarded mismatch reads: {discard_num}"
                    failed_file.write(f"{line}\t{reason}\n")
            else:
                failed_file.write(f"{line}\tfailed_totalcover\n")


def main():
    parser = argparse.ArgumentParser(
        description="BLAT filter for RNA editing candidates. Python replacement for pblat_candidates_ln.pl.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--infile", required=True, help="File containing a list of variants to be filtered.")
    parser.add_argument("--outfile", required=True, help="Output file for filtered variants.")
    parser.add_argument("--bamfile", required=True, help="BAM file with mapped reads.")
    parser.add_argument("--refgenome", required=True, help="FASTA file of the reference genome.")
    parser.add_argument("--minbasequal", type=int, default=25, help="Minimum base quality of a mismatch (default: 25).")
    parser.add_argument("--minmismatch", type=int, default=1, help="Minimum number of mismatches supported by correctly mapped reads (default: 1).")
    parser.add_argument("--scorelimit", type=float, default=0.95, help="Max score of secondary alignments relative to best (default: 0.95).")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads for pblat (default: 1).")
    parser.add_argument("--pblat-path", default="pblat", help="Path to the pblat executable (default: 'pblat').")
    parser.add_argument("--samtools-path", default="samtools", help="Path to the samtools executable (default: 'samtools').")
    parser.add_argument("--illumina", action="store_true", help="Reads in the bam file are in Illumina 1.3+ FASTQ-like format (quality offset 64).")

    args = parser.parse_args()
    
    # Set quality offset based on format
    qual_offset = 64 if args.illumina else 33
    
    # Define temporary file names
    base_name = args.infile
    fa_path = f"{base_name}.fatmp"
    psl_path = f"{base_name}.psltmp"

    try:
        # Step 1: Create FASTA from BAM for reads supporting variants
        if not prepare_blat_input(args.infile, args.bamfile, fa_path, args.minbasequal, qual_offset, args.threads, args.samtools_path):
            print("No reads found supporting the variants. Exiting.")
            # Create empty output files to prevent pipeline failure
            open(args.outfile, 'w').close()
            open(f"{args.outfile}_failed", 'w').close()
            sys.exit(0)

        # Check if the generated FASTA file is empty (meaning no reads were extracted)
        if os.path.getsize(fa_path) == 0:
            print("No reads extracted for pblat (empty fasta). Exiting.")
            open(args.outfile, 'w').close()
            open(f"{args.outfile}_failed", 'w').close()
            sys.exit(0)

        # Step 2: Run pblat
        if not run_pblat(args.refgenome, fa_path, psl_path, args.threads, args.pblat_path):
            sys.stderr.write("PBLAT failed. Exiting.\n")
            sys.exit(1)

        # Step 3: Process pblat results
        site_hash, discard_hash = process_pblat_results(psl_path, args.scorelimit)
        
        # Step 4: Write final filtered output
        finalize_output(args.infile, args.outfile, site_hash, discard_hash, args.minmismatch)

    finally:
        # Clean up temporary files
        if os.path.exists(fa_path):
            try:
                os.remove(fa_path)
            except OSError:
                pass
        if os.path.exists(psl_path):
            try:
                os.remove(psl_path)
            except OSError:
                pass
    
    print("PBLAT filtering complete.")


if __name__ == "__main__":
    main()
