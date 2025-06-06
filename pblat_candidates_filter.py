#!/usr/bin/env python
import argparse
import sys
import os
import subprocess
from collections import defaultdict
import pysam

def run_command(command):
    """Runs a command using subprocess and handles errors."""
    try:
        subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"Error executing command: {e.cmd}\n")
        sys.stderr.write(f"Stderr: {e.stderr.decode()}\n")
        sys.exit(1)

def prepare_blat_input(infile, bam_path, fa_path, min_base_qual, qual_offset):
    """
    Reads a variant file, extracts reads from a BAM file supporting the variant,
    and writes them to a FASTA file for PBLAT.
    """
    print("Preparing BLAT input...")
    mismatch_read_count = 0
    try:
        bamfile = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        sys.stderr.write(f"Error opening BAM file {bam_path}: {e}\n")
        sys.exit(1)

    with open(infile, 'r') as sites_file, open(fa_path, 'w') as fa_file:
        for line in sites_file:
            fields = line.strip().split()
            chrom, pos, alt_nuc = fields[0], int(fields[1]), fields[4]
            pos_0based = pos - 1

            # Fetch reads covering the specific site
            try:
                for read in bamfile.fetch(chrom, pos_0based, pos_0based + 1):
                    if read.is_unmapped or read.is_duplicate:
                        continue
                    
                    # Find the read base corresponding to the variant position
                    read_pos = -1
                    ref_pos_lookup = dict(read.get_aligned_pairs(matches_only=True))
                    
                    if pos_0based in ref_pos_lookup:
                        read_pos = ref_pos_lookup[pos_0based]
                        
                        read_base = read.query_sequence[read_pos]
                        base_qual = read.query_qualities[read_pos]

                        if read_base.upper() == alt_nuc.upper() and base_qual >= min_base_qual:
                            # Write read to FASTA file
                            fa_file.write(f">{chrom}-{pos}-{mismatch_read_count}\n{read.query_sequence}\n")
                            mismatch_read_count += 1
            except ValueError as e:
                 sys.stderr.write(f"Warning: Could not fetch reads for region {chrom}:{pos}. Skipping. Pysam error: {e}\n")


    bamfile.close()
    return mismatch_read_count > 0

def run_pblat(refgenome, fa_path, psl_path, threads, pblat_exe):
    """Runs PBLAT on the generated FASTA file."""
    print(f"Blatting reads (threads: {threads})...")
    pblat_params = (
        f"-threads={threads} -stepSize=5 -repMatch=2253 -minScore=20 "
        f"-minIdentity=0 -noHead {refgenome} {fa_path} {psl_path}"
    )
    command = f"{pblat_exe} {pblat_params}"
    run_command(command)

def process_pblat_results(psl_path, score_limit):
    """
    Processes the PSL output to identify uniquely mapped reads.
    Returns two dictionaries: one for valid sites and one for discarded reads.
    """
    print("Processing PBLAT results...")
    psl_hash = defaultdict(list)
    with open(psl_path, 'r') as psl_file:
        for line in psl_file:
            fields = line.strip().split()
            match, t_name, t_start, q_name = fields[0], fields[13], int(fields[15]), fields[9]
            num_blocks, block_sizes, t_starts = fields[17], fields[18], fields[20]
            
            psl_hash[q_name].append({
                'score': int(match),
                't_name': t_name,
                't_start': t_start,
                'num_blocks': int(num_blocks),
                'block_sizes': [int(x) for x in block_sizes.split(',') if x],
                't_starts': [int(x) for x in t_starts.split(',') if x]
            })

    site_hash = defaultdict(int)
    discard_hash = defaultdict(int)

    for read_name, alignments in psl_hash.items():
        chrom, pos_str, _ = read_name.split('-')
        site_key = f"{chrom}_{pos_str}"
        pos = int(pos_str)

        if not alignments:
            continue

        # Sort by score descending
        alignments.sort(key=lambda x: x['score'], reverse=True)
        best_alignment = alignments[0]
        score_best = best_alignment['score']
        
        # Check if secondary alignments are too good
        score_second_best = alignments[1]['score'] if len(alignments) > 1 else 0

        # Check if the best alignment is the correct one and unique enough
        if best_alignment['t_name'] == chrom and score_second_best < (score_best * score_limit):
            # Check if the variant position is within any of the aligned blocks
            overlap_found = False
            for i in range(best_alignment['num_blocks']):
                block_start = best_alignment['t_starts'][i] + 1  # PSL is 0-based, VCF is 1-based
                block_end = block_start + best_alignment['block_sizes'][i] - 1
                if block_start <= pos <= block_end:
                    overlap_found = True
                    break
            
            if overlap_found:
                site_hash[site_key] += 1
            else:
                discard_hash[site_key] += 1
        else:
            discard_hash[site_key] += 1
            
    return site_hash, discard_hash

def finalize_output(infile, outfile, site_hash, discard_hash, min_mismatch):
    """Writes the final filtered output."""
    print("Writing final output...")
    with open(infile, 'r') as sites_file, \
         open(outfile, 'w') as out_file, \
         open(f"{outfile}_failed", 'w') as failed_file:
        for line in sites_file:
            line = line.strip()
            fields = line.split()
            name = f"{fields[0]}_{fields[1]}"
            
            if name in site_hash:
                new_alter = site_hash[name]
                discard_num = discard_hash.get(name, 0)

                if new_alter >= min_mismatch and new_alter > discard_num:
                    # Original script recalculates coverage, but the logic is complex
                    # and depends on fields not clearly defined.
                    # Replicating the simple filtering logic first.
                    # $fields[2] = "new_cov,new_alter"
                    # $fields[5] = "new_edit_freq"
                    # For now, let's just output the passed variants.
                    # The original script does complex recalculation of coverage and frequency.
                    # This is a simplified port focusing on the filtering logic.
                    # The original fields: $1, $2, $3(cov,alter), $4, $5, $6(freq)
                    old_cov, old_alter = map(int, fields[2].split(','))
                    new_cov = old_cov - (old_alter - new_alter)
                    new_freq = f"{new_alter/new_cov:.3f}" if new_cov > 0 else "0.000"
                    
                    out_fields = [
                        fields[0], fields[1], f"{new_cov},{new_alter}",
                        fields[3], fields[4], new_freq
                    ]
                    out_file.write("\t".join(out_fields) + '\n')
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

    args = parser.parse_args()
    
    # Define temporary file names
    base_name = os.path.splitext(args.infile)[0]
    fa_path = f"{base_name}.pblat_temp.fa"
    psl_path = f"{base_name}.pblat_temp.psl"

    try:
        # Step 1: Create FASTA from BAM for reads supporting variants
        if not prepare_blat_input(args.infile, args.bamfile, fa_path, args.minbasequal, 33):
            print("No reads found supporting the variants. Exiting.")
            # Create empty output files to prevent pipeline failure
            open(args.outfile, 'w').close()
            open(f"{args.outfile}_failed", 'w').close()
            sys.exit(0)

        # Step 2: Run pblat
        run_pblat(args.refgenome, fa_path, psl_path, args.threads, args.pblat_path)

        # Step 3: Process pblat results
        site_hash, discard_hash = process_pblat_results(psl_path, args.scorelimit)
        
        # Step 4: Write final filtered output
        finalize_output(args.infile, args.outfile, site_hash, discard_hash, args.minmismatch)

    finally:
        # Clean up temporary files
        if os.path.exists(fa_path):
            os.remove(fa_path)
        if os.path.exists(psl_path):
            os.remove(psl_path)
    
    print("PBLAT filtering complete.")


if __name__ == "__main__":
    main()