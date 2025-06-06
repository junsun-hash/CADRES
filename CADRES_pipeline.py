#!/usr/bin/env python
import os
import sys
import argparse
import subprocess
import logging
import re
import json
from glob import glob
from collections import defaultdict

# --- Configuration for Logging ---
# Configure logging to output to stdout
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    stream=sys.stdout)

# --- Helper Functions ---
def run_command(command_parts, shell=False, cwd=None):
    """
    Executes an external command and logs its execution and outcome.
    Args:
        command_parts (list): A list of strings representing the command and its arguments.
        shell (bool): Whether to use the shell for command execution (use with caution).
        cwd (str, optional): The working directory for the command. Defaults to None.
    """
    command_str = ' '.join(command_parts)
    logging.info(f"Executing command: {command_str} (in {cwd or os.getcwd()})")
    try:
        process = subprocess.run(command_parts, check=True, shell=shell, cwd=cwd,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if process.stdout:
            logging.info(f"Command STDOUT:\n{process.stdout.strip()}")
        if process.stderr:
            logging.info(f"Command STDERR:\n{process.stderr.strip()}") # GATK often prints info to stderr
        logging.info(f"Command successful: {command_str}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {command_str}")
        logging.error(f"Return code: {e.returncode}")
        if e.stdout:
            logging.error(f"STDOUT:\n{e.stdout.strip()}")
        if e.stderr:
            logging.error(f"STDERR:\n{e.stderr.strip()}")
        sys.exit(1)

def ensure_dir(path):
    """
    Ensures that a directory exists; if not, it creates it.
    Args:
        path (str): The directory path to check/create.
    """
    if not os.path.exists(path):
        logging.info(f"Creating directory: {path}")
        os.makedirs(path)

def convert_vcf_format(input_vcf, output_file):
    """
    Replicates the functionality of revised_convertVCF.sh to convert VCF to a specific format for SNPiR.
    It processes the first sample in the VCF.
    Args:
        input_vcf (str): Path to the input VCF file.
        output_file (str): Path to the output file in SNPiR format.
    """
    logging.info(f"Converting VCF {input_vcf} to SNPiR format -> {output_file}")
    written_count = 0
    with open(input_vcf, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 10: # Must have CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, and at least one sample
                logging.warning(f"Skipping line due to insufficient fields (expected >=10, got {len(fields)}): {line.strip()}")
                continue
                
            chrom = fields[0]
            pos = fields[1]
            ref_allele = fields[3]
            alt_allele = fields[4]
            
            # Use the first sample's data (fields[9] is 0-indexed for the 10th column)
            sample_data_str = fields[9]
            format_fields = fields[8].split(':')
            
            try:
                ad_index = format_fields.index('AD')
            except ValueError:
                logging.warning(f"AD field not found in FORMAT string for line: {line.strip()}. Skipping.")
                continue
                
            sample_values = sample_data_str.split(':')
            if ad_index >= len(sample_values):
                logging.warning(f"AD index out of bounds for sample data on line: {line.strip()}. Skipping.")
                continue

            allele_depths_str = sample_values[ad_index]
            allele_depths = allele_depths_str.split(',')
            
            if len(allele_depths) < 2:
                logging.warning(f"Allele depths (AD) field does not contain at least two values: '{allele_depths_str}' on line: {line.strip()}. Skipping.")
                continue

            try:
                ref_depth = int(allele_depths[0])
                alt_depth = int(allele_depths[1])
            except ValueError:
                logging.warning(f"Could not parse AD values '{allele_depths_str}' as integers on line: {line.strip()}. Skipping.")
                continue

            total_depth = ref_depth + alt_depth
            depths_field = f"{total_depth},{alt_depth}"
            alt_fraction = "0.01" # Fixed value as in original script
            
            outfile.write(f"{chrom}\t{pos}\t{depths_field}\t{ref_allele}\t{alt_allele}\t{alt_fraction}\n")
            written_count +=1
    logging.info(f"VCF format conversion complete. Wrote {written_count} variants to {output_file}.")


def run_pblat_filter(infile, outfile, bamfile, refgenome, min_base_qual, min_mismatch, score_limit, threads, qual_offset=33):
    """
    Replicates the functionality of pblat_candidates_ln.pl for filtering variants using pblat.
    Args:
        infile (str): Input variant file (SNPiR format).
        outfile (str): Output file for pblat-filtered variants.
        bamfile (str): BAM file containing aligned reads.
        refgenome (str): Reference genome FASTA file.
        min_base_qual (int): Minimum base quality for a mismatch to be considered.
        min_mismatch (int): Minimum number of mismatches supported by correctly mapped reads.
        score_limit (float): Fraction of max read score for duplicate mapping consideration.
        threads (int): Number of threads for samtools and pblat.
        qual_offset (int): Quality score offset (e.g., 33 for Phred+33).
    """
    logging.info("======== Starting pblat filtering process ========")
    
    tmp_dir = os.path.dirname(outfile) or "." # Ensure tmp_dir is valid
    base_outfile_name = os.path.basename(outfile)
    fa_file = os.path.join(tmp_dir, f"{base_outfile_name}.tmp.fa")
    psl_file = os.path.join(tmp_dir, f"{base_outfile_name}.tmp.psl")
    
    # --- Part 1: Prepare pblat input (extract reads from BAM) ---
    logging.info(f"Preparing pblat input: Extracting reads from {bamfile} to {fa_file}")
    read_extraction_count = 0
    with open(fa_file, 'w') as fa_handle, open(infile, 'r') as variant_file:
        for line_idx, line in enumerate(variant_file):
            fields = line.strip().split()
            if len(fields) < 5:
                logging.warning(f"Skipping malformed line in variant file (expected >=5 fields): {line.strip()}")
                continue
            chrom, pos_str, edit_nuc = fields[0], fields[1], fields[4]
            
            try:
                variant_pos = int(pos_str)
            except ValueError:
                logging.warning(f"Skipping variant with non-integer position: {line.strip()}")
                continue

            bam_region = f"{chrom}:{variant_pos}-{variant_pos}"
            mismatch_read_count = 0
            
            cmd_samtools = ['samtools', 'view', '-@', str(threads), bamfile, bam_region]
            try:
                process = subprocess.run(cmd_samtools, capture_output=True, text=True, check=True)
                sam_output = process.stdout
            except subprocess.CalledProcessError as e:
                logging.error(f"Samtools view failed for region {bam_region}: {e.stderr}")
                continue # Skip to next variant if samtools fails for this region
            
            for sam_line in sam_output.strip().split('\n'):
                if not sam_line: continue
                sam_fields = sam_line.split('\t')
                if len(sam_fields) < 11: continue # Basic check for valid SAM line

                read_start_str, cigar, seq, quals = sam_fields[3], sam_fields[5], sam_fields[9], sam_fields[10]
                if cigar == '*': continue # Skip unmapped reads or reads with no CIGAR

                try:
                    read_start = int(read_start_str)
                except ValueError:
                    logging.warning(f"Skipping SAM record with non-integer start: {sam_line}")
                    continue

                current_ref_pos = read_start
                read_seq_pos = 0
                base_at_variant_pos_found = False
                
                cigar_tuples = re.findall(r'(\d+)([MIDNSHPX=])', cigar)
                
                for length_str, op in cigar_tuples:
                    length = int(length_str)
                    if op in 'M=X': # Consumes reference and query
                        for _ in range(length):
                            if current_ref_pos == variant_pos:
                                if read_seq_pos < len(seq) and read_seq_pos < len(quals): # Boundary checks
                                    if seq[read_seq_pos] == edit_nuc and \
                                       (ord(quals[read_seq_pos]) - qual_offset) >= min_base_qual:
                                        base_at_variant_pos_found = True
                                else: # Should not happen with valid CIGAR/SEQ/QUAL
                                     logging.warning(f"Index out of bounds for SEQ/QUAL at pos {read_seq_pos} for read from {bam_region}. Read length: {len(seq)}, Qual length: {len(quals)}")
                                break # Found the variant position
                            current_ref_pos += 1
                            read_seq_pos += 1
                    elif op in 'ISHP': # Consumes query only (I, S, H) or neither (P)
                        if op != 'P': read_seq_pos += length
                    elif op in 'DN': # Consumes reference only
                        current_ref_pos += length
                    if base_at_variant_pos_found: break
                
                if base_at_variant_pos_found:
                    fa_handle.write(f">{chrom}-{variant_pos}-{mismatch_read_count}\n{seq}\n")
                    mismatch_read_count += 1
            read_extraction_count += mismatch_read_count

    if read_extraction_count == 0:
        logging.warning(f"No reads extracted for pblat. Output file {outfile} will be empty or not created if it relies on pblat results.")
        # Create an empty output file to prevent downstream errors if expected
        open(outfile, 'w').close()
        open(f"{outfile}_failed", 'w').close()
        if os.path.exists(fa_file): os.remove(fa_file) # Clean up empty fa_file
        return

    logging.info(f"Extracted {read_extraction_count} reads for pblat.")

    # --- Part 2: Run pblat ---
    logging.info(f"Running pblat with {threads} threads...")
    cmd_pblat = [
        'pblat', f"-threads={threads}", '-stepSize=5', '-repMatch=2253',
        '-minScore=20', '-minIdentity=0', '-noHead', refgenome, fa_file, psl_file
    ]
    run_command(cmd_pblat)

    # --- Part 3: Process pblat results ---
    logging.info(f"Processing pblat results from {psl_file}")
    psl_hash = defaultdict(list)
    with open(psl_file, 'r') as psl_handle:
        for line in psl_handle:
            psl_fields = line.strip().split()
            if len(psl_fields) < 21: continue # Ensure enough fields for parsing
            # PSL format: matches, misMatches, repMatches, nCount, qNumInsert, qBaseInsert, tNumInsert, tBaseInsert, strand, qName, qSize, qStart, qEnd, tName, tSize, tStart, tEnd, blockCount, blockSizes, qStarts, tStarts
            name = psl_fields[9] 
            blat_score_info = [psl_fields[0], psl_fields[13], psl_fields[17], psl_fields[18], psl_fields[20]]
            # score, targetName, numBlocks, blockSizes (comma-sep), targetStarts (comma-sep)
            psl_hash[name].append(blat_score_info)

    site_passed_count = defaultdict(int) # Counts reads passing filter per site
    site_discarded_count = defaultdict(int) # Counts reads discarded per site

    for psl_key, score_entries in psl_hash.items():
        key_parts = psl_key.split('-')
        if len(key_parts) < 2: continue # Should be chr-pos-idx
        chrom_from_key, pos_from_key_str = key_parts[0], key_parts[1]
        
        try:
            pos_from_key = int(pos_from_key_str)
        except ValueError:
            continue

        site_id = f"{chrom_from_key}_{pos_from_key_str}"
        
        # Sort score entries by score (first element of blat_score_info) in descending order
        score_entries.sort(key=lambda x: int(x[0]), reverse=True)
        
        if not score_entries: continue

        best_score_info = score_entries[0]
        best_score_val = int(best_score_info[0])
        best_score_target_name = best_score_info[1] # tName
        
        second_best_score_val = 0
        if len(score_entries) > 1:
            second_best_score_val = int(score_entries[1][0])
            
        alignment_passes_filter = False
        # Check if best alignment is on the correct chromosome and significantly better than second best
        if best_score_target_name == chrom_from_key and \
           second_best_score_val < (best_score_val * score_limit):
            
            num_blocks = int(best_score_info[2])
            block_sizes_str = best_score_info[3].strip(',')
            target_starts_str = best_score_info[4].strip(',')

            if not block_sizes_str or not target_starts_str: continue # Skip if block info is missing

            block_sizes = [int(x) for x in block_sizes_str.split(',')]
            target_starts = [int(x) for x in target_starts_str.split(',')]
            
            # Check if the variant position is covered by any block in the best alignment
            for i in range(num_blocks):
                block_start_on_ref = target_starts[i] + 1 # 1-based
                block_end_on_ref = target_starts[i] + block_sizes[i]
                if block_start_on_ref <= pos_from_key <= block_end_on_ref:
                    alignment_passes_filter = True
                    break
        
        if alignment_passes_filter:
            site_passed_count[site_id] += 1
        else:
            site_discarded_count[site_id] += 1
            
    # --- Part 4: Write final filtered variant sites ---
    logging.info(f"Writing pblat-filtered variants to {outfile} and failed to {outfile}_failed")
    final_written_count = 0
    with open(infile, 'r') as variant_file_again, \
         open(outfile, 'w') as out_handle, \
         open(f"{outfile}_failed", 'w') as failed_handle:
        for line in variant_file_again:
            fields = line.strip().split()
            if len(fields) < 2: continue
            site_id = f"{fields[0]}_{fields[1]}" # chrom_pos
            
            passed_reads_for_site = site_passed_count.get(site_id, 0)
            discarded_reads_for_site = site_discarded_count.get(site_id, 0)
            
            # Original Perl: if ($newalter >= $MINMISMATCH && $newalter > $discardnum)
            # Here, passed_reads_for_site is newalter (number of reads supporting variant after pblat)
            if passed_reads_for_site >= min_mismatch and passed_reads_for_site > discarded_reads_for_site:
                # Outputting only chr, pos, ref, alt for bedtools intersect step later
                # Original input file format for pblat (infile) is: $1\t$2\t$depths_field\t$3\t$4\t$alt_fraction
                # $3 is ref, $4 is alt.
                out_handle.write(f"{fields[0]}\t{fields[1]}\t{fields[3]}\t{fields[4]}\n")
                final_written_count += 1
            else:
                failed_handle.write(f"{line.strip()}\tfailed_pblat_filter (passed_reads:{passed_reads_for_site}, min_mismatch:{min_mismatch}, discarded_reads:{discarded_reads_for_site})\n")
    
    # Clean up temporary files
    if os.path.exists(fa_file): os.remove(fa_file)
    if os.path.exists(psl_file): os.remove(psl_file)
    logging.info(f"pblat filtering process complete. Wrote {final_written_count} variants to {outfile}.")

# --- Main Pipeline Stages ---
def run_boost_stage(config):
    """
    Executes the first part of the CADRES pipeline (BAM preprocessing and initial variant calling).
    Args:
        config (dict): A dictionary containing all necessary configurations.
    """
    logging.info("======== Stage 1: CADRES Boost Preprocessing ========")
    paths = config['paths']
    settings = config['settings']
    samples_config = config['samples']

    script_dir = paths['script_dir']
    output_dir = paths['output_dir']
    
    boost_working_dir = os.path.join(output_dir, "boost_workspace")
    dna_bqsr_dir = os.path.join(output_dir, "dna_bqsr_workspace")
    ensure_dir(boost_working_dir)
    ensure_dir(dna_bqsr_dir)

    # --- 1.1 RNA BAM Preprocessing (SplitNCigarReads via bam_calibration_RNA_boost.py) ---
    logging.info("---- 1.1 Starting RNA BAM preprocessing ----")
    split_bam_files = []
    for rna_sample_name in samples_config['all_rna_samples']:
        output_prefix_boost = os.path.join(boost_working_dir, rna_sample_name)
        # Assuming BAMs are named like {rna_sample_name}.final.bam in rna_bam_dir
        input_rna_bam = os.path.join(paths['rna_bam_dir'], f"{rna_sample_name}.final.bam")
        if not os.path.exists(input_rna_bam):
            logging.error(f"Input RNA BAM not found: {input_rna_bam}. Please check sample names and rna_bam_dir.")
            sys.exit(1)

        cmd_rna_calib_boost = [
            'python', os.path.join(script_dir, 'bam_calibration_RNA_boost.py'),
            '--bam', input_rna_bam,
            '--output', output_prefix_boost, # This script appends _split.bam etc.
            '--genome', paths['genome']
        ]
        run_command(cmd_rna_calib_boost)
        split_bam_files.append(f"{output_prefix_boost}_split.bam") # Path to the output of bam_calibration_RNA_boost.py
    
    # --- 1.2 DNA BAM Calibration (BQSR via bam_calibration_DNA.py) ---
    logging.info("---- 1.2 Starting DNA BAM calibration ----")
    dna_output_prefix_bqsr = os.path.join(dna_bqsr_dir, settings['dna_sample_name'])
    cmd_dna_calib = [
        'python', os.path.join(script_dir, 'bam_calibration_DNA.py'),
        '--bam', paths['dna_bam'],
        '--output', dna_output_prefix_bqsr,
        '--genome', paths['genome'],
        '--known', paths['known_snv_bqsr'] # Use known_snv_bqsr for DNA BQSR
    ]
    run_command(cmd_dna_calib)
    calibrated_dna_bam_path = f"{dna_output_prefix_bqsr}_recalibration.bam"

    # --- 1.3 Initial Variant Calling with Mutect2 ---
    logging.info("---- 1.3 Starting initial variant calling (Mutect2) ----")
    project_prefix = settings['project_prefix']
    mutect_vcf_path = os.path.join(boost_working_dir, f"{project_prefix}.mutect.vcf")
    
    cmd_mutect2_boost = ['gatk', 'Mutect2', '-R', paths['genome']]
    for bam_file in split_bam_files:
        cmd_mutect2_boost.extend(['-I', bam_file])
    cmd_mutect2_boost.extend(['-I', calibrated_dna_bam_path])
    cmd_mutect2_boost.extend(['-normal', settings['dna_sample_name']])
    cmd_mutect2_boost.extend(['-O', mutect_vcf_path])
    run_command(cmd_mutect2_boost, cwd=boost_working_dir)

    # --- 1.4 Filter Mutect2 Calls ---
    logging.info("---- 1.4 Filtering Mutect2 calls ----")
    filtered_vcf_boost = os.path.join(boost_working_dir, f"{project_prefix}_boost.filter.vcf")
    cmd_filter_boost = [
        'gatk', 'FilterMutectCalls',
        '-V', mutect_vcf_path,
        '-R', paths['genome'],
        '-O', filtered_vcf_boost,
        '--max-events-in-region', '4' # As per original script
    ]
    run_command(cmd_filter_boost, cwd=boost_working_dir)

    # --- 1.5 Create 'known.vcf' for DVR stage (PASS variants) ---
    pass_vcf_boost = os.path.join(boost_working_dir, f"{project_prefix}.known.vcf") # This will be input for DVR BQSR
    cmd_bcftools_pass = ['bcftools', 'view', '-f', 'PASS', '-o', pass_vcf_boost, filtered_vcf_boost]
    run_command(cmd_bcftools_pass, cwd=boost_working_dir)
    logging.info(f"Generated known sites for DVR stage: {pass_vcf_boost}")

    # --- 1.6 Index the 'known.vcf' ---
    logging.info("---- 1.6 Indexing VCF file for DVR stage ----")
    cmd_index_known_vcf = ['gatk', 'IndexFeatureFile', '-I', pass_vcf_boost]
    run_command(cmd_index_known_vcf, cwd=boost_working_dir)

    logging.info("======== Stage 1: CADRES Boost Preprocessing Complete ========")
    # Return paths needed for DVR stage
    return {
        "boost_known_vcf": pass_vcf_boost,
        "calibrated_dna_bam": calibrated_dna_bam_path,
        "split_rna_bams": split_bam_files,
        "boost_working_dir": boost_working_dir,
        "dna_bqsr_dir": dna_bqsr_dir
    }

def run_dvr_stage(config, boost_outputs):
    """
    Executes the second part of the CADRES pipeline (Differential Variant Region analysis).
    Args:
        config (dict): A dictionary containing all necessary configurations.
        boost_outputs (dict): Outputs from the run_boost_stage function.
    """
    logging.info("======== Stage 2: CADRES DVR Analysis ========")
    paths = config['paths']
    settings = config['settings']
    samples_config = config['samples']
    
    script_dir = paths['script_dir']
    output_dir = paths['output_dir']
    project_prefix = settings['project_prefix']

    rna_bqsr_dir = os.path.join(output_dir, "rna_bqsr_workspace")
    rv_dvr_working_dir = os.path.join(output_dir, "rv_dvr_workspace")
    ensure_dir(rna_bqsr_dir)
    ensure_dir(rv_dvr_working_dir)

    known_vcf_from_boost = boost_outputs['boost_known_vcf']
    
    # --- 2.1 RNA BAM BQSR ---
    logging.info("---- 2.1 Starting RNA BQSR ----")
    bqsr_reports = []
    calibrated_rna_bams_dvr = [] # Store paths to BQSR'd RNA BAMs for this stage

    for i, rna_sample_name in enumerate(samples_config['all_rna_samples']):
        # Input for BQSR is the split BAM from boost stage
        input_split_bam_for_bqsr = boost_outputs['split_rna_bams'][i] 
        
        report_file = os.path.join(rna_bqsr_dir, f"{rna_sample_name}_recalibration_report.grp")
        cmd_bqsr_recal = [
            'gatk', 'BaseRecalibrator', '-I', input_split_bam_for_bqsr, 
            '-R', paths['genome'],
            '-O', report_file, '--known-sites', known_vcf_from_boost
        ]
        run_command(cmd_bqsr_recal, cwd=rna_bqsr_dir)
        bqsr_reports.append(f"-I {report_file}") # For GatherBQSRReports

        # Apply BQSR (this step was separated in original shell script, doing it per sample)
        # However, original GatherBQSRReports implies a single merged report is used for ApplyBQSR.
        # Let's follow the Gather -> Apply pattern

    # Gather BQSR reports
    gathered_report_path = os.path.join(rna_bqsr_dir, f"{project_prefix}_recalibration_report.grp")
    # shell=True needed if bqsr_reports is a list of strings like '-I file1 -I file2' to be expanded by shell
    # If run_command passes it as a list of arguments, shell=False is fine
    cmd_gather_bqsr = ['gatk', 'GatherBQSRReports'] + [item for sublist in [r.split() for r in bqsr_reports] for item in sublist] + ['-O', gathered_report_path]
    run_command(cmd_gather_bqsr, cwd=rna_bqsr_dir) # shell=True if complex string args

    # Apply BQSR using the gathered report
    for i, rna_sample_name in enumerate(samples_config['all_rna_samples']):
        input_split_bam_for_applybqsr = boost_outputs['split_rna_bams'][i]
        output_calibrated_rna_bam = os.path.join(rna_bqsr_dir, f"{rna_sample_name}_recalibration.bam")
        cmd_apply_bqsr = [
            'gatk', 'ApplyBQSR', '-R', paths['genome'], 
            '-I', input_split_bam_for_applybqsr,
            '--bqsr-recal-file', gathered_report_path, 
            '-O', output_calibrated_rna_bam
        ]
        run_command(cmd_apply_bqsr, cwd=rna_bqsr_dir)
        calibrated_rna_bams_dvr.append(output_calibrated_rna_bam)
    
    # --- 2.2 Calculate Contamination ---
    logging.info("---- 2.2 Calculating contamination ----")
    contamination_table_args = []
    for i, rna_sample_name in enumerate(samples_config['all_rna_samples']):
        current_calibrated_rna_bam = calibrated_rna_bams_dvr[i]
        pileups_table = os.path.join(rna_bqsr_dir, f"{rna_sample_name}_pileups.table")
        contamination_table = os.path.join(rna_bqsr_dir, f"{rna_sample_name}_contamination.table")
        
        cmd_getpileups = [
            'gatk', 'GetPileupSummaries', '-I', current_calibrated_rna_bam,
            '-V', paths['genome_ad'], '-L', paths['genome_ad'], '-O', pileups_table
        ]
        run_command(cmd_getpileups, cwd=rna_bqsr_dir)
        
        cmd_calccontam = [
            'gatk', 'CalculateContamination', '-I', pileups_table, '-O', contamination_table
        ]
        run_command(cmd_calccontam, cwd=rna_bqsr_dir)
        contamination_table_args.append(f"--contamination-table {contamination_table}")

    # --- 2.3 RV & DVR Variant Calling (Mutect2) ---
    logging.info("---- 2.3 Starting RV & DVR variant calling ----")
    # Ensure we are in the rv_dvr_working_dir for subsequent relative paths
    # os.chdir(rv_dvr_working_dir) # Be cautious with os.chdir if not strictly necessary
    
    mutect_temp_vcf = os.path.join(rv_dvr_working_dir, f"{project_prefix}.temp.vcf")
    cmd_mutect2_dvr = ['gatk', 'Mutect2', '-R', paths['genome']]
    for bam in calibrated_rna_bams_dvr: # Use BQSR'd RNA BAMs from this DVR stage
        cmd_mutect2_dvr.extend(['-I', bam])
    cmd_mutect2_dvr.extend(['-I', boost_outputs['calibrated_dna_bam']]) # Calibrated DNA BAM from boost stage
    cmd_mutect2_dvr.extend(['-normal', settings['dna_sample_name']])
    cmd_mutect2_dvr.extend(['--germline-resource', paths['genome_ad']])
    cmd_mutect2_dvr.extend(['-O', mutect_temp_vcf])
    run_command(cmd_mutect2_dvr, cwd=rv_dvr_working_dir)

    # Filter Mutect Calls
    mutect_temp2_vcf = os.path.join(rv_dvr_working_dir, f"{project_prefix}.temp2.vcf")
    cmd_filter_dvr = ['gatk', 'FilterMutectCalls', '-V', mutect_temp_vcf, '-R', paths['genome']] \
                   + [item for sublist in [ct.split() for ct in contamination_table_args] for item in sublist] \
                   + ['--min-median-base-quality', '12', '--max-events-in-region', '4', '-O', mutect_temp2_vcf]
    run_command(cmd_filter_dvr, cwd=rv_dvr_working_dir)
    
    # bcftools view for PASS variants
    final_raw_vcf = os.path.join(rv_dvr_working_dir, f"{project_prefix}.vcf")
    cmd_bcftools_dvr = ['bcftools', 'view', '-f', 'PASS', '-o', final_raw_vcf, mutect_temp2_vcf]
    run_command(cmd_bcftools_dvr, cwd=rv_dvr_working_dir)

    # Select SNVs
    snv_vcf_path = os.path.join(rv_dvr_working_dir, f"{project_prefix}.SNV.vcf")
    cmd_select_snv = [
        'gatk', 'SelectVariants', '-R', paths['genome'], '-V', final_raw_vcf,
        '--select-type-to-include', 'SNP', '-O', snv_vcf_path
    ]
    run_command(cmd_select_snv, cwd=rv_dvr_working_dir)

    # --- 2.4 VCF Conversion, Homopolymer, and pBLAT filtering ---
    logging.info("---- 2.4 Starting VCF post-processing and pBLAT filtering ----")
    snv_for_snipir_path = os.path.join(rv_dvr_working_dir, f"{project_prefix}_SNPIR.SNV.vcf")
    convert_vcf_format(snv_vcf_path, snv_for_snipir_path) # Uses first sample in VCF

    # Filter homopolymers (external Perl script)
    homopolymer_filtered_vcf = os.path.join(rv_dvr_working_dir, f"{project_prefix}.homo.vcf")
    cmd_filter_homo = [
        'perl', os.path.join(script_dir, 'filter_homopolymer_nucleotides.pl'), # Ensure this script is in script_dir
        '-infile', snv_for_snipir_path,
        '-outfile', homopolymer_filtered_vcf,
        '-refgenome', paths['genome']
    ]
    # Ensure PERL5LIB is set if this script has dependencies, or modify script to use FindBin
    logging.info("Attempting to run filter_homopolymer_nucleotides.pl. Ensure Perl is installed and script is executable.")
    run_command(cmd_filter_homo, cwd=rv_dvr_working_dir)

    # Merge RNA BAMs for pBLAT
    merged_rna_bam_pblat = os.path.join(rv_dvr_working_dir, f"{project_prefix}_all_rna_for_pblat.bam")
    sorted_rna_bam_pblat = os.path.join(rv_dvr_working_dir, f"{project_prefix}_all_rna_for_pblat_sorted.bam")
    
    cmd_merge_bams = ['samtools', 'merge', '-f', '-@', str(settings['threads']), merged_rna_bam_pblat] + calibrated_rna_bams_dvr
    run_command(cmd_merge_bams)
    cmd_sort_merged = ['samtools', 'sort', '-@', str(settings['threads']), '-o', sorted_rna_bam_pblat, merged_rna_bam_pblat]
    run_command(cmd_sort_merged)
    cmd_index_merged = ['samtools', 'index', '-@', str(settings['threads']), sorted_rna_bam_pblat]
    run_command(cmd_index_merged)

    # Run pBLAT filtering
    pblat_filtered_variants = os.path.join(rv_dvr_working_dir, f"{project_prefix}.pblat.vcf")
    run_pblat_filter(
        infile=homopolymer_filtered_vcf, # Input from homopolymer filtering
        outfile=pblat_filtered_variants,
        bamfile=sorted_rna_bam_pblat,
        refgenome=paths['genome'],
        min_base_qual=settings.get('min_base_qual_pblat', 25),
        min_mismatch=settings.get('min_mismatch_pblat', 1),
        score_limit=settings.get('score_limit_pblat', 0.95),
        threads=settings['threads']
    )
    # Clean up merged BAMs for pBLAT
    os.remove(merged_rna_bam_pblat)
    os.remove(sorted_rna_bam_pblat)
    os.remove(f"{sorted_rna_bam_pblat}.bai")

    # Create BED file from pblat output for bedtools
    pblat_output_bed = os.path.join(rv_dvr_working_dir, f"{project_prefix}.pblat_passed.bed")
    with open(pblat_filtered_variants, 'r') as pblat_in, open(pblat_output_bed, 'w') as bed_out:
        for line in pblat_in: # pblat_filtered_variants contains: chr, pos, ref, alt
            f = line.strip().split()
            if len(f) == 4: # Expecting chr, pos, ref, alt
                 # BED is 0-indexed start, 1-indexed end
                bed_out.write(f"{f[0]}\t{int(f[1])-1}\t{f[1]}\t{f[2]}_{f[3]}\n") # name field with ref_alt
    
    # Intersect with SNV VCF
    final_vcf_for_mats = os.path.join(rv_dvr_working_dir, f"{project_prefix}.final.vcf")
    # Use original SNV.vcf as -a, and pblat_output_bed as -b
    cmd_bed_intersect = [
        'bedtools', 'intersect', '-a', snv_vcf_path, '-b', pblat_output_bed, 
        '-wa', '-header' 
    ]
    with open(final_vcf_for_mats, "w") as f_intersect_out:
        run_command(cmd_bed_intersect, cwd=rv_dvr_working_dir) # Output directly to file

    # --- 2.5 Differential Analysis (mpileup, MATS_LRT, FDR, Annotation) ---
    logging.info("---- 2.5 Starting differential editing analysis ----")
    
    # Prepare for mpileup: separate first and second reads (as per original script)
    rna_bam_sep_list_mpileup = []
    temp_sep_bams_for_mpileup = []
    for calib_bam_path in calibrated_rna_bams_dvr:
        base_name = os.path.splitext(os.path.basename(calib_bam_path))[0] # e.g. RNA_Sample1_Control_recalibration
        first_end_bam = os.path.join(rv_dvr_working_dir, f"{base_name}_sep_firstEnd.bam")
        second_end_bam = os.path.join(rv_dvr_working_dir, f"{base_name}_sep_secondEnd.bam")
        
        run_command(['samtools', 'view', '-@', str(settings['threads']), '-f', '64', '-b', calib_bam_path, '-o', first_end_bam])
        run_command(['samtools', 'view', '-@', str(settings['threads']), '-f', '128', '-b', calib_bam_path, '-o', second_end_bam])
        
        rna_bam_sep_list_mpileup.extend([first_end_bam, second_end_bam])
        temp_sep_bams_for_mpileup.extend([first_end_bam, second_end_bam])

    # Samtools mpileup
    mpileup_output_file = os.path.join(rv_dvr_working_dir, f"{project_prefix}.pileup")
    cmd_mpileup = [
        'samtools', 'mpileup', '-B', '-d', '100000', '-f', paths['genome'],
        '-l', final_vcf_for_mats, '-q', '30', '-Q', '17', '-a', 
        '-o', mpileup_output_file
    ] + rna_bam_sep_list_mpileup
    run_command(cmd_mpileup, cwd=rv_dvr_working_dir)

    # Clean up temporary separated BAMs for mpileup
    for temp_bam in temp_sep_bams_for_mpileup:
        if os.path.exists(temp_bam): os.remove(temp_bam)

    # vcf_to_mats_input_For_Mutect2.py
    inc_txt_file = os.path.join(rv_dvr_working_dir, f"{project_prefix}.inc.txt")
    # Get comma-separated BAM lists for group1 and group2 (full calibrated BAMs)
    group1_bams_str = ",".join([os.path.join(rna_bqsr_dir, f"{s}_recalibration.bam") for s in samples_config['group1_rna_samples']])
    group2_bams_str = ",".join([os.path.join(rna_bqsr_dir, f"{s}_recalibration.bam") for s in samples_config['group2_rna_samples']])
    
    vcf_to_mats_params = settings['vcf_to_mats']
    cmd_vcf2mats = [
        'python', os.path.join(script_dir, 'vcf_to_mats_input_For_Mutect2.py'),
        final_vcf_for_mats, inc_txt_file,
        group1_bams_str, group2_bams_str,
        str(vcf_to_mats_params['min_qual']), str(vcf_to_mats_params['min_depth']),
        vcf_to_mats_params['paired_reads'], mpileup_output_file, vcf_to_mats_params['stranded_reads']
    ]
    run_command(cmd_vcf2mats, cwd=rv_dvr_working_dir)

    # MATS_LRT.py
    mats_lrt_results_prefix = os.path.join(rv_dvr_working_dir, f"{project_prefix}_rMATS-DVR_results")
    mats_lrt_params = settings['mats_lrt']
    cmd_mats_lrt = [
        'python', os.path.join(script_dir, 'MATS_LRT.py'),
        inc_txt_file, mats_lrt_results_prefix,
        str(mats_lrt_params['multiprocessor']), mats_lrt_params['cutoff']
    ]
    run_command(cmd_mats_lrt, cwd=rv_dvr_working_dir)

    # FDR.py
    fdr_input_file = f"{mats_lrt_results_prefix}_rMATS_Result_P.txt"
    fdr_output_file = f"{mats_lrt_results_prefix}_rMATS_Result_FDR.txt"
    cmd_fdr = ['python', os.path.join(script_dir, 'FDR.py'), fdr_input_file, fdr_output_file]
    run_command(cmd_fdr, cwd=rv_dvr_working_dir)

    # snv_annotation.py
    final_annotated_output = os.path.join(rv_dvr_working_dir, f"{project_prefix}_rMATS-DVR_Result.txt")
    summary_stats_output = os.path.join(rv_dvr_working_dir, f"{project_prefix}_rMATS-DVR_Result_summary.txt")
    cmd_snv_anno = [
        'python', os.path.join(script_dir, 'snv_annotation.py'),
        '--input', fdr_output_file,
        '--output', final_annotated_output,
        '--summary', summary_stats_output,
        '--label1', samples_config['group1_label'],
        '--label2', samples_config['group2_label'],
        '--snp', paths['known_snv_annotation'], # Use specific SNV for annotation
        '--repeat', "NA", # As per original script
        '--editing', paths['known_editing_sites'],
        '--gene', paths['gene_anno']
    ]
    run_command(cmd_snv_anno, cwd=rv_dvr_working_dir)
    
    logging.info("======== Stage 2: CADRES DVR Analysis Complete ========")


def main():
    parser = argparse.ArgumentParser(description="CADRES Unified Pipeline for Differential RNA Editing Analysis.")
    parser.add_argument('--config', required=True, help='Path to the JSON configuration file.')
    args = parser.parse_args()

    # Load configuration from JSON file
    try:
        with open(args.config, 'r') as f:
            config = json.load(f)
    except FileNotFoundError:
        logging.error(f"Configuration file not found: {args.config}")
        sys.exit(1)
    except json.JSONDecodeError:
        logging.error(f"Error decoding JSON from configuration file: {args.config}")
        sys.exit(1)

    logging.info(f"Successfully loaded configuration from {args.config}")

    # Create main output directory if it doesn't exist
    ensure_dir(config['paths']['output_dir'])

    # Run Stage 1: Boost
    boost_stage_outputs = run_boost_stage(config)

    # Run Stage 2: DVR
    run_dvr_stage(config, boost_stage_outputs)

    logging.info("======== CADRES Unified Pipeline Successfully Completed ========")

if __name__ == '__main__':
    main()
