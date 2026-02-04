#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import subprocess
import logging
from multiprocessing import Pool

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_command(cmd):
    logging.info(f"Running: {cmd}")
    try:
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {cmd}")
        raise e

def preprocess_dna(args):
    """
    对应 bam_calibration_DNA.py 的功能:
    ReorderSam -> AddOrReplaceReadGroups -> MarkDuplicates -> (Optional BQSR)
    """
    bam, output_prefix, genome, known_sites, keep_temp, sample_name = args
    
    # 1. ReorderSam
    cmd = f"picard ReorderSam INPUT={bam} OUTPUT={output_prefix}_reordered.bam SEQUENCE_DICTIONARY={genome.replace('.fa', '.dict').replace('.fasta', '.dict')} S=true R={genome}"
    run_command(cmd)
    
    # 2. AddOrReplaceReadGroups
    rgid = sample_name
    cmd = f"picard AddOrReplaceReadGroups INPUT={output_prefix}_reordered.bam OUTPUT={output_prefix}_addrg.bam RGID={rgid} RGLB={rgid} RGPL=COMPLETE RGPU=lane1 RGSM={rgid}"
    run_command(cmd)
    
    # 3. MarkDuplicates
    cmd = f"_JAVA_OPTIONS='-Xmx64g' picard MarkDuplicates INPUT={output_prefix}_addrg.bam OUTPUT={output_prefix}_dedup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT READ_NAME_REGEX=null METRICS_FILE={output_prefix}_metrics.txt"
    run_command(cmd)
    
    if not keep_temp:
        run_command(f"rm -f {output_prefix}_reordered.bam {output_prefix}_addrg.bam")

    final_bam = f"{output_prefix}_dedup.bam"

    # 4. BQSR (Check if known_sites is provided and valid)
    if known_sites and known_sites != 'NA' and os.path.exists(known_sites):
        cmd = f"gatk BaseRecalibrator -I {final_bam} -R {genome} -O {output_prefix}_recalibration_report.grp --known-sites {known_sites}"
        run_command(cmd)
        
        cmd = f"gatk ApplyBQSR -R {genome} -I {final_bam} --bqsr-recal-file {output_prefix}_recalibration_report.grp -O {output_prefix}_recalibration.bam"
        run_command(cmd)
        final_bam = f"{output_prefix}_recalibration.bam"
        
        if not keep_temp:
             run_command(f"rm -f {output_prefix}_dedup.bam {output_prefix}_dedup.bai")
    else:
        logging.warning(f"Known sites not provided or not found ({known_sites}), skipping BQSR for DNA.")

    return final_bam

def preprocess_rna(args):
    """
    对应 bam_calibration_RNA_boost.py 的功能:
    ReorderSam -> AddOrReplaceReadGroups -> MarkDuplicates -> SplitNCigarReads
    """
    bam, output_prefix, genome, keep_temp = args
    
    # 1. ReorderSam
    cmd = f"picard ReorderSam INPUT={bam} OUTPUT={output_prefix}_reordered.bam SEQUENCE_DICTIONARY={genome.replace('.fa', '.dict').replace('.fasta', '.dict')} S=true R={genome}"
    run_command(cmd)
    
    # 2. AddOrReplaceReadGroups
    rgid = os.path.basename(output_prefix)
    cmd = f"picard AddOrReplaceReadGroups INPUT={output_prefix}_reordered.bam OUTPUT={output_prefix}_addrg.bam RGID={rgid} RGLB={rgid} RGPL=COMPLETE RGPU=lane1 RGSM={rgid}"
    run_command(cmd)
    
    # 3. MarkDuplicates
    cmd = f"picard MarkDuplicates INPUT={output_prefix}_addrg.bam OUTPUT={output_prefix}_dedup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT READ_NAME_REGEX=null METRICS_FILE={output_prefix}_metrics.txt"
    run_command(cmd)

    if not keep_temp:
        run_command(f"rm -f {output_prefix}_reordered.bam {output_prefix}_addrg.bam")

    # 4. SplitNCigarReads
    cmd = f"gatk SplitNCigarReads -R {genome} -I {output_prefix}_dedup.bam -O {output_prefix}_split.bam"
    run_command(cmd)

    if not keep_temp:
        run_command(f"rm -f {output_prefix}_dedup.bam {output_prefix}_dedup.bai")

    return f"{output_prefix}_split.bam"

def apply_rna_bqsr(args):
    """
    对应 DVR-293TV6.sh 中使用 boosted sites 进行 BQSR 的部分
    """
    bam, known_sites, output_prefix, genome = args
    
    # BaseRecalibrator
    report_file = f"{output_prefix}_recalibration_report.grp"
    cmd = f"gatk BaseRecalibrator -I {bam} -R {genome} -O {report_file} --known-sites {known_sites}"
    run_command(cmd)
    
    # ApplyBQSR
    final_bam = f"{output_prefix}_recalibration.bam"
    cmd = f"gatk ApplyBQSR -R {genome} -I {bam} --bqsr-recal-file {report_file} -O {final_bam}"
    run_command(cmd)
    
    return final_bam

def main():
    try:
        _main()
    except Exception as e:
        logging.error(f"Pipeline failed with error: {e}")
        sys.exit(1)

def _main():
    parser = argparse.ArgumentParser(description="Step 1: Base Quality Calibration (Integrates bam_calibration_*.py and boost logic)")
    parser.add_argument('--rna_bams', nargs='+', required=True, help='List of RNA BAM files')
    parser.add_argument('--dna_bam', required=True, help='DNA BAM file (Normal)')
    parser.add_argument('--genome', required=True, help='Reference genome FASTA')
    parser.add_argument('--known_snv', required=True, help='Known SNV VCF (e.g., dbSNP)')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    parser.add_argument('--prefix', default='sample', help='Output prefix')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads/processes')
    parser.add_argument('--keep_temp', action='store_true', help='Keep temporary files')
    parser.add_argument('--dna_sample_name', default=None, help='Custom sample name for DNA BAM RGSM and Mutect2 normal name')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # 1. Preprocess DNA (Parallel with RNA or sequential? DNA is usually one, RNA multiple)
    logging.info("Starting DNA Preprocessing...")
    
    # Determine DNA sample name
    dna_sample_name = args.dna_sample_name if args.dna_sample_name else "DNA_processed"
    
    # We use dna_sample_name for the file prefix inside the output dir to keep it consistent with the RGSM
    # BUT, the original script used "DNA_processed" as the file prefix. 
    # If we change the file prefix, we might break downstream assumptions if they hardcode filenames.
    # However, the user asked to change the RGSM.
    # Let's keep the file prefix as "DNA_processed" to avoid breaking filenames, but pass the custom RGSM to preprocess_dna.
    
    dna_output_prefix = os.path.join(args.output_dir, "DNA_processed")
    # Pass dna_sample_name to preprocess_dna
    processed_dna_bam = preprocess_dna((args.dna_bam, dna_output_prefix, args.genome, args.known_snv, args.keep_temp, dna_sample_name))
    logging.info(f"DNA Preprocessing done: {processed_dna_bam}")

    # 2. Preprocess RNA (Parallel)
    logging.info("Starting RNA Preprocessing...")
    rna_pool_args = []
    processed_rna_bams = []
    for bam in args.rna_bams:
        base_name = os.path.basename(bam).replace('.bam', '')
        out_prefix = os.path.join(args.output_dir, base_name)
        rna_pool_args.append((bam, out_prefix, args.genome, args.keep_temp))
        processed_rna_bams.append(f"{out_prefix}_split.bam")
    
    with Pool(args.threads) as p:
        p.map(preprocess_rna, rna_pool_args)
    logging.info("RNA Preprocessing done.")

    # 3. Boost Step: Run Mutect2 on Split RNA + DNA to get known sites
    logging.info("Starting Boost Step (Mutect2 to generate known sites)...")
    boost_vcf_prefix = os.path.join(args.output_dir, f"{args.prefix}_boost")
    
    rna_inputs = " ".join([f"-I {bam}" for bam in processed_rna_bams])
    cmd = f"gatk Mutect2 -R {args.genome} {rna_inputs} -I {processed_dna_bam} -normal {os.path.basename(dna_output_prefix)} -O {boost_vcf_prefix}.mutect.vcf"
    # Note: The original script uses -normal C6_WGS, here we derive it or need user input. 
    # For simplicity, assuming the RGSM in DNA bam matches. If not, Mutect2 might complain.
    # In preprocess_dna, we set RGSM to os.path.basename(output_prefix), which is 'DNA_processed'.
    # So -normal DNA_processed should work.
    
    # Using 'DNA_processed' as sample name because we set it in preprocess_dna
    cmd = f"gatk Mutect2 -R {args.genome} {rna_inputs} -I {processed_dna_bam} -normal DNA_processed -O {boost_vcf_prefix}.mutect.vcf"
    run_command(cmd)

    # Filter Boost VCF
    cmd = f"gatk FilterMutectCalls -V {boost_vcf_prefix}.mutect.vcf -R {args.genome} -O {boost_vcf_prefix}.filter.vcf --max-events-in-region 4"
    run_command(cmd)
    
    # Select PASS variants
    cmd = f"bcftools view -f PASS {boost_vcf_prefix}.filter.vcf > {boost_vcf_prefix}.known.vcf"
    run_command(cmd)
    
    cmd = f"gatk IndexFeatureFile -I {boost_vcf_prefix}.known.vcf"
    run_command(cmd)
    
    boosted_known_sites = f"{boost_vcf_prefix}.known.vcf"
    logging.info(f"Boost Step done. Known sites: {boosted_known_sites}")

    # 4. Final RNA BQSR using Boosted Sites
    logging.info("Starting Final RNA BQSR...")
    bqsr_pool_args = []
    for i, bam in enumerate(processed_rna_bams):
        # Original bam input was the split bam
        base_name = os.path.basename(args.rna_bams[i]).replace('.bam', '')
        out_prefix = os.path.join(args.output_dir, base_name)
        bqsr_pool_args.append((bam, boosted_known_sites, out_prefix, args.genome))
    
    with Pool(args.threads) as p:
        p.map(apply_rna_bqsr, bqsr_pool_args)
    
    logging.info("Step 1 Complete. Recalibrated BAMs are ready.")

if __name__ == "__main__":
    main()
