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

def revised_convert_vcf(input_vcf, output_file):
    """
    Python implementation of revised_convertVCF.sh logic.
    """
    logging.info(f"Converting VCF: {input_vcf} -> {output_file}")
    try:
        with open(input_vcf, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 10:
                    continue
                
                # VCF columns: CHROM(0) POS(1) ID(2) REF(3) ALT(4) QUAL(5) FILTER(6) INFO(7) FORMAT(8) SAMPLE(9)...
                chrom = parts[0]
                pos = parts[1]
                ref = parts[3]
                alt = parts[4]
                sample_data = parts[9] # Assuming first sample is the one we want
                
                # Format usually GT:AD:AF:DP...
                # We need AD (Allele Depth)
                format_fields = parts[8].split(':')
                sample_fields = sample_data.split(':')
                
                try:
                    ad_index = format_fields.index('AD')
                    ad_val = sample_fields[ad_index]
                    ref_depth, alt_depth = map(int, ad_val.split(','))
                    total_depth = ref_depth + alt_depth
                    
                    depths_field = f"{total_depth},{alt_depth}"
                    alt_fraction = 0.01 # Hardcoded in original script
                    
                    # Output: CHROM POS DepthField REF ALT AltFraction
                    outfile.write(f"{chrom}\t{pos}\t{depths_field}\t{ref}\t{alt}\t{alt_fraction}\n")
                except (ValueError, IndexError) as e:
                    logging.warning(f"Skipping line due to parse error: {line.strip()} - {e}")
                    continue
    except Exception as e:
        logging.error(f"Failed to convert VCF: {e}")
        raise e

def calculate_contamination_task(args):
    bam, gnomad, output_dir, genome = args
    base_name = os.path.basename(bam).replace('.bam', '')
    pileup_table = os.path.join(output_dir, f"{base_name}_pileups.table")
    contam_table = os.path.join(output_dir, f"{base_name}_contamination.table")
    
    # Limit memory to 128G as requested
    cmd = f"gatk --java-options '-Xmx128G' GetPileupSummaries -R {genome} -I {bam} -V {gnomad} -L {gnomad} -O {pileup_table}"
    run_command(cmd)
    
    cmd = f"gatk --java-options '-Xmx128G' CalculateContamination -I {pileup_table} -O {contam_table}"
    run_command(cmd)
    
    return contam_table

def main():
    try:
        _main()
    except Exception as e:
        logging.error(f"Pipeline failed with error: {e}")
        sys.exit(1)

def _main():
    parser = argparse.ArgumentParser(description="Step 2: Variant Discovery (Mutect2, Filtering, Custom Filters)")
    parser.add_argument('--rna_bams', nargs='+', required=True, help='List of Recalibrated RNA BAM files (from Step 1)')
    parser.add_argument('--dna_bam', required=True, help='Recalibrated DNA BAM file (from Step 1)')
    parser.add_argument('--genome', required=True, help='Reference genome FASTA')
    parser.add_argument('--gnomad', required=True, help='gnomAD VCF for germline resource')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    parser.add_argument('--prefix', default='sample', help='Output prefix')
    parser.add_argument('--scripts_dir', help='Directory containing filter_homopolymer_nucleotides.py, etc. Defaults to script location.')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads/processes')
    
    args = parser.parse_args()

    if args.scripts_dir is None:
        args.scripts_dir = os.path.dirname(os.path.abspath(__file__))

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
        
    prefix_path = os.path.join(args.output_dir, args.prefix)

    # 1. Calculate Contamination
    logging.info(f"Calculating Contamination (Parallel, max {args.threads} processes)...")
    
    contam_pool_args = []
    for bam in args.rna_bams:
        contam_pool_args.append((bam, args.gnomad, args.output_dir, args.genome))
        
    with Pool(processes=args.threads) as p:
        contamination_tables = p.map(calculate_contamination_task, contam_pool_args)

    # 2. Mutect2 (Final Run)
    logging.info("Running Mutect2...")
    rna_inputs = " ".join([f"-I {bam}" for bam in args.rna_bams])
    contam_args = " ".join([f"--contamination-table {ct}" for ct in contamination_tables])
    
    # Automatically detect DNA sample name from BAM header
    logging.info(f"Detecting sample name from DNA BAM: {args.dna_bam}")
    try:
        # Use samtools to get the header and grep for RG line
        # Assuming there is one RG line or we take the first one's SM tag
        # output format: @RG ID:xxx SM:sample_name ...
        cmd_detect = f"samtools view -H {args.dna_bam} | grep '^@RG' | head -n 1"
        rg_line = subprocess.check_output(cmd_detect, shell=True, executable='/bin/bash').decode('utf-8').strip()
        
        # Parse SM tag
        import re
        match = re.search(r'SM:([^\s]+)', rg_line)
        if match:
            normal_sample_name = match.group(1)
            logging.info(f"Detected DNA Normal Sample Name: {normal_sample_name}")
        else:
            logging.warning("Could not find SM tag in RG line, defaulting to 'DNA_processed'")
            normal_sample_name = "DNA_processed"
    except Exception as e:
        logging.warning(f"Failed to detect sample name: {e}. Defaulting to 'DNA_processed'")
        normal_sample_name = "DNA_processed"
    
    cmd = f"gatk Mutect2 -R {args.genome} {rna_inputs} -I {args.dna_bam} -normal {normal_sample_name} --germline-resource {args.gnomad} -O {prefix_path}.temp.vcf"
    run_command(cmd)

    # 3. Filter Mutect Calls
    logging.info("Filtering Mutect Calls...")
    cmd = f"gatk FilterMutectCalls -V {prefix_path}.temp.vcf -R {args.genome} {contam_args} --min-median-base-quality 12 --max-events-in-region 4 -O {prefix_path}.temp2.vcf"
    run_command(cmd)

    # 4. Select Variants (PASS and SNP)
    cmd = f"bcftools view -f PASS {prefix_path}.temp2.vcf > {prefix_path}.vcf"
    run_command(cmd)
    
    cmd = f"gatk SelectVariants -R {args.genome} -V {prefix_path}.vcf --select-type-to-include SNP -O {prefix_path}.SNV.vcf"
    run_command(cmd)

    # 5. Custom Conversions and Filters
    logging.info("Running Custom Filters...")
    
    # Convert VCF
    revised_convert_vcf(f"{prefix_path}.SNV.vcf", f"{prefix_path}_SNPIR.SNV.vcf")
    
    # Filter Homopolymers (Calling external script)
    filter_script = os.path.join(args.scripts_dir, "filter_homopolymer_nucleotides.pl")
    cmd = f"perl {filter_script} -infile {prefix_path}_SNPIR.SNV.vcf -outfile {prefix_path}.homo.vcf -refgenome {args.genome}"
    run_command(cmd)
    
    # PBLAT Steps
    # Merge RNA BAMs for PBLAT reference
    merged_bam = os.path.join(args.output_dir, f"{args.prefix}all.bam")
    sorted_bam = os.path.join(args.output_dir, f"{args.prefix}all_sorted.bam")
    
    cmd = f"samtools merge -f {merged_bam} {' '.join(args.rna_bams)}"
    run_command(cmd)
    
    cmd = f"samtools sort -@ {args.threads} -o {sorted_bam} {merged_bam}"
    run_command(cmd)
    
    cmd = f"samtools index {sorted_bam}"
    run_command(cmd)
    
    # Run PBLAT script
    pblat_script = os.path.join(args.scripts_dir, "pblat_candidates_ln.pl")
    pblat_vcf = f"{prefix_path}.pblat.vcf"
    cmd = f"perl {pblat_script} -infile {prefix_path}.homo.vcf -outfile {pblat_vcf} -bamfile {sorted_bam} -refgenome {args.genome} -minbasequal 5 -threads {args.threads}"
    run_command(cmd)
    
    # Cleanup merged bams
    run_command(f"rm {merged_bam} {sorted_bam} {sorted_bam}.bai")
    
    # Create knownedit.vcf (awk step)
    known_edit_vcf = f"{prefix_path}.knownedit.vcf"
    # awk '{OFS="\t";$2=$2-1"\t"$2;print $0}'
    # This awk command converts 1-based VCF pos to 0-based interval (BED-like) but keeps other columns?
    # Wait, the original script says: awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' $prefix.pblat.vcf > $prefix.knownedit.vcf
    # This looks like it inserts a column or modifies column 2 to be "start end". 
    # VCF is 1-based. bedtools intersect -a VCF -b custom_bed.
    # The output of PBLAT script might not be a standard VCF?
    # Assuming the original awk logic is correct for the output of pblat script.
    
    cmd = f"awk '{{OFS=\"\\t\";$2=$2-1\"\\t\"$2;print $0}}' {pblat_vcf} > {known_edit_vcf}"
    run_command(cmd)

    # Final Intersection
    # bedtools intersect -a $prefix.SNV.vcf -b $prefix.knownedit.vcf -wa -header > $prefix.final.vcf
    final_vcf = f"{prefix_path}.final.vcf"
    cmd = f"bedtools intersect -a {prefix_path}.SNV.vcf -b {known_edit_vcf} -wa -header > {final_vcf}"
    run_command(cmd)
    
    logging.info(f"Step 2 Complete. Final VCF: {final_vcf}")

if __name__ == "__main__":
    main()
