#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
run_dvr_pipeline.py

This script is a Python implementation of the `cardes_DVR.sh` pipeline.
It performs the main steps of the CADRES workflow to identify Differential
Variant Rate (DVR) sites between two sample groups. The pipeline includes:
1. Base Quality Score Recalibration (BQSR) for RNA BAMs.
2. Contamination estimation.
3. RNA variant calling with Mutect2.
4. A series of custom filters (homopolymer, PBLAT) on the called variants.
5. Differential analysis using a custom MATS-like workflow.
6. Final annotation of significant DVR sites.

This script is designed to be modular and can be called from a master workflow.
"""

import os
import sys
import subprocess
import logging
import argparse
from typing import List, Dict, Optional, Union, Tuple

# --- Configuration for Logging ---
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] - %(message)s",
    stream=sys.stdout,
    datefmt="%Y-%m-%d %H:%M:%S"
)

# --- Utility Functions (can be shared with run_boost_pipeline.py) ---

def run_command(command: List[str], cwd: str, outfile: Optional[str] = None) -> bool:
    """
    Executes a command, logs it, captures output, and checks for errors.
    """
    logging.info(f"Executing command: {' '.join(command)}")
    logging.info(f"Working directory: {os.path.abspath(cwd)}")
    
    try:
        # Using shell=True for commands that might need it (like those with pipes in shell scripts)
        # but it's safer to pass args as a list. For simplicity and consistency with the
        # original shell script's behavior, we'll construct the string command.
        cmd_str = ' '.join(command)
        if outfile:
            with open(outfile, "w") as f_out:
                result = subprocess.run(
                    cmd_str,
                    cwd=cwd,
                    check=True,
                    stdout=f_out,
                    stderr=subprocess.PIPE,
                    text=True,
                    shell=True
                )
        else:
            result = subprocess.run(
                cmd_str,
                cwd=cwd,
                check=True,
                capture_output=True,
                text=True,
                shell=True
            )
        
        if not outfile and result.stdout:
            logging.info(f"STDOUT:\n{result.stdout.strip()}")
        if result.stderr:
            logging.info(f"STDERR:\n{result.stderr.strip()}")
            
        return True
    except FileNotFoundError:
        logging.error(f"Command not found: {command[0]}. Please ensure it is in your PATH.")
        return False
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed with exit code {e.returncode}")
        logging.error(f"Command: {e.cmd}")
        if hasattr(e, 'stdout') and e.stdout:
            logging.error(f"STDOUT:\n{e.stdout.strip()}")
        if hasattr(e, 'stderr') and e.stderr:
            logging.error(f"STDERR:\n{e.stderr.strip()}")
        return False

def ensure_dir_exists(path: str):
    """Creates a directory if it doesn't exist."""
    if not os.path.exists(path):
        logging.info(f"Creating directory: {path}")
        os.makedirs(path)

# --- Pipeline Step Functions ---

def bqsr_rna(config: Dict) -> Optional[Tuple[List[str], List[str]]]:
    """
    Step 1: Perform Base Quality Score Recalibration (BQSR) on RNA BAMs.
    """
    logging.info("--- Step 1: Performing BQSR on RNA BAMs ---")
    rna_bqsr_dir = config['rna_bqsr_dir']
    ensure_dir_exists(rna_bqsr_dir)

    recal_reports = []
    recalibrated_bams = []
    
    # BaseRecalibrator for each sample
    for sample in config['sample_list']:
        input_bam = os.path.join(config['boost_working_dir'], f"{sample}_split.bam")
        report_grp = os.path.join(rna_bqsr_dir, f"{sample}_recalibration_report.grp")
        cmd_br = [
            "gatk", "BaseRecalibrator",
            "-I", input_bam,
            "-R", config['genome'],
            "-O", report_grp,
            "--known-sites", config['known_vcf_from_boost']
        ]
        if not run_command(cmd_br, cwd=rna_bqsr_dir):
            return None
        recal_reports.append(report_grp)

    # GatherBQSRReports
    gathered_report = os.path.join(rna_bqsr_dir, f"{config['prefix']}_recalibration_report.grp")
    cmd_gather = ["gatk", "GatherBQSRReports"]
    for report in recal_reports:
        cmd_gather.extend(["-I", report])
    cmd_gather.extend(["-O", gathered_report])
    if not run_command(cmd_gather, cwd=rna_bqsr_dir):
        return None

    # ApplyBQSR for each sample
    for sample in config['sample_list']:
        input_bam = os.path.join(config['boost_working_dir'], f"{sample}_split.bam")
        output_bam = os.path.join(rna_bqsr_dir, f"{sample}_recalibration.bam")
        cmd_apply = [
            "gatk", "ApplyBQSR",
            "-R", config['genome'],
            "-I", input_bam,
            "--bqsr-recal-file", gathered_report,
            "-O", output_bam
        ]
        if not run_command(cmd_apply, cwd=rna_bqsr_dir):
            return None
        recalibrated_bams.append(output_bam)

    logging.info("RNA BQSR completed successfully.")
    return recalibrated_bams, recal_reports


def calculate_contamination(recal_bams: List[str], config: Dict) -> Optional[List[str]]:
    """
    Step 2: Calculate contamination for each RNA sample.
    """
    logging.info("--- Step 2: Calculating Contamination ---")
    rna_bqsr_dir = config['rna_bqsr_dir']
    contamination_tables = []
    
    for i, recal_bam in enumerate(recal_bams):
        sample_name = config['sample_list'][i]
        pileups_table = os.path.join(rna_bqsr_dir, f"{sample_name}_pileups.table")
        contam_table = os.path.join(rna_bqsr_dir, f"{sample_name}_contamination.table")
        
        cmd_pileup = [
            "gatk", "--java-options", "'-Xmx200G'", "GetPileupSummaries",
            "-I", recal_bam,
            "-V", config['genome_ad'],
            "-L", config['genome_ad'],
            "-O", pileups_table
        ]
        if not run_command(cmd_pileup, cwd=rna_bqsr_dir):
            return None
            
        cmd_contam = [
            "gatk", "--java-options", "'-Xmx200G'", "CalculateContamination",
            "-I", pileups_table,
            "-O", contam_table
        ]
        if not run_command(cmd_contam, cwd=rna_bqsr_dir):
            return None
        contamination_tables.append(contam_table)

    logging.info("Contamination calculation completed successfully.")
    return contamination_tables


def call_and_filter_rv(recal_bams: List[str], contam_tables: List[str], config: Dict) -> Optional[str]:
    """
    Step 3: Call RNA Variants (RV) with Mutect2 and filter them.
    """
    logging.info("--- Step 3: Calling and Filtering RNA Variants ---")
    rv_dir = config['rv_working_dir']
    ensure_dir_exists(rv_dir)
    prefix = config['prefix']

    # Mutect2
    temp_vcf = os.path.join(rv_dir, f"{prefix}.temp.vcf")
    cmd_mutect2 = ["gatk", "Mutect2", "-R", config['genome']]
    for bam in recal_bams:
        cmd_mutect2.extend(["-I", bam])
    cmd_mutect2.extend([
        "-I", config['dna_bam'],
        "-normal", config['normal_sample_name'],
        "--germline-resource", config['genome_ad'],
        "-O", temp_vcf
    ])
    if not run_command(cmd_mutect2, cwd=rv_dir):
        return None

    # FilterMutectCalls
    temp2_vcf = os.path.join(rv_dir, f"{prefix}.temp2.vcf")
    cmd_filter = [
        "gatk", "FilterMutectCalls",
        "-V", temp_vcf,
        "-R", config['genome']
    ]
    for table in contam_tables:
        cmd_filter.extend(["--contamination-table", table])
    cmd_filter.extend([
        "--min-median-base-quality", "12",
        "--max-events-in-region", "4",
        "-O", temp2_vcf
    ])
    if not run_command(cmd_filter, cwd=rv_dir):
        return None

    # bcftools view
    pass_vcf = os.path.join(rv_dir, f"{prefix}.vcf")
    cmd_view = ["bcftools", "view", "-f", "PASS", temp2_vcf]
    if not run_command(cmd_view, cwd=rv_dir, outfile=pass_vcf):
        return None

    # SelectVariants
    snv_vcf = os.path.join(rv_dir, f"{prefix}.SNV.vcf")
    cmd_select = [
        "gatk", "SelectVariants",
        "-R", config['genome'],
        "-V", pass_vcf,
        "--select-type-to-include", "SNP",
        "-O", snv_vcf
    ]
    if not run_command(cmd_select, cwd=rv_dir):
        return None
        
    logging.info(f"Variant calling and filtering complete. Final SNV VCF: {snv_vcf}")
    return snv_vcf


def custom_filter_variants(snv_vcf: str, recal_bams: List[str], config: Dict) -> Optional[str]:
    """
    Step 4: Apply a series of custom filters (homopolymer, pblat).
    """
    logging.info("--- Step 4: Applying Custom Filters ---")
    rv_dir = config['rv_working_dir']
    prefix = config['prefix']
    script_dir = config['script_dir']
    
    # revised_convertVCF
    snv_snipr_vcf = os.path.join(rv_dir, f"{prefix}_SNPIR.SNV.vcf")
    cmd_convert = [
        "python", os.path.join(script_dir, "revised_convertVCF.py"),
        snv_vcf, snv_snipr_vcf
    ]
    if not run_command(cmd_convert, cwd=rv_dir):
        return None

    # filter_homopolymer_nucleotides
    homo_vcf = os.path.join(rv_dir, f"{prefix}.homo.vcf")
    cmd_homo = [
        "python", os.path.join(script_dir, "filter_homopolymer_nucleotides.py"),
        "--infile", snv_snipr_vcf,
        "--outfile", homo_vcf,
        "--refgenome", config['genome']
    ]
    if not run_command(cmd_homo, cwd=rv_dir):
        return None

    # Prepare for pblat: merge, sort, index BAMs
    merged_bam = os.path.join(rv_dir, f"{prefix}all.bam")
    sorted_bam = os.path.join(rv_dir, f"{prefix}all_sorted.bam")
    cmd_merge = ["samtools", "merge", "-f", merged_bam] + recal_bams
    cmd_sort = ["samtools", "sort", "-@", str(config['threads']), "-o", sorted_bam, merged_bam]
    cmd_index = ["samtools", "index", sorted_bam]
    
    if not run_command(cmd_merge, cwd=rv_dir) or \
       not run_command(cmd_sort, cwd=rv_dir) or \
       not run_command(cmd_index, cwd=rv_dir):
        return None

    # pblat_candidates
    pblat_vcf = os.path.join(rv_dir, f"{prefix}.pblat.vcf")
    cmd_pblat = [
        "python", os.path.join(script_dir, "pblat_candidates_filter.py"),
        "--infile", homo_vcf,
        "--outfile", pblat_vcf,
        "--bamfile", sorted_bam,
        "--refgenome", config['genome'],
        "--minbasequal", "5",
        "--threads", str(config['threads'])
    ]
    if not run_command(cmd_pblat, cwd=rv_dir):
        return None
    
    # Cleanup temp BAMs
    os.remove(merged_bam)
    os.remove(sorted_bam)
    os.remove(sorted_bam + ".bai")

    # awk and bedtools intersect
    knownedit_vcf = os.path.join(rv_dir, f"{prefix}.knownedit.vcf")
    cmd_awk = ["awk", "'{OFS=\"\\t\";$2=$2-1\"\\t\"$2;print $0}'", pblat_vcf]
    if not run_command(cmd_awk, cwd=rv_dir, outfile=knownedit_vcf):
        return None

    final_vcf = os.path.join(rv_dir, f"{prefix}.final.vcf")
    cmd_intersect = [
        "bedtools", "intersect",
        "-a", snv_vcf,
        "-b", knownedit_vcf,
        "-wa", "-header"
    ]
    if not run_command(cmd_intersect, cwd=rv_dir, outfile=final_vcf):
        return None

    logging.info(f"Custom filtering complete. Final filtered VCF: {final_vcf}")
    return final_vcf


def run_dvr_analysis(final_vcf: str, recal_bams: List[str], config: Dict) -> bool:
    """
    Step 5: Run the core DVR analysis using mpileup and custom scripts.
    """
    logging.info("--- Step 5: Running DVR Analysis ---")
    rv_dir = config['rv_working_dir']
    prefix = config['prefix']
    script_dir = config['script_dir']
    
    # Prepare separated BAM files
    all_sep_bams = []
    for bam in recal_bams:
        base_name = os.path.splitext(bam)[0]
        first_end_bam = f"{base_name}_seperate_firstEnd.bam"
        second_end_bam = f"{base_name}_seperate_secondEnd.bam"
        run_command(["samtools", "view", "-@", str(config['threads']), "-f", "64", "-b", bam], cwd=rv_dir, outfile=first_end_bam)
        run_command(["samtools", "view", "-@", str(config['threads']), "-f", "128", "-b", bam], cwd=rv_dir, outfile=second_end_bam)
        all_sep_bams.extend([first_end_bam, second_end_bam])

    # samtools mpileup
    pileup_file = os.path.join(rv_dir, f"{prefix}.pileup")
    cmd_mpileup = [
        "samtools", "mpileup",
        "-B", "-d", "100000",
        "-f", config['genome'],
        "-l", final_vcf,
        "-q", "30", "-Q", "17", "-a",
        "-o", pileup_file
    ] + all_sep_bams
    if not run_command(cmd_mpileup, cwd=rv_dir):
        return False
    
    # Cleanup separated bams
    for bam in all_sep_bams:
        os.remove(bam)

    # vcf_to_mats_input
    inc_file = os.path.join(rv_dir, f"{prefix}.inc.txt")
    cmd_vcf2mats = [
        "python", os.path.join(script_dir, "vcf_to_mats_input_For_Mutect2.py"),
        final_vcf,
        inc_file,
        config['rna_bam_sample1_str'],
        config['rna_bam_sample2_str'],
        "20", "5", "T", pileup_file, "T"
    ]
    if not run_command(cmd_vcf2mats, cwd=rv_dir):
        return False

    # MATS_LRT and FDR
    mats_results_prefix = os.path.join(rv_dir, f"{prefix}_rMATS-DVR_results")
    mats_p_file = f"{mats_results_prefix}_rMATS_Result_P.txt"
    mats_fdr_file = f"{mats_results_prefix}_rMATS_Result_FDR.txt"
    cmd_mats = [
        "python", os.path.join(script_dir, "MATS_LRT.py"),
        inc_file, mats_results_prefix, "4", "0.0001"
    ]
    cmd_fdr = [
        "python", os.path.join(script_dir, "FDR.py"),
        mats_p_file, mats_fdr_file
    ]
    if not run_command(cmd_mats, cwd=rv_dir) or not run_command(cmd_fdr, cwd=rv_dir):
        return False

    # snv_annotation
    final_annotated_results = os.path.join(rv_dir, f"{prefix}_rMATS-DVR_Result.txt")
    summary_file = os.path.join(rv_dir, f"{prefix}rMATS-DVR_Result_summary.txt")
    cmd_anno = [
        "python", os.path.join(script_dir, "snv_annotation.py"),
        "--input", mats_fdr_file,
        "--output", final_annotated_results,
        "--summary", summary_file,
        "--label1", config['label1'],
        "--label2", config['label2'],
        "--snp", config['known_snv_anno'],
        "--editing", config['known_editing'],
        "--gene", config['gene_anno']
    ]
    if not run_command(cmd_anno, cwd=rv_dir):
        return False

    logging.info(f"DVR analysis complete. Final results: {final_annotated_results}")
    return True

def main():
    """
    Main function to parse arguments and run the DVR pipeline.
    """
    parser = argparse.ArgumentParser(description="Python runner for the CADRES DVR pipeline.")
    
    # --- Directory and Path Arguments ---
    parser.add_argument("--script_dir", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--genome_ad", required=True, help="Path to gnomAD VCF for contamination check.")
    parser.add_argument("--boost_working_dir", required=True)
    parser.add_argument("--rna_bqsr_dir", required=True)
    parser.add_argument("--dna_bqsr_dir", required=True)
    parser.add_argument("--rv_working_dir", required=True)
    parser.add_argument("--dna_bam", required=True)
    parser.add_argument("--known_vcf_from_boost", required=True, help="Path to the known.vcf generated by the boost pipeline.")

    # --- Annotation File Arguments ---
    parser.add_argument("--known_snv_anno", required=True, help="dbSNP VCF for annotation.")
    parser.add_argument("--known_editing", required=True, help="Known RNA editing sites file (e.g., REDIportal).")
    parser.add_argument("--gene_anno", required=True, help="Gene annotation file (refGene.txt format).")

    # --- Naming, Sample, and Config Arguments ---
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--normal_sample_name", required=True, help="Name of the normal DNA sample.")
    parser.add_argument("--label1", required=True, help="Label for sample group 1.")
    parser.add_argument("--label2", required=True, help="Label for sample group 2.")
    parser.add_argument("--sample_list", required=True, nargs='+', help="Full list of RNA sample names.")
    parser.add_argument("--rna_bam_sample1", required=True, nargs='+', help="List of recalibrated BAMs for sample group 1.")
    parser.add_argument("--rna_bam_sample2", required=True, nargs='+', help="List of recalibrated BAMs for sample group 2.")

    args = parser.parse_args()

    config = vars(args)
    config['rna_bam_sample1_str'] = ','.join(args.rna_bam_sample1)
    config['rna_bam_sample2_str'] = ','.join(args.rna_bam_sample2)

    # --- Execute Pipeline ---
    logging.info("=== Starting CADRES DVR Pipeline ===")
    
    recal_bams, _ = bqsr_rna(config)
    if not recal_bams:
        sys.exit("Pipeline aborted at BQSR step.")

    contam_tables = calculate_contamination(recal_bams, config)
    if not contam_tables:
        sys.exit("Pipeline aborted at contamination calculation step.")
        
    snv_vcf = call_and_filter_rv(recal_bams, contam_tables, config)
    if not snv_vcf:
        sys.exit("Pipeline aborted at variant calling step.")

    final_vcf = custom_filter_variants(snv_vcf, recal_bams, config)
    if not final_vcf:
        sys.exit("Pipeline aborted at custom filtering step.")
        
    if not run_dvr_analysis(final_vcf, recal_bams, config):
        sys.exit("Pipeline aborted at DVR analysis step.")
        
    logging.info("=== CADRES DVR Pipeline Completed Successfully! ===")


if __name__ == "__main__":
    main()
