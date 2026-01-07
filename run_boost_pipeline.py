#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
run_boost_pipeline.py

This script is a Python implementation of the `cardes_boost.sh` pipeline.
It performs the initial steps of the CADRES workflow to generate a reliable
set of known variant sites. This involves:
1. Calibrating RNA BAM files.
2. Calibrating the DNA (WGS) BAM file.
3. Running GATK Mutect2 with both RNA and DNA samples to call variants.
4. Filtering these variants to create a high-confidence VCF file for later use.

This script is designed to be modular, with each major step encapsulated
in its own function.
"""

import os
import sys
import subprocess
import logging
import argparse
from typing import List, Dict, Optional, Union

# --- Configuration for Logging ---
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] - %(message)s",
    stream=sys.stdout,
    datefmt="%Y-%m-%d %H:%M:%S"
)

# --- Utility Functions ---

def run_command(command: List[str], cwd: str, outfile: Optional[str] = None) -> bool:
    """
    Executes a command, logs it, captures output, and checks for errors.

    Args:
        command (List[str]): The command and its arguments as a list.
        cwd (str): The working directory to run the command in.
        outfile (Optional[str]): Path to a file to redirect stdout to.

    Returns:
        bool: True if the command succeeded, False otherwise.
    """
    logging.info(f"Executing command: {' '.join(command)}")
    logging.info(f"Working directory: {os.path.abspath(cwd)}")
    
    try:
        if outfile:
            with open(outfile, "w") as f_out:
                result = subprocess.run(
                    command,
                    cwd=cwd,
                    check=True,
                    stdout=f_out,
                    stderr=subprocess.PIPE,
                    text=True
                )
        else:
            result = subprocess.run(
                command,
                cwd=cwd,
                check=True,
                capture_output=True,
                text=True
            )
        
        # Log stdout/stderr if they exist
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
        logging.error(f"Command: {' '.join(e.cmd)}")
        if e.stdout:
            logging.error(f"STDOUT:\n{e.stdout.strip()}")
        if e.stderr:
            logging.error(f"STDERR:\n{e.stderr.strip()}")
        return False

def ensure_dir_exists(path: str):
    """Creates a directory if it doesn't exist."""
    if not os.path.exists(path):
        logging.info(f"Creating directory: {path}")
        os.makedirs(path)

# --- Pipeline Step Functions ---

def calibrate_rna_bams(config: Dict[str, Union[str, List[str]]]) -> Optional[List[str]]:
    """
    Calibrates RNA BAM files using the bam_calibration_RNA_boost.py script.
    This step prepares the RNA reads by reordering, marking duplicates, and
    splitting N-cigar reads.

    Args:
        config (Dict): A dictionary containing configuration parameters.
            Required keys: 'rna_bam_files', 'boost_working_dir', 'script_dir', 'genome'.

    Returns:
        Optional[List[str]]: A list of paths to the generated split BAM files,
                             or None if the step fails.
    """
    logging.info("--- Step 1: Calibrating RNA BAMs ---")
    ensure_dir_exists(str(config['boost_working_dir']))
    split_bam_files = []

    for bam_file in config['rna_bam_files']:
        sample_name = os.path.basename(bam_file).replace('.bam', '')
        output_prefix = os.path.join(str(config['boost_working_dir']), str(sample_name))
        split_bam_path = f"{output_prefix}_split.bam"
        
        cmd = [
            "python",
            os.path.join(str(config['script_dir']), "bam_calibration_RNA_boost.py"),
            "--bam", bam_file,
            "--output", output_prefix,
            "--genome", str(config['genome'])
        ]
        
        if not run_command(cmd, cwd=str(config['boost_working_dir'])):
            logging.error(f"RNA BAM calibration failed for sample {sample_name}.")
            return None
            
        split_bam_files.append(split_bam_path)
    
    logging.info("All RNA BAMs calibrated successfully.")
    return split_bam_files


def calibrate_dna_bam(config: Dict[str, str]) -> Optional[str]:
    """
    Calibrates the DNA BAM file, including Base Quality Score Recalibration (BQSR).

    Args:
        config (Dict): A dictionary containing configuration parameters.
            Required keys: 'dna_bsqr_dir', 'dna_bam', 'dna_name', 'script_dir',
                           'genome', 'known_snv'.

    Returns:
        Optional[str]: The path to the recalibrated DNA BAM file, or None on failure.
    """
    logging.info("--- Step 2: Calibrating DNA BAM ---")
    ensure_dir_exists(config['dna_bsqr_dir'])
    
    cmd = [
        "python",
        os.path.join(config['script_dir'], "bam_calibration_DNA.py"),
        "--bam", config['dna_bam'],
        "--output", os.path.join(config['dna_bsqr_dir'], config['dna_name']),
        "--genome", config['genome'],
        "--known", config['known_snv']
    ]

    if not run_command(cmd, cwd=config['dna_bsqr_dir']):
        logging.error("DNA BAM calibration failed.")
        return None
        
    recalibrated_bam = os.path.join(config['dna_bsqr_dir'], f"{config['dna_name']}_recalibration.bam")
    logging.info(f"DNA BAM calibrated successfully: {recalibrated_bam}")
    return recalibrated_bam


def run_mutect2_boost(rna_bam_list: List[str], recal_dna_bam: str, config: Dict[str, str]) -> Optional[str]:
    """
    Runs GATK Mutect2 using the calibrated RNA and DNA BAMs to generate
    a VCF file for creating a panel of normals.

    Args:
        rna_bam_list (List[str]): List of paths to calibrated RNA BAM files.
        recal_dna_bam (str): Path to the calibrated DNA BAM file.
        config (Dict): Configuration dictionary.
            Required keys: 'boost_working_dir', 'prefix', 'genome', 'dna_name'.

    Returns:
        Optional[str]: The path to the raw Mutect2 VCF output, or None on failure.
    """
    logging.info("--- Step 3: Running GATK Mutect2 for Boosting ---")
    working_dir = config['boost_working_dir']
    output_vcf = os.path.join(working_dir, f"{config['prefix']}.mutect.vcf")

    cmd = ["gatk", "Mutect2", "-R", config['genome']]
    
    for rna_bam in rna_bam_list:
        cmd.extend(["-I", rna_bam])
        
    cmd.extend([
        "-I", recal_dna_bam,
        "-normal", config['dna_name'],
        "-O", output_vcf
    ])

    if not run_command(cmd, cwd=working_dir):
        logging.error("GATK Mutect2 failed.")
        return None
        
    logging.info(f"GATK Mutect2 completed. Output: {output_vcf}")
    return output_vcf


def filter_and_index_vcf(mutect_vcf: str, config: Dict[str, str]) -> Optional[str]:
    """
    Filters the raw Mutect2 VCF, selects PASSing variants, and indexes the result.

    Args:
        mutect_vcf (str): Path to the raw Mutect2 VCF file.
        config (Dict): Configuration dictionary.
            Required keys: 'boost_working_dir', 'prefix', 'genome'.

    Returns:
        Optional[str]: The path to the final, filtered, and indexed VCF, or None on failure.
    """
    logging.info("--- Step 4: Filtering and Indexing VCF ---")
    working_dir = config['boost_working_dir']
    prefix = config['prefix']
    
    # 1. FilterMutectCalls
    filtered_vcf = os.path.join(working_dir, f"{prefix}_boost.filter.vcf")
    cmd_filter = [
        "gatk", "FilterMutectCalls",
        "-V", mutect_vcf,
        "-R", config['genome'],
        "-O", filtered_vcf,
        "--max-events-in-region", "4"
    ]
    if not run_command(cmd_filter, cwd=working_dir):
        logging.error("GATK FilterMutectCalls failed.")
        return None

    # 2. Select PASSing variants with bcftools
    known_vcf = os.path.join(working_dir, f"{prefix}_boost.known.vcf")
    cmd_view = ["bcftools", "view", "-f", "PASS", filtered_vcf]
    if not run_command(cmd_view, cwd=working_dir, outfile=known_vcf):
        logging.error("bcftools view failed.")
        return None
        
    # 3. Index the final VCF
    cmd_index = ["gatk", "IndexFeatureFile", "-I", known_vcf]
    if not run_command(cmd_index, cwd=working_dir):
        logging.error("GATK IndexFeatureFile failed.")
        return None
        
    logging.info(f"VCF filtering and indexing complete. Final VCF: {known_vcf}")
    return known_vcf

def run_boost_pipeline(config: Dict[str, Union[str, List[str]]]) -> Optional[str]:
    """
    Main function to run the boost pipeline programmatically.
    
    Args:
        config: Dictionary containing all configuration parameters
        
    Returns:
        Path to the final VCF file if successful, None otherwise
    """
    logging.info("=== Starting CADRES Boost Pipeline ===")
    
    # 1. Calibrate RNA BAMs
    rna_split_bams = calibrate_rna_bams(config)
    if not rna_split_bams:
        logging.error("Pipeline aborted due to RNA calibration failure.")
        return None
        
    # 2. Calibrate DNA BAM
    recal_dna_bam = calibrate_dna_bam(config)
    if not recal_dna_bam:
        logging.error("Pipeline aborted due to DNA calibration failure.")
        return None
        
    # 3. Run Mutect2
    mutect_vcf = run_mutect2_boost(rna_split_bams, recal_dna_bam, config)
    if not mutect_vcf:
        logging.error("Pipeline aborted due to GATK Mutect2 failure.")
        return None

    # 4. Filter and Index VCF
    final_vcf = filter_and_index_vcf(mutect_vcf, config)
    if not final_vcf:
        logging.error("Pipeline aborted due to VCF filtering/indexing failure.")
        return None
        
    logging.info("=== CADRES Boost Pipeline Completed Successfully! ===")
    logging.info(f"Final known sites VCF is available at: {final_vcf}")
    
    return final_vcf


def main():
    """
    Main function to parse arguments and run the boost pipeline.
    """
    parser = argparse.ArgumentParser(description="Python runner for the CADRES boost pipeline.")
    
    # --- Path Arguments ---
    parser.add_argument("--script_dir", required=True, help="Directory containing the helper python scripts (e.g., bam_calibration_*.py).")
    parser.add_argument("--genome", required=True, help="Path to the reference genome FASTA file.")
    parser.add_argument("--known_snv", required=True, help="Path to a VCF file of known SNVs (e.g., dbSNP).")
    parser.add_argument("--dna_bam", required=True, help="Full path to the input DNA (WGS) BAM file.")
    parser.add_argument("--boost_working_dir", required=True, help="Main working directory for this boost pipeline.")
    parser.add_argument("--dna_bsqr_dir", required=True, help="Working directory for DNA BQSR steps.")
    
    # --- Naming and Sample Arguments ---
    parser.add_argument("--prefix", required=True, help="A general prefix for output files (e.g., 'ProjectX_V1').")
    parser.add_argument("--dna_name", required=True, help="The sample name for the DNA sample as used in the BAM's RG tag.")
    parser.add_argument("--rna_bam_files", required=True, nargs='+', help="A space-separated list of full paths to RNA BAM files.")

    args = parser.parse_args()

    # --- Configuration Dictionary ---
    config = {
        "script_dir": args.script_dir,
        "genome": args.genome,
        "known_snv": args.known_snv,
        "dna_bam": args.dna_bam,
        "boost_working_dir": args.boost_working_dir,
        "dna_bsqr_dir": args.dna_bsqr_dir,
        "prefix": args.prefix,
        "dna_name": args.dna_name,
        "rna_bam_files": args.rna_bam_files
    }

    # --- Execute Pipeline ---
    logging.info("=== Starting CADRES Boost Pipeline ===")
    
    # 1. Calibrate RNA BAMs
    rna_split_bams = calibrate_rna_bams(config)
    if not rna_split_bams:
        logging.error("Pipeline aborted due to RNA calibration failure.")
        sys.exit(1)
        
    # 2. Calibrate DNA BAM
    recal_dna_bam = calibrate_dna_bam(config)
    if not recal_dna_bam:
        logging.error("Pipeline aborted due to DNA calibration failure.")
        sys.exit(1)
        
    # 3. Run Mutect2
    mutect_vcf = run_mutect2_boost(rna_split_bams, recal_dna_bam, config)
    if not mutect_vcf:
        logging.error("Pipeline aborted due to GATK Mutect2 failure.")
        sys.exit(1)

    # 4. Filter and Index VCF
    final_vcf = filter_and_index_vcf(mutect_vcf, config)
    if not final_vcf:
        logging.error("Pipeline aborted due to VCF filtering/indexing failure.")
        sys.exit(1)
        
    logging.info("=== CADRES Boost Pipeline Completed Successfully! ===")
    logging.info(f"Final known sites VCF is available at: {final_vcf}")


if __name__ == "__main__":
    main()

