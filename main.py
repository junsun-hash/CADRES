#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
main.py

Master script for the CADRES pipeline.
This script orchestrates the execution of the two main phases:
1. The 'boost' phase to generate a high-confidence panel of known sites.
2. The 'DVR' phase to identify differentially edited sites between sample groups.

It reads all parameters from a central 'config.json' file.
"""

import json
import logging
import sys
import os
import argparse
from typing import Dict, Any

# Import the main runner functions from the pipeline scripts
# Assumes they are in the same directory or in the PYTHONPATH
import run_boost_pipeline as boost
import run_dvr_pipeline as dvr

def setup_logging(log_file: str):
    """
    Configures logging to write to both a file and the console.
    """
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )

def load_config(config_path: str) -> Dict[str, Any]:
    """
    Loads and parses the JSON configuration file.
    """
    try:
        with open(config_path, 'r') as f:
            config = json.load(f)
        logging.info(f"Configuration loaded successfully from {config_path}")
        return config
    except FileNotFoundError:
        logging.error(f"Configuration file not found at: {config_path}")
        return None
    except json.JSONDecodeError:
        logging.error(f"Error decoding JSON from the configuration file: {config_path}")
        return None

def build_pipeline_configs(master_config: Dict) -> Dict:
    """
    Flattens the master config into a single-level dictionary for easier access
    by the pipeline functions.
    """
    config = {}
    config.update(master_config.get('general', {}))
    config.update(master_config.get('paths', {}))
    config.update(master_config.get('work_dirs', {}))
    
    # BAM file paths for each group
    group1_bam_files = master_config.get('paths', {}).get('group1_rna_bam_files', [])
    group2_bam_files = master_config.get('paths', {}).get('group2_rna_bam_files', [])
    config['group1_rna_bam_files'] = group1_bam_files
    config['group2_rna_bam_files'] = group2_bam_files
    
    # Automatically extract sample names from BAM file paths
    def extract_sample_name(bam_path: str) -> str:
        """Extract sample name from BAM file path by removing .bam extension"""
        return os.path.basename(bam_path).replace('.bam', '')
    
    group1_names = [extract_sample_name(bam) for bam in group1_bam_files]
    group2_names = [extract_sample_name(bam) for bam in group2_bam_files]
    all_rna_samples = group1_names + group2_names
    
    config['sample_list'] = all_rna_samples
    config['group1_names'] = group1_names
    config['group2_names'] = group2_names
    
    # Specific sample config processing (use config values if provided, otherwise use defaults)
    sample_config = master_config.get('samples', {})
    
    # Automatically extract dna_name from input_dna_bam path
    input_dna_bam = config.get('input_dna_bam', '')
    dna_name = os.path.basename(input_dna_bam).replace('.bam', '') if input_dna_bam else 'dna_sample'
    config['dna_name'] = dna_name
    
    config['label1'] = sample_config.get('group1_label', 'Group1')
    config['label2'] = sample_config.get('group2_label', 'Group2')
    config['normal_sample_name'] = config['dna_name'] # for DVR step
    
    # Rename for clarity in boost step
    config['known_snv'] = config.get('known_snv_boost')
    
    return config

def main_workflow(config_path: str):
    """
    The main orchestrator for the entire CADRES workflow.
    """
    master_config = load_config(config_path)
    if not master_config:
        sys.exit(1)

    # Setup logging using the file path from the config
    log_file = master_config.get('general', {}).get('log_file', 'pipeline.log')
    setup_logging(log_file)
    
    # Build a flat config dictionary to pass to functions
    pipeline_config = build_pipeline_configs(master_config)

    # --- Phase 1: Run Boost Pipeline ---
    logging.info(">>> STARTING PHASE 1: BOOST PIPELINE <<<")
    
    # Rename keys to match boost script function signature
    boost_config = pipeline_config.copy()
    boost_config['dna_bsqr_dir'] = boost_config['dna_bqsr_dir']
    boost_config['rna_bam_files'] = boost_config['group1_rna_bam_files'] + boost_config['group2_rna_bam_files']
    boost_config['dna_bam'] = boost_config['input_dna_bam']
    boost_config['script_dir'] = boost_config['script_dir']
    boost_config['genome'] = boost_config['genome']
    boost_config['known_snv'] = boost_config['known_snv']
    boost_config['prefix'] = boost_config['prefix']
    boost_config['dna_name'] = boost_config['dna_name']
    
    known_vcf_from_boost = boost.run_boost_pipeline(boost_config)
    
    if not known_vcf_from_boost:
        logging.error("Boost pipeline failed. Aborting the main workflow.")
        sys.exit(1)
    logging.info(f"Boost pipeline completed. Generated known sites VCF: {known_vcf_from_boost}")


    # --- Phase 2: Run DVR Pipeline ---
    logging.info(">>> STARTING PHASE 2: DVR PIPELINE <<<")
    
    # Prepare config for the DVR pipeline
    dvr_config = pipeline_config.copy()
    dvr_config['known_vcf_from_boost'] = known_vcf_from_boost
    dvr_config['dna_bam'] = os.path.join(
        dvr_config['dna_bqsr_dir'], 
        f"{dvr_config['dna_name']}_recalibration.bam"
    ) # Use the recalibrated DNA bam
    dvr_config['rna_bam_sample1'] = dvr_config['group1_rna_bam_files']
    dvr_config['rna_bam_sample2'] = dvr_config['group2_rna_bam_files']
    dvr_config['rna_bam_sample1_str'] = ','.join(dvr_config['rna_bam_sample1'])
    dvr_config['rna_bam_sample2_str'] = ','.join(dvr_config['rna_bam_sample2'])
    dvr_config['script_dir'] = dvr_config['script_dir']
    dvr_config['genome'] = dvr_config['genome']
    dvr_config['genome_ad'] = dvr_config['genome_ad']
    dvr_config['boost_working_dir'] = dvr_config['boost_working_dir']
    dvr_config['rna_bqsr_dir'] = dvr_config['rna_bqsr_dir']
    dvr_config['dna_bqsr_dir'] = dvr_config['dna_bqsr_dir']
    dvr_config['rv_working_dir'] = dvr_config['rv_working_dir']
    dvr_config['known_snv_anno'] = dvr_config['known_snv_anno']
    dvr_config['known_editing'] = dvr_config['known_editing']
    dvr_config['gene_anno'] = dvr_config['gene_anno']
    dvr_config['prefix'] = dvr_config['prefix']
    dvr_config['threads'] = dvr_config['threads']
    dvr_config['normal_sample_name'] = dvr_config['normal_sample_name']
    dvr_config['label1'] = dvr_config['label1']
    dvr_config['label2'] = dvr_config['label2']
    dvr_config['sample_list'] = dvr_config['sample_list']

    success = dvr.run_dvr_pipeline(dvr_config)

    if not success:
        logging.error("DVR pipeline failed.")
        sys.exit(1)

    logging.info(">>> CADRES Main Workflow Completed Successfully! <<<")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Master runner for the CADRES pipeline.")
    parser.add_argument(
        "-c", "--config",
        required=True,
        help="Path to the master JSON configuration file."
    )
    args = parser.parse_args()
    
    main_workflow(args.config)
