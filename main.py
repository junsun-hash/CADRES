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
    
    # Specific sample config processing
    sample_config = master_config.get('samples', {})
    config['dna_name'] = sample_config.get('dna_name')
    config['sample_list'] = sample_config.get('all_rna_samples')
    config['group1_names'] = sample_config.get('group1_names')
    config['group2_names'] = sample_config.get('group2_names')
    config['label1'] = sample_config.get('group1_label')
    config['label2'] = sample_config.get('group2_label')
    config['normal_sample_name'] = sample_config.get('dna_name') # for DVR step
    
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
    boost_config['rna_bam_dir'] = boost_config['input_rna_bam_dir']
    boost_config['dna_bam'] = boost_config['input_dna_bam']
    
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
