#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e
# Print commands and their arguments as they are executed
set -x

# Initialize Conda (Adjust path if necessary, based on provided scripts)
source $(conda info --base)/etc/profile.d/conda.sh
conda activate CADRES

# ==========================================
# Configuration (Based on 293_boost_V6.sh & DVR-293TV6.sh)
# ==========================================

# Base Directories
SCRIPT_DIR="/home/sunjun/cadres2/src"
OUTPUT_ROOT="/home/sunjun/cadres2/pipeline_output"

# Reference Files
GENOME="/home/sunjun/hdd/sunjun_file/genome/plsmid/Homo_sapiens.GRCh38.dna.primary_assembly_plsmid.fa"
KNOWN_SNV="/home/sunjun/hdd/sunjun_file/genome/dbsnp150.vcf"
# Using the full hg38 gnomAD file as in 293_boost_V6.sh, assuming chr22 data needs full or specific chr resource
# DVR-293TV6.sh used 'af-only-gnomad_chr_3.hg38.vcf' but samples are chr22. Using the safer full version.
GNOMAD="/home/sunjun/hdd/sunjun_file/genome/germline_resource/af-only-gnomad_chr_3.hg38.vcf"
GENE_ANNO="/home/sunjun/hdd/sunjun_file/genome/DVR_ref/refGene.txt"
KNOWN_EDITING="/home/sunjun/hdd/sunjun_file/genome/DVR_ref/REDIportal_AtoI.txt"

# Input Data
RNA_BAM_DIR="/home/sunjun/cadres2/2_rna_bam_22"
DNA_BAM_INPUT="/home/sunjun/cadres2/1_dna_bam/test_sorted_chr22.bam"

# Sample Configuration
PREFIX="293T_V6"
THREADS=16

# Sample Names (Order matters for grouping in Step 3)
# Group 1 (DMSO / Control): C6-DOX
# Group 2 (DOX / Treated): C6+DOX
SAMPLES=(
    "C6-DOX-1A_chr22"
    "C6-DOX-2A_chr22"
    "C6-DOX-3A_chr22"
    "C6+DOX-1A_chr22"
    "C6+DOX-2A_chr22"
    "C6+DOX-3A_chr22"
)

# Create Output Directory
mkdir -p "$OUTPUT_ROOT"

# ==========================================
# Preparation
# ==========================================

# Construct Input RNA BAM List
RNA_BAMS_INPUT=""
for sample in "${SAMPLES[@]}"; do
    RNA_BAMS_INPUT="$RNA_BAMS_INPUT $RNA_BAM_DIR/${sample}.bam"
done

# ==========================================
# Step 1: Calibration
# ==========================================
echo "========================================"
echo "Starting Step 1: Calibration"
echo "========================================"

STEP1_OUT="$OUTPUT_ROOT/step1_calibration"
<<EOF
python "$SCRIPT_DIR/pipeline_step1_calibration.py" \
    --rna_bams $RNA_BAMS_INPUT \
    --dna_bam "$DNA_BAM_INPUT" \
    --genome "$GENOME" \
    --known_snv "$KNOWN_SNV" \
    --output_dir "$STEP1_OUT" \
    --prefix "$PREFIX" \
    --threads "$THREADS" \
    --dna_sample_name "$DNA_SAMPLE_NAME"
EOF
# ==========================================
# Step 2: Variant Calling
# ==========================================
echo "========================================"
echo "Starting Step 2: Variant Calling"
echo "========================================"

STEP2_OUT="$OUTPUT_ROOT/step2_variant_calling"

# Identify Recalibrated BAMs from Step 1
# DNA: The script produces DNA_processed_recalibration.bam if BQSR runs, or dedup.bam otherwise.
# Assuming BQSR runs since known_snv is provided.
DNA_RECAL="$STEP1_OUT/DNA_processed_recalibration.bam"
if [ ! -f "$DNA_RECAL" ]; then
    echo "Warning: $DNA_RECAL not found, checking for dedup.bam..."
    DNA_RECAL="$STEP1_OUT/DNA_processed_dedup.bam"
fi

# RNA: The script produces {sample}_recalibration.bam
RNA_RECAL_BAMS=""
for sample in "${SAMPLES[@]}"; do
    RNA_RECAL_BAMS="$RNA_RECAL_BAMS $STEP1_OUT/${sample}_recalibration.bam"
done

python "$SCRIPT_DIR/pipeline_step2_variant_calling.py" \
    --rna_bams $RNA_RECAL_BAMS \
    --dna_bam "$DNA_RECAL" \
    --genome "$GENOME" \
    --gnomad "$GNOMAD" \
    --output_dir "$STEP2_OUT" \
    --prefix "$PREFIX" \
    --scripts_dir "$SCRIPT_DIR" \
    --threads "$THREADS"

# ==========================================
# Step 3: Statistical Test
# ==========================================
echo "========================================"
echo "Starting Step 3: Statistical Test"
echo "========================================"

STEP3_OUT="$OUTPUT_ROOT/step3_statistical_test"
FINAL_VCF="$STEP2_OUT/$PREFIX.final.vcf"

# Group Indices (0-based)
# Group 1: Indices 0, 1, 2 (C6-DOX)
# Group 2: Indices 3, 4, 5 (C6+DOX)

GROUP1_BAMS=""
for sample in "${SAMPLES[@]:0:3}"; do
    GROUP1_BAMS="$GROUP1_BAMS $STEP1_OUT/${sample}_recalibration.bam"
done

GROUP2_BAMS=""
for sample in "${SAMPLES[@]:3:3}"; do
    GROUP2_BAMS="$GROUP2_BAMS $STEP1_OUT/${sample}_recalibration.bam"
done

python "$SCRIPT_DIR/pipeline_step3_statistical_test.py" \
    --group1_rna_bams $GROUP1_BAMS \
    --group2_rna_bams $GROUP2_BAMS \
    --final_vcf "$FINAL_VCF" \
    --genome "$GENOME" \
    --output_dir "$STEP3_OUT" \
    --prefix "$PREFIX" \
    --scripts_dir "$SCRIPT_DIR" \
    --labels DMSO DOX \
    --known_snv "$KNOWN_SNV" \
    --known_editing "$KNOWN_EDITING" \
    --gene_anno "$GENE_ANNO" \
    --threads "$THREADS"

echo "========================================"
echo "Pipeline Completed Successfully!"
echo "Outputs are in: $OUTPUT_ROOT"
echo "========================================"
