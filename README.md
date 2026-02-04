# CADRES: RNA Editing Site Discovery and Differential Analysis Pipeline

## 1. Introduction

The Calibrated Differential RNA Editing Scanner (CADRES) is a comprehensive bioinformatics pipeline designed to precisely identify differential RNA editing sites from both RNA-Seq and whole-genome sequencing (WGS) data. Detecting these modifications, known as Differential Variants on RNA (DVRs), is challenging due to interference from genetic variants (SNVs) and sequencing errors. While millions of Adenosine-to-Inosine (A>I) editing sites have been found, there is a significant gap in identifying Cytidine-to-Uridine (C>U) RNA editing sites. This is largely because the enzymes responsible, cytidine deaminases like APOBEC3B (A3B), can edit both DNA and RNA, making it difficult to distinguish true RNA edits from DNA mutations.

To address this, CADRES integrates sophisticated joint DNA/RNA variant calling with rigorous statistical analysis to accurately detect RNA editing events, with a particular focus on improving the identification of C>U sites. The pipeline's effectiveness has been validated using inducible cell models of the A3B deaminase, where it demonstrated improved accuracy and specificity over existing methods by successfully filtering out sequencing artifacts and A3B-mediated DNA mutations.

This repository contains a refactored version of the original CADRES pipeline. The core logic is preserved, while the execution has been modernized from monolithic shell scripts into a modular Python workflow.

### Key Features
- **High Precision**: Employs a dual-phase analysis—RNA-DNA Difference (RDD) and RNA-RNA Difference (RRD)—to reliably distinguish RNA edits from genomic variants.
- **Boost Recalibration**: Implements an innovative 'boost recalibration' strategy. This involves an initial joint DNA-RNA variant call to create a de novo library of high-confidence RNA editing sites, which is then used to refine Base Quality Score Recalibration (BQSR). This minimizes the loss of sensitivity for true RNA variants.
- **Statistical Rigor**: Utilizes the robust Generalized Linear Mixed Model (GLMM) from the rMATS statistical framework to accurately quantify and test for differential editing levels.
- **Modular Design**: The pipeline is split into three main steps (Calibration, Variant Calling, Statistical Testing), allowing for flexible execution and easier debugging.
- **Robust Execution**: Python-based scripts with integrated logging and error handling.
- **Reproducibility**: A `environment.yml` file is provided to create a consistent Conda environment.

---

## 2. Installation & Setup

### Prerequisites
- [Conda](https://docs.conda.io/en/latest/miniconda.html) package manager must be installed.

### Environment Setup

1.  **Clone the repository or download the files** into a dedicated project directory.

2.  **Create the Conda Environment**: Navigate to the project directory in your terminal and use the provided `environment.yml` file to create the isolated environment.
    ```bash
    conda env create -f environment.yml
    ```

3.  **Activate the Environment**:
    ```bash
    conda activate CADRES
    ```

---

## 3. Usage

The pipeline is designed to be run in three sequential steps. You can execute each python script independently or use the provided `run.sh` script as a template to orchestrate the workflow.

### Option 1: Using the Wrapper Script (Recommended)

The `run.sh` script provides a complete example of how to run the pipeline. You should edit this file to set your specific file paths, sample names, and configuration.

1.  Open `run.sh` in a text editor.
2.  Update the `Configuration` section with your paths (Genome, BAM files, Output directory).
3.  Update the `SAMPLES` array with your RNA sample names.
4.  Run the script:
    ```bash
    bash run.sh
    ```

### Option 2: Running Step-by-Step

You can run the Python scripts for each step manually. This gives you more control over the arguments and execution.

#### Step 1: Calibration & Boost
This step preprocesses the DNA and RNA BAM files, performs an initial "Boost" variant calling to identify high-confidence editing sites, and applies Base Quality Score Recalibration (BQSR).

```bash
python pipeline_step1_calibration.py \
    --rna_bams /path/to/rna1.bam /path/to/rna2.bam ... \
    --dna_bam /path/to/dna.bam \
    --genome /path/to/genome.fa \
    --known_snv /path/to/dbsnp.vcf \
    --output_dir ./output/step1 \
    --prefix my_project \
    --threads 16
```

#### Step 2: Variant Calling
This step calculates contamination profiles and performs the final joint variant calling using Mutect2 on the recalibrated BAM files from Step 1.

```bash
python pipeline_step2_variant_calling.py \
    --rna_bams ./output/step1/rna1_split_recalibration.bam ... \
    --dna_bam ./output/step1/DNA_processed_recalibration.bam \
    --genome /path/to/genome.fa \
    --gnomad /path/to/gnomad.vcf \
    --output_dir ./output/step2 \
    --prefix my_project \
    --threads 16
```

#### Step 3: Statistical Testing & Annotation
This step performs statistical testing (using GLMM/MATS) to identify differential editing sites and annotates the results.

```bash
python pipeline_step3_statistical_test.py \
    --group1_rna_bams ./output/step1/rna1_split_recalibration.bam ... \
    --group2_rna_bams ./output/step1/rna3_split_recalibration.bam ... \
    --final_vcf ./output/step2/my_project.final.vcf \
    --genome /path/to/genome.fa \
    --known_snv /path/to/dbsnp.vcf \
    --known_editing /path/to/rediportal.txt \
    --gene_anno /path/to/refGene.txt \
    --output_dir ./output/step3 \
    --labels Control Treated \
    --threads 4
```

---

## 4. Input Files

The CADRES pipeline requires several reference and data files.

### Required Reference Files

| Parameter | File Type | Description |
|-----------|-----------|-------------|
| `genome` | FASTA | Reference genome sequence (e.g., GRCh38). |
| `known_snv` | VCF | Known SNV database (e.g., dbSNP) for BQSR and filtering. |
| `gnomad` | VCF | Germline resource (e.g., gnomAD) for Mutect2 contamination estimation. |
| `known_editing` | TXT | Known RNA editing site database (e.g., REDIportal) for annotation. |
| `gene_anno` | TXT | Gene annotation file (e.g., refGene.txt) for functional annotation. |

### Sequencing Data

- **DNA BAM**: Whole-genome sequencing (WGS) alignment file. Must be sorted and indexed.
- **RNA BAMs**: RNA-seq alignment files. Must be sorted and indexed.
