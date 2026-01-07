# CADRES: RNA Editing Site Discovery and Differential Analysis Pipeline

## 1. Introduction

The Calibrated Differential RNA Editing Scanner (CADRES) is a comprehensive bioinformatics pipeline designed to precisely identify differential RNA editing sites from both RNA-Seq and whole-genome sequencing (WGS) data. Detecting these modifications, known as Differential Variants on RNA (DVRs), is challenging due to interference from genetic variants (SNVs) and sequencing errors. While millions of Adenosine-to-Inosine (A>I) editing sites have been found, there is a significant gap in identifying Cytidine-to-Uridine (C>U) RNA editing sites. This is largely because the enzymes responsible, cytidine deaminases like APOBEC3B (A3B), can edit both DNA and RNA, making it difficult to distinguish true RNA edits from DNA mutations.
To address this, CADRES integrates sophisticated joint DNA/RNA variant calling with rigorous statistical analysis to accurately detect RNA editing events, with a particular focus on improving the identification of C>U sites. The pipeline's effectiveness has been validated using inducible cell models of the A3B deaminase, where it demonstrated improved accuracy and specificity over existing methods by successfully filtering out sequencing artifacts and A3B-mediated DNA mutations.
This repository contains a refactored version of the original CADRES pipeline. The core logic is preserved, while the execution has been modernized from shell scripts into a modular, robust, and configurable Python workflow.
Key Features
High Precision: Employs a dual-phase analysis—RNA-DNA Difference (RDD) and RNA-RNA Difference (RRD)—to reliably distinguish RNA edits from genomic variants. This stringent approach significantly increases the precision of detected editing events.
Boost Recalibration: Implements an innovative 'boost recalibration' strategy. This involves an initial joint DNA-RNA variant call to create a de novo library of high-confidence RNA editing sites, which is then used to refine Base Quality Score Recalibration (BQSR). This unique step minimizes the loss of sensitivity for true RNA variants, which can be mistaken for sequencing artifacts by standard BQSR methods.
Statistical Rigor: Utilizes the robust Generalized Linear Mixed Model (GLMM) from the rMATS statistical framework to accurately quantify and test for differential editing levels between experimental conditions.
Modular Design: The pipeline is split into two main phases (Boost and DVR), each composed of functional steps for clarity and maintainability.
Centralized Configuration: All paths, parameters, and sample information are managed in a single config.json file, eliminating the need to edit scripts.
Robust Execution: The workflow is managed by a master Python script (main.py) with integrated logging and error handling.
Reproducibility: A environment.yml file is provided to create a consistent Conda environment with all necessary dependencies.

### Key Features
- **Modular Design**: The pipeline is split into two main phases (Boost and DVR), each composed of functional steps.
- **Centralized Configuration**: All paths, parameters, and sample information are managed in a single `config.json` file, eliminating the need to edit scripts.
- **Robust Execution**: The workflow is managed by a master Python script (`main.py`) with integrated logging and error handling.
- **Reproducibility**: A `environment.yml` file is provided to create a consistent Conda environment with all necessary dependencies.

---

## 2. Installation & Setup

### Prerequisites
- [Conda](https://docs.conda.io/en/latest/miniconda.html) package manager must be installed.

### Environment Setup
1.  **Clone the repository or download the files** into a dedicated project directory.

2.  **Create the Conda Environment**: Navigate to the project directory in your terminal and use the provided `environment.yml` file to create the isolated environment. This will install all required bioinformatics tools and Python libraries.
    ```bash
    conda env create -f environment.yml
    ```

3.  **Activate the Environment**: Before running any part of the pipeline, you must activate the newly created environment.
    ```bash
    conda activate CADRES
    ```

---

## 3. Configuration

The entire pipeline is controlled by the `config.json` file. Before running the analysis, you must edit this file to match your file paths, working directories, and sample details.

#### `config.json` Structure:

-   **`general`**: Global settings like project prefix, number of threads, and log file name.
-   **`paths`**: Absolute paths to all input data, including the reference genome, dbSNP files, annotation files, and raw BAM files.
-   **`work_dirs`**: Paths to the directories where intermediate and final results for each major step will be stored. These directories will be created automatically if they do not exist.
-   **`samples`**: Information about your samples.
    -   `dna_name`: The sample name of the WGS/DNA data.
    -   `all_rna_samples`: A list of all RNA sample names (base names, without `.bam`).
    -   `group1_names` & `group2_names`: Lists of RNA sample names belonging to each experimental condition.
    -   `group1_label` & `group2_label`: The labels for the two conditions (e.g., "Control", "Treatment").

---

## 4. Input Files

The CADRES pipeline requires the following input files. All file paths must be specified in the `paths` section of `config.json`.

### Required Input Files

| Configuration Parameter | File Type | Description |
|------------------------|-----------|-------------|
| `genome` | FASTA | Reference genome sequence file (e.g., GRCh38). Used for sequence alignment and variant calling. |
| `known_snv_boost` | VCF | Known single nucleotide variant (SNV) database for variant detection and filtering in the Boost phase. Typically uses dbSNP database. |
| `known_snv_anno` | VCF | Known SNV annotation file for final result variant annotation. Should match the reference genome version. |
| `genome_ad` | VCF | Allele frequency database (e.g., gnomAD) for BQSR (Base Quality Score Recalibration) step. |
| `known_editing` | TXT | Known RNA editing site database (e.g., REDIportal) for reference and validation. |
| `gene_anno` | TXT | Gene annotation file (e.g., refGene.txt) for functional annotation of detected editing sites. |

### Sequencing Data Files

| Configuration Parameter | File Type | Description |
|------------------------|-----------|-------------|
| `input_dna_bam` | BAM | Whole-genome sequencing (WGS) data alignment file. Must be a sorted and indexed BAM file (.bam and .bai). |
| `input_rna_bam_dir` | Directory | Directory containing all RNA-seq sample alignment files. Each sample's BAM file should be named `{sample_name}.bam` with a corresponding index file `.bai`. |

### Input File Requirements

1. **BAM Files**:
   - Must be sorted (using samtools sort)
   - Must have corresponding index files (.bai)
   - Recommended to process all samples with the same alignment tool and parameters

2. **Reference Genome**:
   - FASTA format
   - Recommended to use a genome version that matches the dbSNP database version (e.g., GRCh38)
   - Should contain all required chromosome sequences

3. **VCF Files**:
   - Must match the reference genome version
   - Recommended to use standardized VCF format
   - Should contain necessary INFO and FORMAT fields

4. **Sample Naming**:
   - RNA-seq BAM file names should match the sample names in `config.json`
   - Sample names without `.bam` suffix are used for configuration

---

## 5. Running the Pipeline

The workflow is executed via the `main.py` script, which orchestrates all the steps.

1.  **Verify Configuration**: Double-check that all paths and names in `config.json` are correct.

2.  **Execute the Master Script**: Run the following command from your terminal within the project directory (ensure the `CADRES` conda environment is active).
    ```bash
    python main.py --config config.json
    ```

The script will then:
-   Set up logging to both the console and the file specified in `config.json`.
-   Run the **Boost Phase**: Calibrate DNA/RNA BAMs and run Mutect2 to create a high-confidence panel of known variant sites.
-   Run the **DVR Phase**: Use the output from the boost phase to perform BQSR, contamination checks, variant calling, custom filtering, and finally, the differential editing analysis.

The progress will be printed to the console, and detailed logs will be saved for debugging and record-keeping.

---

## 6. File Structure

A recommended structure for your project directory:

```
/path/to/your_project/
├── main.py
├── run_boost_pipeline.py
├── run_dvr_pipeline.py
├── config.json
├── environment.yml
├── README.md
├── scripts/                  # Directory for all helper python scripts
├── bam_calibration_DNA.py
├── bam_calibration_RNA_boost.py
├── revised_convertVCF.py
├── filter_homopolymer_nucleotides.py
├── pblat_candidates_filter.py
└── ... (other helper .py scripts)
└── data/                     # (Optional) For storing input data
└── results/                  # (Optional) For storing top-level results
```
*Note: The `config.json` provided assumes `script_dir` is ".", so all scripts should be in the main directory. If you place them in a `scripts/` subdirectory, update the `script_dir` path in `config.json` to `"./scripts/"`.*

---

## 7. Output

The most critical output file will be generated by the final step of the DVR pipeline inside your specified `rv_working_dir`:

-   **`{prefix}_rMATS-DVR_Result.txt`**: A tab-delimited file containing the final annotated list of DVR sites, including their location, p-values, FDR, editing level differences, and gene annotations.
-   **`{prefix}rMATS-DVR_Result_summary.txt`**: A summary of the different types of variants found.

## 8. How to Cite
If you use the CADRES pipeline in your research, please cite the following publication:

Sun, J., Zhang, C. & Li, X. Precise detection of differential RNA editing sites across varied biological conditions using the CADRES pipeline. Sci Rep 15, 19683 (2025). https://doi.org/10.1038/s41598-025-04957-7