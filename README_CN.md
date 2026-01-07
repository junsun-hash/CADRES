# CADRES：RNA编辑位点发现和差异分析流程

## 1. 简介

校准差异RNA编辑扫描器（CADRES）是一个全面的生物信息学流程，旨在从RNA测序（RNA-Seq）和全基因组测序（WGS）数据中精确识别差异RNA编辑位点。检测这些被称为RNA上的差异变异（DVRs）的修饰具有挑战性，因为会受到遗传变异（SNVs）和测序误差的干扰。虽然已经发现了数百万个腺苷到肌苷（A>I）的编辑位点，但在识别胞苷到尿苷（C>U）的RNA编辑位点方面存在显著差距。这主要是因为负责的酶（如APOBEC3B（A3B）等胞苷脱氨酶）可以同时编辑DNA和RNA，使得难以区分真正的RNA编辑和DNA突变。

为了解决这个问题，CADRES集成了复杂的联合DNA/RNA变异检测和严格的统计分析，以准确检测RNA编辑事件，特别关注提高C>U位点的识别。该流程的有效性已在A3B脱氨酶的可诱导细胞模型中得到验证，通过成功过滤测序伪影和A3B介导的DNA突变，证明了比现有方法更高的准确性和特异性。

本仓库包含原始CADRES流程的重构版本。保留了核心逻辑，同时将执行从shell脚本现代化为模块化、健壮且可配置的Python工作流。

### 核心特性

- **高精度**：采用双阶段分析——RNA-DNA差异（RDD）和RNA-RNA差异（RRD）——以可靠地区分RNA编辑和基因组变异。这种严格的方法显著提高了检测到的编辑事件的精度。
- **Boost重新校准**：实施创新的"boost重新校准"策略。这包括初始的联合DNA-RNA变异检测，以创建高置信度RNA编辑位点的新库，然后用于改进碱基质量分数重新校准（BQSR）。这一独特步骤最大限度地减少了对真实RNA变异的敏感性损失，这些变异可能被标准BQSR方法误认为是测序伪影。
- **统计严谨性**：利用rMATS统计框架中强大的广义线性混合模型（GLMM）来准确量化和测试实验条件之间的差异编辑水平。
- **模块化设计**：流程分为两个主要阶段（Boost和DVR），每个阶段由功能步骤组成，以提高清晰度和可维护性。
- **集中配置**：所有路径、参数和样本信息都在单个config.json文件中管理，无需编辑脚本。
- **健壮执行**：工作流由主Python脚本（main.py）管理，具有集成的日志记录和错误处理。
- **可重复性**：提供environment.yml文件以创建一致的Conda环境，包含所有必要的依赖项。

---

## 2. 安装与设置

### 前置要求
- 必须安装 [Conda](https://docs.conda.io/en/latest/miniconda.html) 包管理器。

### 环境设置

1.  **克隆仓库或下载文件**到专用项目目录。

2.  **创建Conda环境**：在终端中导航到项目目录，使用提供的`environment.yml`文件创建隔离环境。这将安装所有必需的生物信息学工具和Python库。
    ```bash
    conda env create -f environment.yml
    ```

3.  **激活环境**：在运行流程的任何部分之前，必须激活新创建的环境。
    ```bash
    conda activate CADRES
    ```

---

## 3. 配置

整个流程由`config.json`文件控制。在运行分析之前，必须编辑此文件以匹配您的文件路径、工作目录和样本详细信息。

#### `config.json` 结构：

-   **`general`**：全局设置，如项目前缀、线程数和日志文件名。
-   **`paths`**：所有输入数据的绝对路径，包括参考基因组、dbSNP文件、注释文件和原始BAM文件。
-   **`work_dirs`**：将存储每个主要步骤的中间和最终结果的目录路径。如果这些目录不存在，将自动创建。
-   **`samples`**：关于样本的信息。
    -   `dna_name`：WGS/DNA数据的样本名称。
    -   `all_rna_samples`：所有RNA样本名称的列表（基本名称，不含`.bam`）。
    -   `group1_names` & `group2_names`：属于每个实验条件的RNA样本名称列表。
    -   `group1_label` & `group2_label`：两个条件的标签（例如，"Control"、"Treatment"）。

---

## 4. 输入

CADRES流程需要以下输入文件。所有文件路径必须在`config.json`的`paths`部分指定。

### 必需的输入文件

| 配置参数 | 文件类型 | 说明 |
|---------|---------|------|
| `genome` | FASTA | 参考基因组序列文件（如GRCh38）。用于序列比对和变异检测。 |
| `known_snv_boost` | VCF | 已知单核苷酸变异（SNV）数据库，用于Boost阶段的变异检测和过滤。通常使用dbSNP数据库。 |
| `known_snv_anno` | VCF | 已知SNV注释文件，用于最终结果的变异注释。应与参考基因组版本匹配。 |
| `genome_ad` | VCF | 等位基因频率数据库（如gnomAD），用于BQSR（碱基质量分数重新校准）步骤。 |
| `known_editing` | TXT | 已知RNA编辑位点数据库（如REDIportal），用于参考和验证。 |
| `gene_anno` | TXT | 基因注释文件（如refGene.txt），用于对检测到的编辑位点进行基因功能注释。 |

### 测序数据文件

| 配置参数 | 文件类型 | 说明 |
|---------|---------|------|
| `input_dna_bam` | BAM | 全基因组测序（WGS）数据的比对文件，必须是排序和索引的BAM文件（.bam和.bai）。 |
| `input_rna_bam_dir` | 目录 | 包含所有RNA-seq样本比对文件的目录。每个样本的BAM文件应命名为`{sample_name}.bam`，并有对应的索引文件`.bai`。 |

### 输入文件要求

1. **BAM文件**：
   - 必须经过排序（使用samtools sort）
   - 必须有对应的索引文件（.bai）
   - 建议使用相同的比对工具和参数处理所有样本

2. **参考基因组**：
   - FASTA格式
   - 建议使用与dbSNP数据库版本匹配的基因组版本（如GRCh38）
   - 应包含所有必需的染色体序列

3. **VCF文件**：
   - 必须与参考基因组版本匹配
   - 建议使用标准化的VCF格式
   - 应包含必要的INFO和FORMAT字段

4. **样本命名**：
   - RNA-seq BAM文件名应与`config.json`中的样本名称一致
   - 不包含`.bam`后缀的样本名称用于配置

---

## 5. 运行流程

工作流通过`main.py`脚本执行，该脚本协调所有步骤。

1.  **验证配置**：再次检查`config.json`中的所有路径和名称是否正确。

2.  **执行主脚本**：在项目目录内的终端中运行以下命令（确保`CADRES` conda环境处于活动状态）。
    ```bash
    python main.py --config config.json
    ```

脚本将：
-   设置日志记录到控制台和`config.json`中指定的文件。
-   运行**Boost阶段**：校准DNA/RNA BAM文件并运行Mutect2以创建已知变异位点的高置信度面板。
-   运行**DVR阶段**：使用boost阶段的输出来执行BQSR、污染检查、变异检测、自定义过滤，最后进行差异编辑分析。

进度将打印到控制台，详细日志将保存用于调试和记录。

---

## 6. 文件结构

项目目录的推荐结构：

```
/path/to/your_project/
├── main.py
├── run_boost_pipeline.py
├── run_dvr_pipeline.py
├── config.json
├── environment.yml
├── README.md
├── scripts/                  # 所有辅助python脚本的目录
│   ├── bam_calibration_DNA.py
│   ├── bam_calibration_RNA_boost.py
│   ├── revised_convertVCF.py
│   ├── filter_homopolymer_nucleotides.py
│   ├── pblat_candidates_filter.py
│   └── ... (其他辅助 .py 脚本)
└── data/                     # （可选）用于存储输入数据
└── results/                  # （可选）用于存储顶层结果
```
*注意：提供的`config.json`假设`script_dir`是"."，因此所有脚本都应该在主目录中。如果将它们放在`scripts/`子目录中，请将`config.json`中的`script_dir`路径更新为`"./scripts/"`。*

---

## 7. 输出

最重要的输出文件将在DVR流程的最后一步在您指定的`rv_working_dir`内生成：

-   **`{prefix}_rMATS-DVR_Result.txt`**：包含最终注释的DVR位点列表的制表符分隔文件，包括它们的位置、p值、FDR、编辑水平差异和基因注释。
-   **`{prefix}rMATS-DVR_Result_summary.txt`**：发现的变异类型摘要。

## 8. 如何引用

如果您在研究中使用CADRES流程，请引用以下出版物：

Sun, J., Zhang, C. & Li, X. Precise detection of differential RNA editing sites across varied biological conditions using the CADRES pipeline. Sci Rep 15, 19683 (2025). https://doi.org/10.1038/s41598-025-04957-7
