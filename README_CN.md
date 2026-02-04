# CADRES：RNA编辑位点发现和差异分析流程

## 1. 简介

校准差异RNA编辑扫描器（CADRES）是一个全面的生物信息学流程，旨在从RNA测序（RNA-Seq）和全基因组测序（WGS）数据中精确识别差异RNA编辑位点。检测这些被称为RNA上的差异变异（DVRs）的修饰具有挑战性，因为会受到遗传变异（SNVs）和测序误差的干扰。虽然已经发现了数百万个腺苷到肌苷（A>I）的编辑位点，但在识别胞苷到尿苷（C>U）的RNA编辑位点方面存在显著差距。这主要是因为负责的酶（如APOBEC3B（A3B）等胞苷脱氨酶）可以同时编辑DNA和RNA，使得难以区分真正的RNA编辑和DNA突变。

为了解决这个问题，CADRES集成了复杂的联合DNA/RNA变异检测和严格的统计分析，以准确检测RNA编辑事件，特别关注提高C>U位点的识别。该流程的有效性已在A3B脱氨酶的可诱导细胞模型中得到验证，通过成功过滤测序伪影和A3B介导的DNA突变，证明了比现有方法更高的准确性和特异性。

本仓库包含原始CADRES流程的重构版本。保留了核心逻辑，同时将执行从单一的shell脚本现代化为模块化的Python工作流。

### 核心特性
- **高精度**：采用双阶段分析——RNA-DNA差异（RDD）和RNA-RNA差异（RRD）——以可靠地区分RNA编辑和基因组变异。
- **Boost重新校准**：实施创新的"boost重新校准"策略。这包括初始的联合DNA-RNA变异检测，以创建高置信度RNA编辑位点的新库，然后用于改进碱基质量分数重新校准（BQSR）。这一独特步骤最大限度地减少了对真实RNA变异的敏感性损失。
- **统计严谨性**：利用rMATS统计框架中强大的广义线性混合模型（GLMM）来准确量化和测试实验条件之间的差异编辑水平。
- **模块化设计**：流程分为三个主要步骤（校准、变异检测、统计测试），允许灵活执行和更简单的调试。
- **健壮执行**：基于Python的脚本，具有集成的日志记录和错误处理。
- **可重复性**：提供`environment.yml`文件以创建一致的Conda环境。

---

## 2. 安装与设置

### 前置要求
- 必须安装 [Conda](https://docs.conda.io/en/latest/miniconda.html) 包管理器。

### 环境设置

1.  **克隆仓库或下载文件**到专用项目目录。

2.  **创建Conda环境**：在终端中导航到项目目录，使用提供的`environment.yml`文件创建隔离环境。
    ```bash
    conda env create -f environment.yml
    ```

3.  **激活环境**：
    ```bash
    conda activate CADRES
    ```

---

## 3. 使用方法

该流程设计为三个顺序步骤运行。您可以单独执行每个Python脚本，或使用提供的`run.sh`脚本作为模板来编排工作流。

### 选项 1：使用封装脚本（推荐）

`run.sh`脚本提供了如何运行流程的完整示例。您应该编辑此文件以设置您的特定文件路径、样本名称和配置。

1.  在文本编辑器中打开 `run.sh`。
2.  更新 `Configuration` 部分的路径（基因组、BAM文件、输出目录）。
3.  更新 `SAMPLES` 数组以包含您的RNA样本名称。
4.  运行脚本：
    ```bash
    bash run.sh
    ```

### 选项 2：分步运行

您可以手动运行每个步骤的Python脚本。这使您可以更好地控制参数和执行。

#### 步骤 1：校准与Boost
此步骤预处理DNA和RNA BAM文件，执行初始的"Boost"变异检测以识别高置信度编辑位点，并应用碱基质量分数重新校准（BQSR）。

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

#### 步骤 2：变异检测
此步骤计算污染概况，并使用Mutect2对步骤1中的重新校准BAM文件执行最终的联合变异检测。

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

#### 步骤 3：统计测试与注释
此步骤执行统计测试（使用GLMM/MATS）以识别差异编辑位点并注释结果。

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

## 4. 输入文件

CADRES流程需要几个参考和数据文件。

### 必需的参考文件

| 参数 | 文件类型 | 说明 |
|-----------|-----------|-------------|
| `genome` | FASTA | 参考基因组序列（如GRCh38）。 |
| `known_snv` | VCF | 已知SNV数据库（如dbSNP），用于BQSR和过滤。 |
| `gnomad` | VCF | 胚系资源（如gnomAD），用于Mutect2污染估计。 |
| `known_editing` | TXT | 已知RNA编辑位点数据库（如REDIportal），用于注释。 |
| `gene_anno` | TXT | 基因注释文件（如refGene.txt），用于功能注释。 |

### 测序数据

- **DNA BAM**：全基因组测序（WGS）比对文件。必须排序并建立索引。
- **RNA BAMs**：RNA-seq比对文件。必须排序并建立索引。
