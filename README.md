## Installation
Clone the repository and create the virtual environment using `conda` to ensure all necessary dependencies are installed correctly:
```bash
conda env create -f environment.yml
conda activate prs_env
```
---
## How to Use
The workflow is divided into three main steps: Quality Control (QC) of the SumStats, project parameter configuration, and running the orchestrator.

### Step 1: Standardizing Summary Statistics (`qc_sumstat.py`)
Before running the main analysis, the Summary Statistics files need to be standardized.

**Basic usage example:**
```bash
python qc_sumstat.py --input gwas_raw.txt --output gwas_qced.txt --snp-col rsID --effect-col A1 --other-col A2 --chr-col CHR --bp-col BP --beta-col BETA
```

**Available Parameters:**

*Required:*
* `--input`: Input file.
* `--output`: Output file.
* `--snp-col`: SNP column name (rsID).
* `--effect-col`: Effect allele column name.
* `--other-col`: Other allele column name (reference allele).
* `--chr-col`: Chromosome column name.
* `--bp-col`: Position (basepair) column name.
* `--beta-col`: Beta column name (effect size).

*Optional:*
* `--or-col`: Odds ratio column name.
* `--se-col`: Standard error column name.
* `--p-col`: P-value column name.
* `--maf-col`: MAF/frequency column name.
* `--info-col`: Imputation INFO column name.
* `--n-col`: Sample size column name.
* `--maf`: MAF filter (Default: `0.01`).
* `--info`: INFO filter (Default: `0.8`).

### Step 2: Configuration (`config.ini`)
Fill in the `config.ini` file with your data paths and column definitions.

> **Important Notes:**
> * In the `VCF_FILE` and `MSP_FILE` fields, use the `@` character in place of the chromosome number. The orchestrator will handle the substitution automatically.
> * In the `SUMSTATS_FILE_ANCESTRY` field, the mapped population name (e.g., `AFR`, `EUR`) **must** be identical to what appears in the LAI files' header. The syntax is `file_path:POP`, separated by commas for multiple ancestries.

**Example `config.ini`:**
```ini
[Project]
PROJECT_NAME = Epigen
[Paths]
VCF_FILE = /media/lucasf/NAS/PRS_GP2/Database/GeneticData/Phased_vcf/Phased_1KGP_Peru_Chr@.vcf
MSP_FILE = /media/lucasf/NAS/PRS_GP2/Database/LA/@_query_results.msp
# The population name must match the LAI header
SUMSTATS_FILE_ANCESTRY = /media/lucasf/NAS/PRS_GP2/Sumstat/Rizig_et_al_2023_AFR_AAC_metaGWAS_no23andMe_hg38_QCed.txt:AFR,/media/lucasf/NAS/PRS_GP2/Sumstat/GP2_EUR_ONLY_HG38_12162024_rsID_QCed.tsv:EUR
OUTPUT_DIR = /home/lucasf/PRS-main/Resultado
[Columns]
SNP_COL = rsid
```

### Step 3: Running the Orchestrator (`orquestrador.py`)
With the data cleaned and the `.ini` configured, start the processing by defining which chromosomes you want to analyze.

**To process specific chromosomes (space-separated):**
```bash
python orquestrador.py --chroms 1 3 21 22
```

**To process a continuous range of chromosomes (e.g., from 1 to 22):**
```bash
python orquestrador.py --chrom-range 1 22
```

### Step 4: Running Merge Harmonized (`merge_harmonized.py`)
```bash
python merge_harmonized.py
```

### Step 5: Running PRSice2 (`prs_run.py`)
```bash
python prs_run.py
```

### Step 6: Running Standardize (`standardize_prs.py`)
```bash
python standardize_prs.py
```

### Step 7: Running Regression (`regression_multi_prs.py`)
```bash
python regression_multi_prs.py
```
