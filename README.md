## Instalação

Clone o repositório e crie o ambiente virtual utilizando o `conda` para garantir que todas as dependências necessárias sejam instaladas corretamente:

```bash
conda env create -f environment.yml
conda activate prs_env
```

---

## Como Usar

O fluxo de trabalho é dividido em três etapas principais: Controle de Qualidade (QC) das SumStats, configuração dos parâmetros do projeto e execução do orquestrador.

### Passo 1: Padronização das Summary Statistics (`qc_sumstat.py`)

Antes de rodar a análise principal, é necessário padronizar os arquivos de Summary Statistics.

**Exemplo de uso básico:**
```bash
python qc_sumstat.py --input gwas_raw.txt --output gwas_qced.txt --snp-col rsID --effect-col A1 --other-col A2 --chr-col CHR --bp-col BP --beta-col BETA
```

**Parâmetros Disponíveis:**

*Obrigatórios:*
* `--input`: Arquivo de entrada.
* `--output`: Arquivo de saída.
* `--snp-col`: Nome da coluna de SNP (rsID).
* `--effect-col`: Nome da coluna do alelo de efeito.
* `--other-col`: Nome da coluna do outro alelo (alelo de referência).
* `--chr-col`: Nome da coluna de cromossomo.
* `--bp-col`: Nome da coluna de posição (basepair).
* `--beta-col`: Nome da coluna de beta (tamanho do efeito).

*Opcionais:*
* `--or-col`: Nome da coluna de odds ratio.
* `--se-col`: Nome da coluna de erro padrão.
* `--p-col`: Nome da coluna de p-valor.
* `--maf-col`: Nome da coluna de MAF/frequência.
* `--info-col`: Nome da coluna de INFO de imputação.
* `--n-col`: Nome da coluna de tamanho de amostra.
* `--maf`: Filtro de MAF (Valor padrão: `0.01`).
* `--info`: Filtro de INFO (Valor padrão: `0.8`).

### Passo 2: Configuração (`config.ini`)

Complete o arquivo `config.ini` com os caminhos dos seus dados e definições de colunas. 

> **Avisos Importantes:**
> * Nos campos `VCF_FILE` e `MSP_FILE`, utilize o caractere `@` no lugar do número do cromossomo. O orquestrador fará a substituição automaticamente.
> * No campo `SUMSTATS_FILE_ANCESTRY`, o nome da população mapeada (ex: `AFR`, `EUR`) **precisa** ser idêntico ao que consta no cabeçalho dos arquivos de LAI. A sintaxe é `caminho_do_arquivo:POP`, separados por vírgula para múltiplas ancestralidades.

**Exemplo de `config.ini`:**
```ini
[Project]
PROJECT_NAME = Epigen

[Paths]
VCF_FILE = /media/lucasf/NAS/PRS_GP2/Database/GeneticData/Phased_vcf/Phased_1KGP_Peru_Chr@.vcf
MSP_FILE = /media/lucasf/NAS/PRS_GP2/Database/LA/@_query_results.msp

# O nome da pop tem que estar igual ao do header LAI
SUMSTATS_FILE_ANCESTRY = /media/lucasf/NAS/PRS_GP2/Sumstat/Rizig_et_al_2023_AFR_AAC_metaGWAS_no23andMe_hg38_QCed.txt:AFR,/media/lucasf/NAS/PRS_GP2/Sumstat/GP2_EUR_ONLY_HG38_12162024_rsID_QCed.tsv:EUR
OUTPUT_DIR = /home/lucasf/PRS-main/Resultado

[Columns]
SNP_COL = rsid
```

### Passo 3: Execução do Orquestrador (`orquestrador.py`)

Com os dados limpos e o `.ini` configurado, inicie o processamento definindo quais cromossomos você deseja analisar. 

**Para processar cromossomos específicos (separados por espaço):**
```bash
python orquestrador.py --chroms 1 3 21 22
```

**Para processar um intervalo contínuo de cromossomos (ex: do 1 ao 22):**
```bash
python orquestrador.py --chrom-range 1 22
```

### Passo 4: Execução Merged Hamonized (`merge_harmonized.py`)

```bash
python merge_harmonized.py
```

### Passo 5: Execução PRSice2 (`prs_run.py`)

```bash
python prs_run.py
```

### Passo 6: Execução Standarize (`standardize_prs.py`)

```bash
python standardize_prs.py
```
