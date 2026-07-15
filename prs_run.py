import os
from pathlib import Path
import subprocess
from configparser import ConfigParser
import gzip
import shutil
import pandas as pd
import argparse
import pysam
from collections import defaultdict

def executa(comando):
    print(f"{comando}\n")
    subprocess.run(comando, shell=True, check=True)

script_dir = os.getcwd()
config_path = os.path.join(script_dir, "config.ini")
config = ConfigParser()
config.read(config_path)

project_name = config.get("Project", "PROJECT_NAME")
output_dir = config.get("Paths", "OUTPUT_DIR")
sumstats = config.get("Paths", "SUMSTATS_FILE_ANCESTRY").split(",")
run_prs = config.get("Paths", "PRSICE")
print(run_prs)

ld_pop = config.get("PRS", "LD_POP")
pheno_file = config.get("PRS", "PHENO_FILE")
binary_target = config.get("PRS", "BINARY_TARGET")
prs_missing = config.get("PRS", "PRS_MISSING")
score_prs= config.get("PRS", "SCORE_PRS")

if ld_pop != 'None':
    ld_pop = config.get("PRS", "LD_POP").split(",")
    for info_sumstat in sumstats:
        ancestry_name = info_sumstat.split(":")

        ancestry_sumstat= f'{output_dir}/{ancestry_name[-1]}/{project_name}_{ancestry_name[-1]}_merged_sumstats.txt'

        ancestry_file= f"{output_dir}/{ancestry_name[-1]}/{project_name}_{ancestry_name[-1]}_All_Chr_merged_filtered_harmonized"

        ld_file = next((item.split(':')[0] for item in ld_pop if item.split(':')[1] == ancestry_name[-1]),None)
        
        if ld_file is None:
            raise ValueError(f"LD path não encontrado para pop {ancestry_name[-1]}")
        
        prsice= f'{run_prs}/PRSice_linux  \
            --base {ancestry_sumstat} \
            --chr chr \
            --A1 effect_allele \
            --A2 other_allele \
            --stat beta \
            --snp rsid \
            --bp pos \
            --pvalue p-value \
            --beta \
            --target {ancestry_file} \
            --print-snp --out {output_dir}/{ancestry_name[-1]}/{project_name}_{ancestry_name[-1]}_PRS_run \
            --ld {ld_file} \
            --pheno  {pheno_file} \
            --score {score_prs} \
            --missing {prs_missing}\
            --pheno-col status \
            --binary-target {binary_target}'
        executa(prsice)


else:
    for info_sumstat in sumstats:
        ancestry_name = info_sumstat.split(":")
        

        ancestry_sumstat= f'{output_dir}/{ancestry_name[-1]}/{project_name}_{ancestry_name[-1]}_merged_sumstats.txt'

        ancestry_file= f"{output_dir}/{ancestry_name[-1]}/{project_name}_{ancestry_name[-1]}_All_Chr_merged_filtered_harmonized"
        
        prsice= f'{run_prs}/PRSice_linux  \
            --base {ancestry_sumstat} \
            --chr chr \
            --A1 effect_allele \
            --A2 other_allele \
            --stat beta \
            --snp rsid \
            --bp pos \
            --pvalue p-value \
            --beta \
            --target {ancestry_file} \
            --print-snp --out {output_dir}/{ancestry_name[-1]}/{project_name}_{ancestry_name[-1]}_PRS_run \
            --pheno  {pheno_file} \
            --score {score_prs} \
            --pheno-col status \
            --missing {prs_missing} \
            --binary-target {binary_target}'

        executa(prsice)
