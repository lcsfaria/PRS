import os
from pathlib import Path
import subprocess
from configparser import ConfigParser
import gzip
import shutil
import pandas as pd
import argparse

script_dir = os.getcwd()
config_path = os.path.join(script_dir, "config.ini")
config = ConfigParser()
config.read(config_path)

project_name = config.get("Project", "PROJECT_NAME")

#ancestry_splitter = config.get("Scripts", "ANCESTRY_SPLITTER")
ancestry_splitter = os.path.join(script_dir, "vcf_maker.py")


#filtrar_variantes = config.get("Scripts", "FILTRAR_VARIANTES")
filtrar_variantes = os.path.join(script_dir, "filtrar_variantes.py")

base_vcf_file = config.get("Paths", "VCF_FILE")
msp_file = config.get("Paths", "MSP_FILE")
output_dir = config.get("Paths", "OUTPUT_DIR")
sumstats = config.get("Paths", "SUMSTATS_FILE_ANCESTRY").split(",")

temp_dir = os.path.join(output_dir, "temp")

os.makedirs(output_dir, exist_ok=True)
os.makedirs(temp_dir, exist_ok=True)

#aqui junta arquivos de sumstats e VCFs por ancestralidade
for info_sumstat in sumstats:
    ancestry_name = info_sumstat.split(":")

    vcf_files = sorted(Path(f'{output_dir}/{ancestry_name[-1]}/').glob("*harmonized*.vcf.gz"))
    sumstats_files = sorted(Path(f'{output_dir}/{ancestry_name[-1]}/').glob("*filtered_sumstats.txt"))
    #print(vcf_files)
    #print(sumstats_files)

    dfs = []
    for f in sumstats_files:
        print(f"  -> {f.name}")
        dfs.append(pd.read_csv(f, sep="\t"))

    merged_sumstats = pd.concat(dfs, ignore_index=True)
    merged_sumstats_path = f'{output_dir}/{ancestry_name[-1]}/{project_name}_{ancestry_name[-1]}_merged_sumstats.txt'
    merged_sumstats.to_csv(merged_sumstats_path, sep="\t", index=False)

    vcf_list_file = f"{output_dir}/{ancestry_name[-1]}/{project_name}_{ancestry_name[-1]}_vcf_list.txt"
    output_vcf = f"{output_dir}/{ancestry_name[-1]}/{project_name}_{ancestry_name[-1]}_All_Chr_merged.vcf.gz"

    with open(vcf_list_file, "w") as f:
        for vcf in vcf_files:
            f.write(str(vcf) + "\n")

    print("Concatenando VCFs\n")
    print(f"bcftools concat -f {vcf_list_file} -Oz -o {output_vcf}\n")
    subprocess.run(["bcftools", "concat", "-f", vcf_list_file, "-Oz", "-o", output_vcf], check=True)

    print("Indexando VCF final...")
    print(f"bcftools index -t {output_vcf}\n")
    subprocess.run(["bcftools", "index", "-t", output_vcf], check=True)

    print(f"Convertendo VCF final para plink bfile\n")
    print(f'plink --vcf {output_vcf} --make-bed --vcf-half-call r --double-id --out {output_dir}/{ancestry_name[-1]}/{project_name}_{ancestry_name[-1]}_All_Chr_merged\n')
    subprocess.run(['plink', '--vcf', output_vcf, '--make-bed', '--vcf-half-call', 'r',
        '--double-id', '--out', f'{output_dir}/{ancestry_name[-1]}/{project_name}_{ancestry_name[-1]}_All_Chr_merged_filtered_harmonized'], check=True)
    