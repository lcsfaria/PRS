
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
ancestry_splitter = os.path.join(script_dir, "vcf_maker.py")
filtrar_variantes = os.path.join(script_dir, "filtrar_variantes.py")
output_dir = config.get("Paths", "OUTPUT_DIR")
sumstats = config.get("Paths", "SUMSTATS_FILE_ANCESTRY").split(",")



for info_sumstat in sumstats:
    ancestry_label = info_sumstat.split(":")[-1]
    prs_chr_pos_used = pd.read_csv(f"{output_dir}/{ancestry_label}/{project_name}_{ancestry_label}_PRS_run.snp", sep="\t")

    region_file = f"{output_dir}/{ancestry_label}/{project_name}_{ancestry_label}_PRS_chr_pos_list.txt"
    with open(region_file, "w") as fh:
        for _, row in prs_chr_pos_used.iterrows():
            chrom = row["CHR"]
            pos = int(row["BP"])
            fh.write(f"{chrom}\t{pos}\t{pos}\n")
    print(f"Extracting regions and saving in em: {region_file}\n")
    print(f'bcftools view -R {region_file} -Oz -o {output_dir}/{ancestry_label}/{project_name}_{ancestry_label}_All_Chr_merged_PRS_chr_pos.vcf.gz {output_dir}/{ancestry_label}/{ancestry_label}_All_Chr_merged.vcf.gz\n')
    subprocess.run(['bcftools', 'view', '-R', region_file,'-Oz', '-o', f"{output_dir}/{ancestry_label}/{project_name}_{ancestry_label}_All_Chr_merged_PRS_chr_pos.vcf.gz",
        f"{output_dir}/{ancestry_label}/{project_name}_{ancestry_label}_All_Chr_merged.vcf.gz"
    ], check=True)

    print(f'bcftools index -t {output_dir}/{ancestry_label}/{project_name}_{ancestry_label}_All_Chr_merged_PRS_chr_pos.vcf.gz\n')
    subprocess.run(['bcftools', 'index', '-t',
        f"{output_dir}/{ancestry_label}/{project_name}_{ancestry_label}_All_Chr_merged_PRS_chr_pos.vcf.gz"
    ], check=True)



for info_sumstat in sumstats:
    ancestry_label = info_sumstat.split(":")[-1]
    vcf_harmonized = f"{output_dir}/{ancestry_label}/{project_name}_{ancestry_label}_All_Chr_merged_PRS_chr_pos.vcf.gz"

    print(f"\nCounting missing alleles per sample to standardize ({ancestry_label})\n")
    vcf = pysam.VariantFile(vcf_harmonized)

    missing_alleles = defaultdict(int)
    present_alleles = defaultdict(int)

    for record in vcf:
        for sample in record.samples:
            gt = record.samples[sample]["GT"]

            if gt is None:
                continue

            missing = sum(a is None for a in gt)
            present = sum(a is not None for a in gt)

            missing_alleles[sample] += missing
            present_alleles[sample] += present

    df = pd.DataFrame({
        "Sample": list(missing_alleles.keys()),
        "Missing_Alleles": [missing_alleles[s] for s in missing_alleles],
        "Present_Alleles": [present_alleles[s] for s in present_alleles],
        "Ancestry": ancestry_label,
    })

    df["Total_Alleles"] = df["Missing_Alleles"] + df["Present_Alleles"]
    df["Missing_%"] = 100 * df["Missing_Alleles"] / df["Total_Alleles"]
    df["Present_%"] = 100 * df["Present_Alleles"] / df["Total_Alleles"]
    
    out_path = os.path.join(output_dir, ancestry_label, f"{project_name}_{ancestry_label}_missing_alleles_per_sample.tsv")
    print(f"Saving missing allele summary to: {out_path}")
    
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    
    
    
    df = df.rename(columns={"Sample": "IID"})
    df.to_csv(out_path, sep="\t", index=False)

    print(f"Saved missing allele summary to: {out_path}")
    

    prs_df=pd.read_csv(f"{output_dir}/{ancestry_label}/{project_name}_{ancestry_label}_PRS_run.best", sep="\s+")
    merged_df = pd.merge(prs_df, df, on="IID", how="inner")
    merged_df['standardized_PRS'] = merged_df['PRS'] / merged_df['Present_Alleles']
    merged_df.to_csv(f"{output_dir}/{ancestry_label}/{project_name}_{ancestry_label}_PRS_run.best.standardized", sep="\t", index=False)
