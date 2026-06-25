#!/usr/bin/env python3

import pandas as pd
import gzip
import argparse
import subprocess
import os

'''
Como executar: 

python filtrar_variantes.py \
    --vcf dados.vcf.gz \
    --sumstats gwas.txt \
    --snp-col SNP \
    --out-vcf vcf_filtrado.vcf \
    --out-sumstats sumstats_filtrado.txt

'''

def carregar_ids_vcf(vcf_file):
    ids = set()
    opener = gzip.open if vcf_file.endswith(".gz") else open

    with opener(vcf_file, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            campos = line.strip().split("\t")
            ids.add(campos[2].strip())  # coluna ID do VCF

    return {str(x) for x in ids}


def normalize_dataframe_headers(df):
    df.columns = [str(col).strip() for col in df.columns]
    return df


def normalize_id_values(series):
    return series.astype(str).str.strip()


def filtrar_summary(summary_file, snp_col, ids_vcf, output_file):
    df = pd.read_csv(summary_file, sep=None, engine="python")
    df = normalize_dataframe_headers(df)

    if snp_col not in df.columns:
        raise ValueError(
            f"Coluna '{snp_col}' não encontrada no sumstats. Colunas disponíveis: {list(df.columns)}"
        )

    df[snp_col] = normalize_id_values(df[snp_col])
    df_filtrado = df[df[snp_col].isin(ids_vcf)].copy()
    df_filtrado.drop_duplicates(subset=[snp_col], inplace=True)

    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    df_filtrado.to_csv(output_file, sep="\t", index=False)

    return set(df_filtrado[snp_col])



def filtrar_vcf(vcf_file, ids_validos, output_file):
    opener = gzip.open if vcf_file.endswith(".gz") else open

    with opener(f'{vcf_file}', "rt") as fin, open(f'{output_file}', "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue

            campos = line.strip().split("\t")

            if campos[2] in ids_validos:
                fout.write(line)


def compress_and_index_vcf(vcf_file):
    """
    Comprime um VCF usando bcftools e cria o índice.
    """

    gz_file = f"{vcf_file}.gz"

    subprocess.run(
        ["bcftools", "view", "-Oz", "-o", gz_file, vcf_file],
        check=True
    )

    subprocess.run(
        ["bcftools", "index", "-t", gz_file],
        check=True
    )

    # REMOVIDO/COMENTADO PARA NÃO APAGAR O ARQUIVO .VCF ORIGINAL
    # os.remove(vcf_file)

    return gz_file

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--vcf", required=True)
    parser.add_argument("--sumstats", required=True)
    parser.add_argument("--snp-col", required=True, help="Nome da coluna com os SNPs no summary statistics")
    parser.add_argument("--out-vcf", required=True)
    parser.add_argument("--out-sumstats", required=True)

    args = parser.parse_args()

    # --- TRATAMENTO ROBUSTO DE EXTENSÕES ---
    if args.out_vcf.endswith('.vcf.gz'):
        args.out_vcf = args.out_vcf[:-3]
    elif not args.out_vcf.endswith('.vcf'):
        args.out_vcf += '.vcf'

    if not args.out_sumstats.endswith(('.txt', '.tsv', '.csv')):
        args.out_sumstats += '.txt'
    # ---------------------------------------

    print("Lendo IDs do VCF...")
    ids_vcf = carregar_ids_vcf(args.vcf)

    print("Filtrando summary statistics...")
    ids_comuns = filtrar_summary(
        args.sumstats,
        args.snp_col,
        ids_vcf,
        args.out_sumstats
    )

    print("Filtrando VCF...")
    filtrar_vcf(
        args.vcf,
        ids_comuns,
        args.out_vcf
    )
    
    print("Comprimindo e indexando o VCF...")
    vcf_gz = compress_and_index_vcf(args.out_vcf)
    
    print(f"\n--- Resumo dos Arquivos Salvos ---")
    print(f"Sumstats filtrado: {args.out_sumstats}")
    print(f"VCF filtrado (texto): {args.out_vcf}")
    print(f"VCF comprimido: {vcf_gz}")
    print(f"Índice criado: {vcf_gz}.tbi")
    print(f"----------------------------------\n")

    print(f"{len(ids_comuns)} variantes em comum encontradas.")


if __name__ == "__main__":
    main()