#!/usr/bin/env python3

import pandas as pd
import argparse
import os


'''
exemplo de uso

python qc_sumstats.py \
    --input Rizig_et_al_2023_AFR_AAC_metaGWAS_no23andMe_hg38.txt \
    --output Rizig_QCed.txt \
    --snp-col rsid \
    --effect-col effect_allele \
    --other-col other_allele \
    --chr-col chr \
    --bp-col pos \
    --beta-col beta \
    --or-col OR \
    --se-col SE \
    --p-col P \
    --maf-col effect_allele_frequency \
    --info-col INFO \
    --n-col N

O usuário deve fornecer explicitamente os nomes das colunas no arquivo de entrada.
''' 


def normalize_dataframe_headers(df):
    df.columns = [str(col).strip() for col in df.columns]
    return df


def normalize_id_values(series):
    return series.astype(str).str.strip()


def qc_sumstats(
    infile,
    outfile,
    snp_col,
    effect_col,
    other_col,
    chr_col,
    bp_col,
    beta_col,
    or_col=None,
    maf_col=None,
    info_col=None,
    p_col=None,
    se_col=None,
    n_col=None,
    maf_threshold=0.01,
    info_threshold=0.8
):

    print("Lendo arquivo...")
    df = pd.read_csv(infile, sep=None, engine="python")
    df = normalize_dataframe_headers(df)

    required = [snp_col, effect_col, other_col, chr_col, bp_col, beta_col]
    if or_col is not None:
        required.append(or_col)
    if maf_col is not None:
        required.append(maf_col)
    if info_col is not None:
        required.append(info_col)
    if p_col is not None:
        required.append(p_col)
    if se_col is not None:
        required.append(se_col)
    if n_col is not None:
        required.append(n_col)

    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(
            f"Colunas não encontradas no sumstats: {missing}. Colunas disponíveis: {list(df.columns)}"
        )

    df[snp_col] = normalize_id_values(df[snp_col])
    df[effect_col] = normalize_id_values(df[effect_col]).str.upper()
    df[other_col] = normalize_id_values(df[other_col]).str.upper()

    n0 = len(df)

    df = df.drop_duplicates(subset=snp_col, keep=False)
    n1 = len(df)

    ambiguous = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")}
    df = df[~df[[effect_col, other_col]].apply(tuple, axis=1).isin(ambiguous)]
    n2 = len(df)

    if maf_col is not None:
        df = df[df[maf_col] > maf_threshold]
    n3 = len(df)

    if info_col is not None:
        df = df[df[info_col] > info_threshold]
    n4 = len(df)

    rename_map = {
        snp_col: "rsid",
        chr_col: "chr",
        bp_col: "pos",
        effect_col: "effect_allele",
        other_col: "other_allele",
        beta_col: "beta",
    }
    if or_col is not None:
        rename_map[or_col] = "or"
    if maf_col is not None:
        rename_map[maf_col] = "maf"
    if info_col is not None:
        rename_map[info_col] = "info"
    if p_col is not None:
        rename_map[p_col] = "p"
    if se_col is not None:
        rename_map[se_col] = "se"
    if n_col is not None:
        rename_map[n_col] = "n"

    df = df.rename(columns=rename_map)

    preferred_order = [
        "chr", "pos", "rsid", "effect_allele", "other_allele",
        "beta", "or", "se", "p", "maf", "info", "n"
    ]
    output_columns = [c for c in preferred_order if c in df.columns] + [c for c in df.columns if c not in preferred_order]

    output_dir = os.path.dirname(outfile)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    df.to_csv(outfile, sep="\t", index=False, columns=output_columns)

    print(f"Inicial: {n0:,}")
    print(f"Após remover duplicados: {n1:,}")
    print(f"Após remover AT/TA/CG/GC: {n2:,}")
    print(f"Após filtro de frequência: {n3:,}")
    print(f"Após filtro INFO: {n4:,}")
    print(f"Arquivo salvo em: {outfile}")


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)

    parser.add_argument("--snp-col", required=True, help="Nome da coluna de SNP (rsID)")
    parser.add_argument("--effect-col", required=True, help="Nome da coluna de alelo de efeito")
    parser.add_argument("--other-col", required=True, help="Nome da coluna de outro alelo")
    parser.add_argument("--chr-col", required=True, help="Nome da coluna de cromossomo")
    parser.add_argument("--bp-col", required=True, help="Nome da coluna de posição basepair")
    parser.add_argument("--beta-col", required=True, help="Nome da coluna de beta")
    parser.add_argument("--or-col", required=False, help="Nome da coluna de odds ratio (opcional)")
    parser.add_argument("--se-col", required=False, help="Nome da coluna de erro padrão (opcional)")
    parser.add_argument("--p-col", required=False, help="Nome da coluna de p-valor (opcional)")
    parser.add_argument("--maf-col", required=False, help="Nome da coluna de MAF/frequência (opcional)")
    parser.add_argument("--info-col", required=False, help="Nome da coluna de INFO de imputação (opcional)")
    parser.add_argument("--n-col", required=False, help="Nome da coluna de tamanho de amostra (opcional)")

    parser.add_argument("--maf", type=float, default=0.01, help="Filtro de MAF - Default: 0.01 (Opcional)")
    parser.add_argument("--info", type=float, default=0.8, help="Filtro de INFO - Default: 0.8 (Opcional)")

    args = parser.parse_args()

    qc_sumstats(
        args.input,
        args.output,
        args.snp_col,
        args.effect_col,
        args.other_col,
        args.chr_col,
        args.bp_col,
        args.beta_col,
        args.or_col,
        args.maf_col,
        args.info_col,
        args.p_col,
        args.se_col,
        args.n_col,
        args.maf,
        args.info
    )


if __name__ == "__main__":
    main()