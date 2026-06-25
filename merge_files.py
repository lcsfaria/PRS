import argparse
import subprocess
import pandas as pd
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge harmonized VCFs and filtered sumstats by ancestry."
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Base output directory (same as OUTPUT_DIR in config.ini).",
    )
    parser.add_argument(
        "--sumstats",
        required=True,
        help=(
            "Comma-separated list of sumstats entries in the format "
            "'<sumstats_path>:<ancestry_name>' "
            "(same as SUMSTATS_FILE_ANCESTRY in config.ini)."
        ),
    )
    return parser.parse_args()


def merge_sumstats(sumstats_files, output_path):
    dfs = []
    for f in sumstats_files:
        print(f"  -> {f.name}")
        dfs.append(pd.read_csv(f, sep="\t"))

    merged = pd.concat(dfs, ignore_index=True)
    merged.to_csv(output_path, sep="\t", index=False)
    print(f"Sumstats mergeados salvos em: {output_path}")


def merge_vcfs(vcf_files, vcf_list_file, output_vcf):
    with open(vcf_list_file, "w") as f:
        for vcf in vcf_files:
            f.write(str(vcf) + "\n")

    print(f"\nConcatenando {len(vcf_files)} VCF(s) -> {output_vcf}")
    subprocess.run(
        ["bcftools", "concat", "-f", vcf_list_file, "-Oz", "-o", output_vcf],
        check=True,
    )

    print("Indexando VCF mergeado...")
    subprocess.run(["bcftools", "index", "-t", output_vcf], check=True)
    print(f"VCF mergeado e indexado: {output_vcf}")


def main():
    args = parse_args()

    output_dir = args.output_dir
    sumstats_entries = [s.strip() for s in args.sumstats.split(",")]

    for entry in sumstats_entries:
        parts = entry.split(":")
        ancestry_name = parts[-1]
        ancestry_dir = Path(output_dir) / ancestry_name

        print(f"\n=== Ancestralidade: {ancestry_name} ===")
        print(f"Diretório: {ancestry_dir}")

        vcf_files = sorted(ancestry_dir.glob("*harmonized*.vcf.gz"))
        sumstats_files = sorted(ancestry_dir.glob("*filtered_sumstats.txt"))

        print(f"VCFs encontrados: {[f.name for f in vcf_files]}")
        print(f"Sumstats encontrados: {[f.name for f in sumstats_files]}")

        if not vcf_files:
            print(f"[AVISO] Nenhum VCF harmonizado encontrado para {ancestry_name}. Pulando.")
            continue

        if not sumstats_files:
            print(f"[AVISO] Nenhum sumstats filtrado encontrado para {ancestry_name}. Pulando.")
            continue

        # Merge sumstats
        merged_sumstats_path = ancestry_dir / f"{ancestry_name}_merged_sumstats.txt"
        merge_sumstats(sumstats_files, merged_sumstats_path)

        # Merge VCFs
        vcf_list_file = ancestry_dir / f"{ancestry_name}_vcf_list.txt"
        output_vcf = ancestry_dir / f"{ancestry_name}_All_Chr_merged.vcf.gz"
        merge_vcfs(vcf_files, vcf_list_file, str(output_vcf))


if __name__ == "__main__":
    main()