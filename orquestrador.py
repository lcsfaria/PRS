import os
from pathlib import Path
import subprocess
from configparser import ConfigParser
import gzip
import shutil
import pandas as pd
import argparse
from concurrent.futures import ProcessPoolExecutor

'''
como rodar 
python orquestrador.py --chroms 1 3 21 22
ou
python orquestrador.py --chrom-range 1 22
'''

def executa(comando):
    print(f"{comando}\n")
    subprocess.run(comando, shell=True, check=True)

def ensure_uncompressed_vcf(vcf_file, temp_dir):
    if not vcf_file.endswith(".gz"):
        return vcf_file, False

    os.makedirs(temp_dir, exist_ok=True)

    uncompressed_vcf = os.path.join(temp_dir, os.path.basename(vcf_file).replace(".vcf.gz", ".vcf"))

    if os.path.exists(uncompressed_vcf):
        return uncompressed_vcf, False

    print(f"Decompressing VCF:\n" f"  {vcf_file}\n"  f"  -> {uncompressed_vcf}")

    with open(uncompressed_vcf, "wb") as fout:
        subprocess.run(["gunzip", "-c", vcf_file],stdout=fout,check=True)

    return uncompressed_vcf, True

def resolve_vcf_file(vcf_pattern, temp_dir):

    if vcf_pattern.endswith('.gz'):
        print('Compressed VCF file provided. Attempting to decompress it\n')
        if os.path.exists(vcf_pattern):
            return ensure_uncompressed_vcf(vcf_pattern, temp_dir)
        raise FileNotFoundError(f"Não encontrado: {vcf_pattern}")

    # usuário passou .vcf
    if os.path.exists(vcf_pattern):
        print('Uncompressed VCF file provided\n')
        return vcf_pattern, False

    # existe apenas a versão comprimida
    if os.path.exists(vcf_pattern + '.gz'):
        return ensure_uncompressed_vcf(vcf_pattern + '.gz', temp_dir)

    raise FileNotFoundError(
        f"File not found: {vcf_pattern} ou {vcf_pattern}.gz"
    )

def vcf_maker(ancestry_splitter, vcf_file, msp_file, output_folder, project_name, chrom):
    comando = f'python {ancestry_splitter} --vcf {vcf_file} --msp {msp_file} --output_prefix {output_folder}/{project_name}_Chr{chrom}'
    executa(comando)

def intersection_filter(filtrar_variantes, input_vcf, sumstats, snp_col, out_vcf, out_sumstats):
    comando = f'python {filtrar_variantes} --vcf {input_vcf} --sumstats {sumstats} --snp-col {snp_col} --out-vcf {out_vcf} --out-sumstats {out_sumstats}'
    executa(comando)

def compress_and_index_vcf(vcf_file):
    gz_file = f"{vcf_file}.gz"
    subprocess.run(["bcftools", "view", "-Oz", "-o", gz_file, vcf_file], check=True)
    subprocess.run(["bcftools", "index", "-f", gz_file], check=True)


def process_chromosome(chrom, script_dir, config_path):
    config = ConfigParser()
    config.read(config_path)

    project_name = config.get("Project", "PROJECT_NAME")
    ancestry_splitter = os.path.join(script_dir, "vcf_maker.py")
    filtrar_variantes = os.path.join(script_dir, "filtrar_variantes.py")

    base_vcf_file = config.get("Paths", "VCF_FILE")
    msp_file = config.get("Paths", "MSP_FILE")
    output_dir = config.get("Paths", "OUTPUT_DIR")
    sumstats = config.get("Paths", "SUMSTATS_FILE_ANCESTRY").split(",")

    chrom_temp_dir = os.path.join(output_dir, "temp", f"chrom_{chrom}")
    os.makedirs(chrom_temp_dir, exist_ok=True)

    base_vcf_input = base_vcf_file.replace("@", str(chrom))
    vcf_to_use, _ = resolve_vcf_file(base_vcf_input, chrom_temp_dir)

    msp_file_chrom = msp_file.replace("@", str(chrom))

    print(f"\n=== Chromosome {chrom} ===")
    print(f"Input VCF: {vcf_to_use}\n")
    print('Running vcf_maker:\n')
    vcf_maker(ancestry_splitter, vcf_to_use, msp_file_chrom, chrom_temp_dir, project_name, chrom)

    if os.path.abspath(vcf_to_use).startswith(os.path.abspath(chrom_temp_dir)):
        try:
            os.remove(vcf_to_use)
            print(f"VCF temporário removido: {vcf_to_use}")
        except OSError as e:
            print(f"It was not possible to remove {vcf_to_use}: {e}\n")

    for info_sumstat in sumstats:
        ancestry_name = info_sumstat.split(":")
        ancestry_out = os.path.join(output_dir, ancestry_name[-1])
        os.makedirs(ancestry_out, exist_ok=True)

        print(f'Moving VCF to {ancestry_name[-1]}\n')
        mover = f'mv {chrom_temp_dir}/{project_name}_Chr{chrom}_{ancestry_name[-1]}.vcf {ancestry_out}'
        executa(mover)

        print(f'Compressing VCF {ancestry_out}/{project_name}_Chr{chrom}_{ancestry_name[-1]}.vcf\n')
        zip_1 = f'bgzip -f {ancestry_out}/{project_name}_Chr{chrom}_{ancestry_name[-1]}.vcf'
        executa(zip_1)

        print(f'Creating the intersection between VCF and Sumstats for {ancestry_name[-1]}...\n')
        intersection_filter(
            filtrar_variantes,
            f'{ancestry_out}/{project_name}_Chr{chrom}_{ancestry_name[-1]}.vcf.gz',
            info_sumstat.split(":")[0],
            'rsid',
            f'{ancestry_out}/{project_name}_Chr{chrom}_{ancestry_name[-1]}_filtered.vcf',
            f'{ancestry_out}/{project_name}_Chr{chrom}_{ancestry_name[-1]}_filtered_sumstats.txt'
        )
        print('Intersection Created\n')

        print(f'Removing VCF {ancestry_out}/{project_name}_Chr{chrom}_{ancestry_name[-1]}.vcf.gz to save space\n')
        remover_filtered = f'rm {ancestry_out}/{project_name}_Chr{chrom}_{ancestry_name[-1]}.vcf.gz'
        executa(remover_filtered)

        print('Harmonizing final  VCF\n')
        harmonize_cmd = (
            f'python {os.path.join(script_dir, "harmonize.py")} '
            f'--vcf {ancestry_out}/{project_name}_Chr{chrom}_{ancestry_name[-1]}_filtered.vcf.gz '
            f'--sumstats {ancestry_out}/{project_name}_Chr{chrom}_{ancestry_name[-1]}_filtered_sumstats.txt '
            f'--snp-col rsid --effect-col effect_allele --other-col other_allele '
            f'--out {ancestry_out}/{project_name}_Chr{chrom}_{ancestry_name[-1]}_filtered_harmonized.vcf.gz'
        )
        executa(harmonize_cmd)

        print(f'Removing VCF {ancestry_out}/{project_name}_Chr{chrom}_{ancestry_name[-1]}_filtered.vcf.gz to save space\n')
        remover_filtered = f'rm {ancestry_out}/{project_name}_Chr{chrom}_{ancestry_name[-1]}_filtered.vcf*'
        executa(remover_filtered)


def main():
    parser = argparse.ArgumentParser(description="Pipeline de ancestralidade por cromossomo.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--chroms", nargs="+", type=int, metavar="N", help="Cromossomos específicos. Ex: --chroms 1 3 21 22")
    group.add_argument("--chrom-range", nargs=2, type=int, metavar=("START", "END"), help="Intervalo de cromossomos. Ex: --chrom-range 1 22")
    parser.add_argument("--workers", type=int, default=1, help="Número de cromossomos a processar em paralelo (padrão: 1)")
    args = parser.parse_args()

    if args.chroms:
        chromosomes = args.chroms
    else:
        start, end = args.chrom_range
        if start > end:
            parser.error(f"START ({start}) não pode ser maior que END ({end}).")
        chromosomes = list(range(start, end + 1))

    print(f"Cromossomos analisados: {chromosomes}")
    print(f"Workers: {args.workers}\n")

    script_dir = os.getcwd()
    config_path = os.path.join(script_dir, "config.ini")

    if args.workers > 1 and len(chromosomes) > 1:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            list(executor.map(process_chromosome, chromosomes, [script_dir] * len(chromosomes), [config_path] * len(chromosomes)))
    else:
        for chrom in chromosomes:
            process_chromosome(chrom, script_dir, config_path)


if __name__ == "__main__":
    main()



