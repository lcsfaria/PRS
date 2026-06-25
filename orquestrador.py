import os
from pathlib import Path
import subprocess
from configparser import ConfigParser
import gzip
import shutil
import pandas as pd
import argparse

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
    if not vcf_file.endswith('.gz'):
        return vcf_file, False

    os.makedirs(temp_dir, exist_ok=True)
    uncompressed_vcf = os.path.join(temp_dir, os.path.basename(vcf_file).replace('.vcf.gz', '.vcf'))

    if os.path.exists(uncompressed_vcf):
        return uncompressed_vcf, False

    print(f"Descompactando VCF para uso do vcf_maker: {vcf_file} -> {uncompressed_vcf}")
    with gzip.open(vcf_file, 'rt') as fin, open(uncompressed_vcf, 'wt') as fout:
        shutil.copyfileobj(fin, fout)

    return uncompressed_vcf, True

def resolve_vcf_file(vcf_pattern, temp_dir):
    if not vcf_pattern.endswith('.gz') and os.path.exists(vcf_pattern):
        return vcf_pattern, False
    if os.path.exists(vcf_pattern + '.gz'):
        uncompressed, was_decompressed = ensure_uncompressed_vcf(vcf_pattern + '.gz', temp_dir)
        return uncompressed, was_decompressed
    if os.path.exists(vcf_pattern):
        return vcf_pattern, False
    raise FileNotFoundError(f"Não encontrado: {vcf_pattern} ou {vcf_pattern}.gz")

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


# ─── Argparse ────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Pipeline de ancestralidade por cromossomo.")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("--chroms",nargs="+",type=int,metavar="N",help="Cromossomos específicos. Ex: --chroms 1 3 21 22"
)
group.add_argument("--chrom-range",nargs=2,type=int,metavar=("START", "END"),help="Intervalo de cromossomos. Ex: --chrom-range 1 22"
)
args = parser.parse_args()

if args.chroms:
    chromosomes = args.chroms
else:
    start, end = args.chrom_range
    if start > end:
        parser.error(f"START ({start}) não pode ser maior que END ({end}).")
    chromosomes = list(range(start, end + 1))

print(f"Cromossomos analisados: {chromosomes}\n")
# ─────────────────────────────────────────────────────────────────────────────


script_dir = os.getcwd()
config_path = os.path.join(script_dir, "config.ini")
config = ConfigParser()
config.read(config_path)

project_name = config.get("Project", "PROJECT_NAME")

ancestry_splitter = config.get("Scripts", "ANCESTRY_SPLITTER")
filtrar_variantes = config.get("Scripts", "FILTRAR_VARIANTES")

base_vcf_file = config.get("Paths", "VCF_FILE")
msp_file = config.get("Paths", "MSP_FILE")
output_dir = config.get("Paths", "OUTPUT_DIR")
sumstats = config.get("Paths", "SUMSTATS_FILE_ANCESTRY").split(",")

temp_dir = os.path.join(output_dir, "temp")

os.makedirs(output_dir, exist_ok=True)
os.makedirs(temp_dir, exist_ok=True)


for chrom in chromosomes:
    base_vcf_input = config.get("Paths", "VCF_FILE")
    base_vcf_input = base_vcf_input.replace("@", str(chrom))

    vcf_to_use, was_decompressed = resolve_vcf_file(base_vcf_input, temp_dir)

    msp_file_chrom = config.get("Paths", "MSP_FILE")
    msp_file_chrom = msp_file_chrom.replace("@", str(chrom))

    print(f"\n=== Cromossomo {chrom} ===")
    print(f"VCF de entrada: {vcf_to_use}\n")
    print('Executando vcf_maker:\n')
    vcf_maker(ancestry_splitter, vcf_to_use, msp_file_chrom, os.path.join(output_dir, temp_dir), project_name, chrom)

    if os.path.abspath(vcf_to_use).startswith(os.path.abspath(temp_dir)):
        try:
            os.remove(vcf_to_use)
            print(f"VCF temporário removido: {vcf_to_use}")
        except OSError as e:
            print(f"Aviso: não foi possível remover {vcf_to_use}: {e}")
    '''        
    if was_decompressed:
        try:
            os.remove(vcf_to_use)
        except OSError:
            pass
    '''
    for info_sumstat in sumstats:
        ancestry_name = info_sumstat.split(":")
        os.makedirs(os.path.join(output_dir, ancestry_name[-1]), exist_ok=True)

        print(f'Movendo arquivos VCF para pasta {ancestry_name[-1]}\n')
        mover = f'mv {temp_dir}/{project_name}_Chr{chrom}_{ancestry_name[-1]}.vcf {output_dir}/{ancestry_name[-1]}'
        executa(mover)

        print(f'Comprimindo arquivo VCF {output_dir}/{ancestry_name[-1]}/{project_name}_Chr{chrom}_{ancestry_name[-1]}.vcf\n')
        zip_1 = f'bgzip -f {output_dir}/{ancestry_name[-1]}/{project_name}_Chr{chrom}_{ancestry_name[-1]}.vcf'
        executa(zip_1)

        intersection_filter(
            filtrar_variantes,
            f'{output_dir}/{ancestry_name[-1]}/{project_name}_Chr{chrom}_{ancestry_name[-1]}.vcf.gz',
            info_sumstat.split(":")[0],'rsid',f'{output_dir}/{ancestry_name[-1]}/{project_name}_Chr{chrom}_{ancestry_name[-1]}_filtered.vcf',
            f'{output_dir}/{ancestry_name[-1]}/{project_name}_Chr{chrom}_{ancestry_name[-1]}_filtered_sumstats.txt'
        )

        print(f'Removendo arquivo VCF {output_dir}/{ancestry_name[-1]}/{project_name}_Chr{chrom}_{ancestry_name[-1]}.vcf.gz\n')
        remover_filtered = f'rm {output_dir}/{ancestry_name[-1]}/{project_name}_Chr{chrom}_{ancestry_name[-1]}.vcf.gz'
        executa(remover_filtered)

        print('Harmonizando VCF final com harmonize.py...')
        harmonize_cmd = (
            f'python {os.path.join(script_dir, "harmonize.py")} '
            f'--vcf {output_dir}/{ancestry_name[-1]}/{project_name}_Chr{chrom}_{ancestry_name[-1]}_filtered.vcf.gz '
            f'--sumstats {output_dir}/{ancestry_name[-1]}/{project_name}_Chr{chrom}_{ancestry_name[-1]}_filtered_sumstats.txt '
            f'--snp-col rsid --effect-col effect_allele --other-col other_allele '
            f'--out {output_dir}/{ancestry_name[-1]}/{project_name}_Chr{chrom}_{ancestry_name[-1]}_filtered_harmonized.vcf.gz'
        )
        executa(harmonize_cmd)

        remover_filtered = f'rm {output_dir}/{ancestry_name[-1]}/{project_name}_Chr{chrom}_{ancestry_name[-1]}_filtered.vcf'
        executa(remover_filtered)


#aqui junta arquivos de sumstats e VCFs por ancestralidade
for info_sumstat in sumstats:
    ancestry_name = info_sumstat.split(":")

    vcf_files = sorted(Path(f'{output_dir}/{ancestry_name[-1]}/').glob("*harmonized*.vcf.gz"))
    sumstats_files = sorted(Path(f'{output_dir}/{ancestry_name[-1]}/').glob("*filtered_sumstats.txt"))
    print(vcf_files)
    print(sumstats_files)

    dfs = []
    for f in sumstats_files:
        print(f"  -> {f.name}")
        dfs.append(pd.read_csv(f, sep="\t"))

    merged_sumstats = pd.concat(dfs, ignore_index=True)
    merged_sumstats.to_csv(
        f'{output_dir}/{ancestry_name[-1]}/{ancestry_name[-1]}_merged_sumstats.txt',
        sep="\t", index=False
    )

    vcf_list_file = f"{output_dir}/{ancestry_name[-1]}/{ancestry_name[-1]}_vcf_list.txt"
    OUTPUT_VCF = f"{output_dir}/{ancestry_name[-1]}/{ancestry_name[-1]}_All_Chr_merged.vcf.gz"

    with open(vcf_list_file, "w") as f:
        for vcf in vcf_files:
            f.write(str(vcf) + "\n")

    print("\nConcatenando VCFs...")
    subprocess.run(["bcftools", "concat", "-f", vcf_list_file, "-Oz", "-o", OUTPUT_VCF], check=True)
    subprocess.run(["bcftools", "index", "-t", OUTPUT_VCF], check=True)
    subprocess.run(['plink', '--vcf', OUTPUT_VCF, '--make-bed', '--vcf-half-call', 'r',
        '--double-id', '--out', f"{output_dir}/{ancestry_name[-1]}/{ancestry_name[-1]}_All_Chr_merged"])
