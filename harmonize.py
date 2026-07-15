#!/usr/bin/env python3

import os
import argparse
import subprocess
import pandas as pd
import pysam


COMP = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C"
}


def complement(allele):
    return COMP.get(allele, allele)


def swap_gt(gt):
    """
    Inverte os genótipos:
    0 -> 1
    1 -> 0
    """

    if gt is None:
        return gt

    new_gt = []

    for allele in gt:

        if allele is None:
            new_gt.append(None)

        elif allele == 0:
            new_gt.append(1)

        elif allele == 1:
            new_gt.append(0)

        else:
            new_gt.append(allele)

    return tuple(new_gt)


def load_sumstats(sumstats_file, snp_col, ea_col, oa_col):

    df = pd.read_csv(sumstats_file, sep=None, engine="python")

    return {
        row[snp_col]: (
            str(row[ea_col]).upper(),
            str(row[oa_col]).upper()
        )
        for _, row in df.iterrows()
    }


def compress_and_index(vcf_file):

    gz_file = f"{vcf_file}.gz"
    # remove existing gz if present to avoid prompts/duplicates
    if os.path.exists(gz_file):
        os.remove(gz_file)

    subprocess.run([
        "bcftools",
        "view",
        "-Oz",
        "-o",
        gz_file,
        vcf_file
    ], check=True)

    subprocess.run([
        "bcftools",
        "index",
        "-t",
        gz_file
    ], check=True)

    os.remove(vcf_file)

    return gz_file


def harmonize_vcf(
    vcf_file,
    sumstats_dict,
    output_vcf
):
    vcf_in = pysam.VariantFile(vcf_file)

    vcf_out = pysam.VariantFile(
        output_vcf,
        "w",
        header=vcf_in.header
    )

    stats = {
        "match": 0,
        "swapped": 0,
        "flipped": 0,
        "flipped_swapped": 0,
        "removed": 0
    }

    changes = []

    for record in vcf_in:

        rsid = record.id

        if rsid not in sumstats_dict:
            continue

        ea, oa = sumstats_dict[rsid]

        orig_ref = record.ref.upper()
        orig_alt = record.alts[0].upper()
        chrom = getattr(record, 'chrom', None)
        pos = getattr(record, 'pos', None)

        ref = orig_ref
        alt = orig_alt

        #################################################
        # MATCH
        #################################################

        if ref == oa and alt == ea:

            stats["match"] += 1
            vcf_out.write(record)
            continue

        #################################################
        # SWAPPED
        #################################################

        if ref == ea and alt == oa:

            record.ref = oa
            record.alts = (ea,)

            for sample in record.samples:

                gt = record.samples[sample]["GT"]

                if gt is not None:
                    record.samples[sample]["GT"] = swap_gt(gt)

            stats["swapped"] += 1
            changes.append({
                "chrom": chrom,
                "pos": pos,
                "rsid": rsid,
                "ref_before": orig_ref,
                "alt_before": orig_alt,
                "ref_after": record.ref,
                "alt_after": record.alts[0],
                "change": "swapped"
            })
            vcf_out.write(record)
            continue

        #################################################
        # FLIPPED
        #################################################

        if (
            complement(ref) == oa and
            complement(alt) == ea
        ):

            record.ref = oa
            record.alts = (ea,)

            stats["flipped"] += 1
            changes.append({
                "chrom": chrom,
                "pos": pos,
                "rsid": rsid,
                "ref_before": orig_ref,
                "alt_before": orig_alt,
                "ref_after": record.ref,
                "alt_after": record.alts[0],
                "change": "flipped"
            })
            vcf_out.write(record)
            continue

        #################################################
        # FLIPPED + SWAPPED
        #################################################

        if (
            complement(ref) == ea and
            complement(alt) == oa
        ):

            record.ref = oa
            record.alts = (ea,)

            for sample in record.samples:

                gt = record.samples[sample]["GT"]

                if gt is not None:
                    record.samples[sample]["GT"] = swap_gt(gt)

            stats["flipped_swapped"] += 1
            changes.append({
                "chrom": chrom,
                "pos": pos,
                "rsid": rsid,
                "ref_before": orig_ref,
                "alt_before": orig_alt,
                "ref_after": record.ref,
                "alt_after": record.alts[0],
                "change": "flipped_swapped"
            })
            vcf_out.write(record)
            continue

        #################################################
        # REMOVE
        #################################################

        stats["removed"] += 1
        changes.append({
            "chrom": chrom,
            "pos": pos,
            "rsid": rsid,
            "ref_before": orig_ref,
            "alt_before": orig_alt,
            "ref_after": "",
            "alt_after": "",
            "change": "removed"
        })

    vcf_out.close()
    vcf_in.close()

    return stats, changes


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--vcf", required=True)

    parser.add_argument("--sumstats", required=True)

    parser.add_argument("--snp-col", required=True)

    parser.add_argument("--effect-col", required=True)

    parser.add_argument("--other-col", required=True)

    parser.add_argument("--out", required=True)

    args = parser.parse_args()

    print("Loading summary statistics\n")

    sumstats_dict = load_sumstats(
        args.sumstats,
        args.snp_col,
        args.effect_col,
        args.other_col
    )

    print("Harmonizing VCF\n")

    # Determine temporary uncompressed output path if user requested a .vcf.gz
    out_arg = args.out
    if out_arg.endswith('.vcf.gz'):
        tmp_out = out_arg[:-3]
    else:
        tmp_out = out_arg

    stats, changes = harmonize_vcf(
        args.vcf,
        sumstats_dict,
        tmp_out
    )

    gz_file = None
    # If user requested gzipped output, compress and index the temporary file
    if out_arg.endswith('.vcf.gz'):
        print("Compressing and indexing\n")
        # remove existing final gz if present to ensure overwrite
        if os.path.exists(out_arg):
            os.remove(out_arg)
        gz_file = compress_and_index(tmp_out)
    else:
        gz_file = tmp_out

    print(f"\n===== SUMMARY =====")
    for k, v in stats.items():
        print(f"{k}: {v:,}")

    # write changes file (only positions where something changed)
    change_file = None
    if out_arg.endswith('.vcf.gz'):
        change_file = out_arg.replace('.vcf.gz', '_changes.txt')
    else:
        change_file = f"{out_arg}_changes.txt"

    if changes:
        with open(change_file, 'w') as cf:
            cf.write('\t'.join(['chrom','pos','rsid','ref_before','alt_before','ref_after','alt_after','change']) + '\n')
            for c in changes:
                cf.write('\t'.join([str(c.get(col,'')) for col in ['chrom','pos','rsid','ref_before','alt_before','ref_after','alt_after','change']]) + '\n')

    print(f"\nOutput: {gz_file}")
    if changes:
        print(f"Changes file: {change_file}\n")


if __name__ == "__main__":
    main()