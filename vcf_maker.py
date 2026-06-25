import argparse
import sys
from collections import defaultdict
import tempfile
import os

def parse_msp_header(msp_header_line):
    ancestry_map = {}
    parts = msp_header_line.strip().split(': ')[1].split()
    for part in parts:
        name, code = part.split('=')
        ancestry_map[int(code)] = name
    return ancestry_map

def process_ancestry_files(vcf_path, msp_path, output_prefix):
    output_files = {}
    processed_snps = 0
    try:
        with open(vcf_path, 'r') as vcf_file, open(msp_path, 'r') as msp_file:
            vcf_header_lines = []
            individuals = []
            for line in vcf_file:
                if line.startswith('#'):
                    vcf_header_lines.append(line)
                    if line.startswith("#CHROM"):
                        individuals = line.strip().split('\t')[9:]
                else:
                    # The first non-header line is where data starts
                    break
            
            if not individuals:
                print("ERROR: Could not find VCF header line starting with '#CHROM'.", file=sys.stderr)
                return False

            # Find and parse the MSP header
            ancestry_map = None
            for line in msp_file:
                if line.startswith("#Subpopulation order/codes"):
                    ancestry_map = parse_msp_header(line)
                    break
            if ancestry_map is None:
                print("ERROR: Could not find '#Subpopulation order/codes' header in MSP file.", file=sys.stderr)
                return False

            # Prepare Output VCF Files
            for anc_name in ancestry_map.values():
                filename = f"{output_prefix}_{anc_name}.vcf"
                file_handle = open(filename, 'w')

                for header_line in vcf_header_lines:
                    file_handle.write(header_line)
                output_files[anc_name] = file_handle

            # Synchronized Line-by-Line Processing
            msp_iterator = (line for line in msp_file if not line.startswith('#'))
            
            vcf_file.seek(0)
            vcf_iterator = (line for line in vcf_file if not line.startswith('#'))

            current_msp_fields = next(msp_iterator).strip().split('\t')
            
            for vcf_line in vcf_iterator:
                vcf_fields = vcf_line.strip().split('\t')
                
                # These are the standard VCF columns before the genotype data
                fixed_vcf_cols = vcf_fields[:9]
                snp_pos = int(vcf_fields[1])
                genotypes = [gt.split(':')[0] for gt in vcf_fields[9:]]

                # Advance the MSP reader until its region covers the current SNP
                while snp_pos > int(current_msp_fields[2]):
                    try:
                        current_msp_fields = next(msp_iterator).strip().split('\t')
                    except StopIteration:
                        break
                
                if snp_pos >= int(current_msp_fields[1]):
                    msp_ancestries = current_msp_fields[6:]
                    ancestry_genotypes = defaultdict(list)

                    for i, ind_gt in enumerate(genotypes):
                        haplotype_alleles = ind_gt.split('|')
                        if len(haplotype_alleles) != 2: continue

                        anc_hap1_code = int(msp_ancestries[i*2])
                        anc_hap2_code = int(msp_ancestries[i*2 + 1])
                        
                        anc_hap1_name = ancestry_map.get(anc_hap1_code, None)
                        anc_hap2_name = ancestry_map.get(anc_hap2_code, None)

                        for target_anc_name in ancestry_map.values():
                            new_gt1 = haplotype_alleles[0] if anc_hap1_name == target_anc_name else '.'
                            new_gt2 = haplotype_alleles[1] if anc_hap2_name == target_anc_name else '.'
                            # Re-add the other format fields if they exist (e.g., :DS:HDS)
                            original_format_fields = vcf_fields[9+i].split(':')
                            new_gt_fields = [f"{new_gt1}|{new_gt2}"] + original_format_fields[1:]
                            ancestry_genotypes[target_anc_name].append(":".join(new_gt_fields))

                    # Write the constructed VCF lines to the respective output files
                    for anc_name, genotypes_list in ancestry_genotypes.items():
                        # Join the fixed columns with the modified genotype columns
                        output_line = "\t".join(fixed_vcf_cols + genotypes_list) + "\n"
                        output_files[anc_name].write(output_line)
                    
                    processed_snps += 1
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        return False
    finally:
        for f in output_files.values():
            f.close()
        print(f"INFO: Processing complete. Processed a total of {processed_snps} SNPs.")
    return True

def run_head_test(vcf_path, msp_path, num_lines=100):
    print(f"INFO: Running test on the first {num_lines} data lines of the input files.")
    print("INFO: Test files will be saved in the current directory.")

    test_prefix = "test_run"
    temp_vcf_path = f"{test_prefix}_input.vcf"
    temp_msp_path = f"{test_prefix}_input.msp"
    output_prefix = f"{test_prefix}_output"

    vcf_header_lines = []
    individuals = []
    ancestries = {}

    try:
        print(f"INFO: Creating subset file: {temp_vcf_path}")
        with open(vcf_path, 'r') as real_f, open(temp_vcf_path, 'w') as temp_f:
            line_count = 0
            for line in real_f:
                if line.startswith('#'):
                    vcf_header_lines.append(line)
                    if line.startswith("#CHROM"):
                        individuals = line.strip().split('\t')[9:]
                    temp_f.write(line)
                else:
                    if line_count >= num_lines:
                        break
                    temp_f.write(line)
                    line_count += 1
        
        print(f"INFO: Creating subset file: {temp_msp_path}")
        with open(msp_path, 'r') as real_f, open(temp_msp_path, 'w') as temp_f:
            line_count = 0
            for line in real_f:
                if line.startswith('#'):
                    if line.startswith("#Subpopulation order/codes"):
                        ancestries = parse_msp_header(line)
                    temp_f.write(line)
                else:
                    if line_count >= num_lines:
                        break
                    temp_f.write(line)
                    line_count += 1
    except FileNotFoundError as e:
        print(f"ERROR: Cannot find input file for test: {e}", file=sys.stderr)
        return

    if not individuals or not ancestries:
        print("ERROR: Could not parse required headers from input files. Check file format.", file=sys.stderr)
        return

    print("INFO: Executing main processing logic...")
    success = process_ancestry_files(temp_vcf_path, temp_msp_path, output_prefix)

    if not success:
        print("\nERROR: The script encountered an error during processing.")
        return

    # Validate the output VCF files
    print("INFO: Validating output files...")
    all_tests_passed = True
    for anc_name in ancestries.values():
        output_filename = f"{output_prefix}_{anc_name}.vcf"
        if not os.path.exists(output_filename):
            print(f"Test for {anc_name}: FAIL - Output file not created.")
            all_tests_passed = False
            continue

        with open(output_filename, 'r') as f:
            first_line = f.readline()
            if "##fileformat=VCF" in first_line:
                print(f"Test for {anc_name}: PASS - File is in VCF format.")
            else:
                print(f"Test for {anc_name}: FAIL - File does not have a valid VCF header.")
                all_tests_passed = False
    
    print("-" * 20)
    if all_tests_passed:
        print("Test finished successfully.")
        print("The following files have been saved in your directory:")
        print(f"  - {temp_vcf_path}")
        print(f"  - {temp_msp_path}")
        for anc_name in ancestries.values():
            print(f"  - {output_prefix}_{anc_name}.vcf")
    else:
        print("Some tests failed. Check the errors above.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Split a VCF file into ancestry-specific VCF files using an MSP file for local ancestry inference.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('--vcf', type=str, help="Path to the input VCF file.")
    parser.add_argument('--msp', type=str, help="Path to the input MSP file.")
    parser.add_argument('--output_prefix', type=str, help="Prefix for the output files (e.g., 'ancestry_results').")
    parser.add_argument('--test', action='store_true', help="Run a test on the first N lines of the provided files.")
    parser.add_argument('--test-lines', type=int, default=100, help="Number of data lines to use for the test (default: 100).")

    args = parser.parse_args()

    if args.test:
        if not args.vcf or not args.msp:
            print("ERROR: For --test, you must provide --vcf and --msp arguments.", file=sys.stderr)
            parser.print_help()
        else:
            run_head_test(args.vcf, args.msp, args.test_lines)
    elif args.vcf and args.msp and args.output_prefix:
        print(f"INFO: Starting full run with VCF='{args.vcf}', MSP='{args.msp}', PREFIX='{args.output_prefix}'")
        process_ancestry_files(args.vcf, args.msp, args.output_prefix)
    else:
        parser.print_help()