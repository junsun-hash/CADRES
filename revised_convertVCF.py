#!/usr/bin/env python
import sys
import argparse

def convert_vcf(vcf_path, output_path):
    """
    Parses a VCF file to extract allele depths and reformat specific fields.
    This script replaces the revised_convertVCF.sh shell script.

    Args:
        vcf_path (str): Path to the input VCF file.
        output_path (str): Path to the output file.
    """
    try:
        with open(vcf_path, 'r') as vcf_file, open(output_path, 'w') as out_file:
            for line in vcf_file:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 10:
                    continue

                # Fields from VCF: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE
                chrom, pos, _, ref, alt, *_ , sample_data = fields

                # Expected format field: GT:AD:DP:GQ:PL
                # We need the Allele Depth (AD) which is the 2nd field
                format_keys = fields[8].split(':')
                sample_values = sample_data.split(':')

                try:
                    ad_index = format_keys.index('AD')
                    allele_depths = sample_values[ad_index].split(',')
                    ref_depth = allele_depths[0]
                    alt_depth = allele_depths[1]
                    total_depth = int(ref_depth) + int(alt_depth)
                    
                    # Create the required output fields
                    depths_field = f"{total_depth},{alt_depth}"
                    alt_fraction = "0.01" # This was hardcoded in the original awk script

                    # Write to output file in the specified format
                    out_file.write(f"{chrom}\t{pos}\t{depths_field}\t{ref}\t{alt}\t{alt_fraction}\n")

                except (ValueError, IndexError) as e:
                    # Skip lines where AD is not found or is malformed
                    sys.stderr.write(f"Skipping line due to format error: {line.strip()} | Error: {e}\n")
                    continue

    except IOError as e:
        sys.stderr.write(f"Error processing files: {e}\n")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Convert VCF file format. Python replacement for revised_convertVCF.sh.")
    parser.add_argument("input_vcf", help="Input VCF file path.")
    parser.add_argument("output_file", help="Output file path.")
    # The quality_filter argument from the original .sh script was not used in the awk part, so it's omitted.
    
    args = parser.parse_args()
    
    convert_vcf(args.input_vcf, args.output_file)
    print(f"Conversion complete. Output written to {args.output_file}")

if __name__ == "__main__":
    main()