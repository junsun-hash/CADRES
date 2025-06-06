#!/usr/bin/env python
import argparse
import sys
from pyfaidx import Fasta

def check_homopolymer(sequence, edit_nuc):
    """
    Checks if the variant is in a homopolymer region based on the logic
    from the original Perl script.
    The sequence is 10bp long with the variant at the 5th position (index 4).
    The original logic checks 5 windows of 5 bases.
    """
    s = sequence.upper()
    edit_nuc = edit_nuc.upper()
    
    # The variant base is at index 4. The surrounding sequence is s[0:4] and s[5:9].
    # The original Perl script had a bug in its logic, using splitseq[5] which should be the base after the variant.
    # We will replicate the intended logic which is to check 5-base windows around the variant.
    # Window 1: -4, -3, -2, -1, variant
    if all(base == edit_nuc for base in (s[0], s[1], s[2], s[3])):
        return True
    # Window 2: -3, -2, -1, variant, +1
    if all(base == edit_nuc for base in (s[1], s[2], s[3], s[5])):
        return True
    # Window 3: -2, -1, variant, +1, +2
    if all(base == edit_nuc for base in (s[2], s[3], s[5], s[6])):
        return True
    # Window 4: -1, variant, +1, +2, +3
    if all(base == edit_nuc for base in (s[3], s[5], s[6], s[7])):
        return True
    # Window 5: variant, +1, +2, +3, +4
    if all(base == edit_nuc for base in (s[5], s[6], s[7], s[8])):
        return True
        
    return False


def filter_homopolymers(infile, outfile, refgenome_path):
    """
    Filters out variants located in homopolymer regions.
    Python replacement for filter_homopolymer_nucleotides.pl.
    """
    try:
        # Load reference genome. pyfaidx will automatically look for a .fai index
        # or create one if it doesn't exist.
        genome = Fasta(refgenome_path)
    except Exception as e:
        sys.stderr.write(f"Error: Could not open or index reference genome file {refgenome_path}. Please ensure it is a valid FASTA file.\nError: {e}\n")
        sys.exit(1)

    left_buffer = 4
    right_buffer = 4

    with open(infile, 'r') as sites_file, \
         open(outfile, 'w') as output_file, \
         open(f"{outfile}_failed", 'w') as failed_file:
        
        for line in sites_file:
            line = line.strip()
            if not line:
                continue
            
            fields = line.split()
            chrom, pos, edit_nuc = fields[0], int(fields[1]), fields[4]

            # 1-based position from file, convert to 0-based for pyfaidx
            start = pos - 1 - left_buffer
            end = pos + right_buffer
            
            try:
                # Fetch sequence. pyfaidx handles chromosome names like 'chr1' or '1'.
                sequence = genome[chrom][start:end].seq
            except KeyError:
                sys.stderr.write(f"Warning: Chromosome '{chrom}' not found in reference genome. Skipping variant at {chrom}:{pos}.\n")
                continue
            except Exception as e:
                sys.stderr.write(f"Warning: Could not fetch sequence for {chrom}:{start}-{end}. Skipping. Error: {e}\n")
                continue

            # The variant base itself is not in our fetched sequence, which is s[0:4] and s[5:9]. Let's verify.
            # pos is 1-based. The base *at* pos is what we're testing.
            # Our fetched sequence is from pos-4 to pos+4 (0-based: pos-1-4 to pos+4). The length is 9.
            # Example: pos=100. fetch 95-104. Sequence is for bases 96,97,98,99, 101,102,103,104.
            # This seems to match the perl script logic which reconstructs a 9-char string around the variant.
            
            is_homopolymer = check_homopolymer(sequence, edit_nuc)

            if not is_homopolymer:
                output_file.write(line + '\n')
            else:
                failed_file.write(line + '\n')

def main():
    parser = argparse.ArgumentParser(
        description="Homopolymer filter. Python replacement for filter_homopolymer_nucleotides.pl.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-i", "--infile", required=True, help="File containing a list of variants to be filtered.")
    parser.add_argument("-o", "--outfile", required=True, help="Output file for filtered variants.")
    parser.add_argument("-r", "--refgenome", required=True, help="File in FASTA format containing the reference genome.")
    
    args = parser.parse_args()
    filter_homopolymers(args.infile, args.outfile, args.refgenome)
    print("Homopolymer filtering complete.")

if __name__ == '__main__':
    main()