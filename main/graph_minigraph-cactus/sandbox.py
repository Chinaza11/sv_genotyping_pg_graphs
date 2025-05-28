# -*- coding: utf-8 -*-
"""
Created on Wed May 28 10:51:09 2025

@author: Chinaza
"""

import pandas as pd
import os
import sys

# =============================================================================
# thoroughly examining the filtered VCF files
# =============================================================================

def parse_vcf_to_dataframe(vcf_file_path):
    """
    Parses a VCF file to create a Pandas DataFrame with
    CHROM, POS, ID, REF, ALT, MAX_ALT, and CLASSIFICATION columns.

    Args:
        vcf_file_path (str): Path to the input VCF file.

    Returns:
        pd.DataFrame: A DataFrame containing the parsed data, or None if an error occurs.
    """
    data_rows = []

    try:
        with open(vcf_file_path, 'r') as infile:
            for line in infile:
                # Skip header lines
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')

                # Ensure line has at least the first 5 standard VCF columns
                if len(fields) < 5:
                    print(f"Warning: Skipping malformed line in {vcf_file_path}: {line.strip()}", file=sys.stderr)
                    continue

                chrom = fields[0]
                pos = int(fields[1]) # Convert POS to integer
                # ID can be '.' if not specified, keep as string
                id_val = fields[2] if fields[2] != '.' else None
                ref_seq = fields[3]
                alt_seqs_str = fields[4]

                ref_len = len(ref_seq)
                max_alt_len = 0
                classification = None # Initialize classification

                alt_alleles = alt_seqs_str.split(',')

                for alt_allele in alt_alleles:
                    alt_len = len(alt_allele)
                    if alt_len > max_alt_len:
                        max_alt_len = alt_len

                # Determine classification based on indel criteria
                if max_alt_len > ref_len:
                    classification = "insertion"
                elif max_alt_len < ref_len:
                    classification = "deletion"
                # If max_alt_len == ref_len, it's likely a SNP or complex substitution,
                # and per the original awk logic, we only classify true indels.
                # So 'classification' remains None for these.

                data_rows.append({
                    'CHROM': chrom,
                    'POS': pos,
                    'ID': id_val,
                    'REF': ref_seq,
                    'ALT': alt_seqs_str, # Keep original ALT string for context
                    'NO_of_ALT_ALLELE': len(alt_alleles),
                    'REF_LEN': ref_len,
                    'MAX_ALT_LEN': max_alt_len,
                    'CLASSIFICATION': classification
                })

    except FileNotFoundError:
        print(f"Error: Input VCF file not found at {vcf_file_path}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An error occurred while processing {vcf_file_path}: {e}", file=sys.stderr)
        return None

    df = pd.DataFrame(data_rows)
    return df

# --- Example Usage ---
if __name__ == "__main__":
    # Replace 'path/to/your/input.vcf' with the actual path to your VCF file
    # For demonstration, let's create a dummy VCF file
    dummy_vcf_content = """##fileformat=VCFv4.2
##fileDate=20240528
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	rs123	A	G	100	PASS	.
chr1	101	.	T	TAA	90	PASS	.
chr1	105	indel1	CAT	C	85	PASS	.
chr2	200	.	G	GCCTCCTGGACTT	70	PASS	.
chr2	205	.	AAAAA	A	60	PASS	.
chr3	300	.	T	A,TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	50	PASS	.
chr4	400	.	CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	C	45	PASS	.
chr5	500	snp1	C	G	100	PASS	.
chr6	600	long_indel	AAAAA	A,GCCTTAGGCCTTAGGCCTTAGGCCTTAGGCCTTAGGCCTTAGGCCTTAGGCCTTAGGCCTTAGGCCTTAGG	95	PASS	.
chrX	700	.	AGCT	A	70	PASS	.
"""
    dummy_vcf_path = "dummy_input.vcf"
    with open(dummy_vcf_path, 'w') as f:
        f.write(dummy_vcf_content)

    print(f"Processing VCF file: {dummy_vcf_path}")
    vcf_df = parse_vcf_to_dataframe(dummy_vcf_path)

    if vcf_df is not None:
        print("\n--- Full DataFrame ---")
        print(vcf_df)

        print("\n--- Filtered to show only Classified Indels ---")
        # Filter for rows where classification is not None
        indels_df = vcf_df[vcf_df['CLASSIFICATION'].notna()]
        print(indels_df)

        # Clean up dummy file
        os.remove(dummy_vcf_path)
        
        

vcf_df = parse_vcf_to_dataframe("C:/Users/nnamd/Downloads/Grif16309.filtered.vcf")

# examine function, maybe delete later
vcf_file_path = "C:/Users/nnamd/Downloads/Grif16309.filtered.vcf"

vcf_df_2 = parse_vcf_to_dataframe("C:/Users/nnamd/Downloads/Grif16309.filtered_2.vcf")

# =============================================================================
# 
# =============================================================================

def reverse_sequence(sequence):
  """
  Reverses the order of a given DNA sequence string.

  Args:
    sequence (str): The input DNA sequence.

  Returns:
    str: The reversed DNA sequence.
  """
  return sequence[::-1]

# Your sequence
original_sequence = "ATACTACCACTAGTCAAAACACTAGGGAATCCTACTAGGAGACATTTAGTAGGAAAAGTTTATACT"

# Get the reversed sequence
reversed_seq = reverse_sequence(original_sequence)

print(f"Original Sequence: {original_sequence}")
print(f"Reversed Sequence: {reversed_seq}")