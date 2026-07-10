import sys
import csv
csv.field_size_limit(sys.maxsize)
from ast import literal_eval
import argparse
import os


def parse_args():
    parser = argparse.ArgumentParser(description="Convert UnbeQuant output to the custom ProteoBench format.")
    parser.add_argument("-i", "--input", required=True, help="Input TSV file with raw quantification and identifications as given by UnbeQuant.")
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file for protein header information extraction (mapping back to species)")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file in ProteoBench format.")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    input_file = args.input
    fasta_input = args.fasta
    output_file = args.output

    # Handle case where output is a directory
    if os.path.isdir(output_file):
        output_file = os.path.join(output_file, "proteobench_output.tsv")
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Read fasta first and save all header information:
    fasta_headers = {}
    with open(fasta_input, 'r') as fasta_file:
        current_header = None
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):
                protein, descr = line[1:].split("|")[1:]  # Get the first word after '>'
                fasta_headers[protein] = descr  # Store full header without '>'

    # Read UnbeQuant output and save all relevant information and write output
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        all_headers = next(reader)

        intensity_headers = [idx for idx, h in enumerate(all_headers) if h.endswith("_____intensity")]
        intensity_headers_out = [h[:-14] for h in all_headers if h.endswith("_____intensity")]
        charge_headers = [idx for idx, h in enumerate(all_headers) if h.endswith("_____charge")]
        protein_headers = [idx for idx, h in enumerate(all_headers) if h.endswith("l_prot_ident")]
        peptide_headers = [idx for idx, h in enumerate(all_headers) if h.endswith("l_pep_ident")]
        raw_peptide_headers = [idx for idx, h in enumerate(all_headers) if h.endswith("l_raw_pep_ident")]

        # Prepare for table with the following columns:
        # Sequence, Proteins, Charge, Modified sequence, Intensity1, Intensity2, ...
        result_entries = []
        for row in reader:
            # Check if it is identified:
            raw_peps = set([y for idx in raw_peptide_headers if row[idx] for y in literal_eval(row[idx])])

            if len(raw_peps) == 1:  
                # We only include uniquely identified linked features
                # We do not include linked features with multiple different peptide identifications
                # or unidentified linked features in the final table.
                raw_peps = raw_peps.pop()
                # Handle n-terminal modifications to be changed to the PROFORMA format
                if raw_peps.startswith("n"):  # n terminal mod:
                    mod = raw_peps[raw_peps.index("["):raw_peps.index("]")+1]
                    raw_peps = "" + mod + "-" + raw_peps[raw_peps.index("]")+1:]

                sequence = set([y for idx in peptide_headers if row[idx] for y in literal_eval(row[idx])]).pop()
                proteins = set([z for idx in protein_headers if row[idx] for y in literal_eval(row[idx]) for z in y.split(",")])
                charge = set([row[idx] for idx in charge_headers if row[idx]]).pop()

                # Special Case for iRT peptides: Ignore
                if "iRT_kit_WR_fusion" in proteins and len(proteins) == 1:
                    continue
                if "iRT_kit_WR_fusion" in proteins:
                    proteins.remove("iRT_kit_WR_fusion")

                # Get columns correctly formatted:
                proteins = ";".join([p + "_" + fasta_headers[p].split("_")[1] for p in proteins])
                raw_peptide = raw_peps
                
                result_entries.append(
                    [sequence, proteins, charge, raw_peptide] + [row[idx] for idx in intensity_headers]
                )

        # Write header for ProteoBench format
        writer.writerow(["Sequence", "Proteins", "Charge", "Modified sequence"] + intensity_headers_out)
        writer.writerows(result_entries)
