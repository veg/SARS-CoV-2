import csv
import json
import sys
import argparse
import itertools
from datetime import date, timedelta
from operator import itemgetter
from Bio import SeqIO
import copy

# If one fails, then copy the other to the output. If both fail, then throw an error
def merge(protein_json, nuc_json):
    # For each key in nuc_json, get all values. If one of the values is in protein duplicates, then merge all values in the protein duplicate to nucleotide keys, and remove values from protein duplicates.

    # Get all keys for nuc_json
    protein_vals = [list(values.values()) for values in protein_json.values()]
    new_nuc_json = copy.deepcopy(nuc_json)

    # For each nuc_json, get intersection, if matches, take protein seqs and place them in key and remove
    for k, v in nuc_json.items():
        nuc_seq_names = ['_'.join(x.split('_')[:3]) for x in v.values()]

        for prot_vals in protein_vals:
            prot_vals_trim = ['_'.join(x.split('_')[:3]) for x in prot_vals]
            prot_vals_trim_dict = {'_'.join(x.split('_')[:3]) : x for x in prot_vals}

            # Need to collect *all* nuc_seq_names and intersect with prot_vals_trim
            if(set(nuc_seq_names).intersection(set(prot_vals_trim))):
                # Add copy number count
                if(len(set(prot_vals_trim))):
                    # Get difference between nuc_seq_names and prot_vals, then add diff to nuc_seq
                    difference = set(prot_vals_trim).difference(set(nuc_seq_names))
                    # Get greatest number in nuc_json, and increment from there
                    cnt = max(int(x) for x in new_nuc_json[k].keys()) + 1
                    for x in difference:
                        # Update dictionary if it exists
                        new_nuc_json[k][str(cnt)] = prot_vals_trim_dict[x]
                        cnt+=1

    # Validate that all duplicates count up to original sequence count
    return new_nuc_json

def merge_duplicates(protein_duplicates, nuc_duplicates, output):
    try:
        with open(protein_duplicates, "r") as protein_duplicates_fh:
            protein_json = json.load(protein_duplicates_fh)
    except:
        protein_json = False
    try:
        with open(nuc_duplicates, "r") as protein_duplicates_fh:
            nuc_json = json.load(nuc_duplicates_fh)
    except:
        nuc_json = False

    if not protein_json and not nuc_json:
        output_json = {}
    elif not protein_json:
        output_json = nuc_json
    elif not nuc_json:
        output_json = protein_json
    else:
        output_json = merge(protein_json, nuc_json)

    with open(output, 'w') as output_fh:
        json.dump(output_json, output_fh, indent=4, sort_keys=True)

    return

if __name__ == "__main__":

    arguments = argparse.ArgumentParser(description='Report which dates have full report')
    arguments.add_argument('-p', '--protein-duplicates',   help = 'fasta to filter duplicates', required = True, type = str)
    arguments.add_argument('-n', '--nuc-duplicates',   help = 'duplicates file', required = True, type = str)
    arguments.add_argument('-o', '--output', help = 'write compressed fasta here', type = str, default = sys.stdout)
    args = arguments.parse_args()

    merge_duplicates(args.protein_duplicates, args.nuc_duplicates, args.output)


