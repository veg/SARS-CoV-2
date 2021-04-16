import csv
import json
import sys
import argparse
import itertools
import shutil
import copy
from datetime import date, timedelta
from operator import itemgetter
from Bio import SeqIO

# If one fails, then copy the other to the output. If both fail, then throw an error
def fix(k,v):
    old_key = k
    # chop off everything after null
    chopped_seq_name = list(itertools.takewhile(lambda x: x != 'null', k.split('_')))
    new_key = '_'.join(chopped_seq_name) + '_null_' + str(len(v.values()))
    transform = lambda x: x if x != old_key else new_key
    new_vals = {k: transform(v) for k,v in v.items()}
    return (new_key, new_vals, old_key)

def fix_duplicates(duplicates, map, overwrite):

    print(duplicates)

    with open(duplicates, 'r') as duplicates_fh:
        dupe_json = json.load(duplicates_fh)

    itemize = [fix(*x) for x in dupe_json.items()]
    output_json = {item[0]: item[1] for item in itemize}
    old_name_map = {item[2]: item[0] for item in itemize}

    # Overwrite files
    orig_fn = duplicates

    tmp_fn = orig_fn + '.tmp'

    if map:
        with open(map, 'w') as map_fh:
            json.dump(old_name_map, map_fh, indent=4, sort_keys=True)

    if overwrite:
        # Fix FASTA headers
        with open(tmp_fn, 'w') as tmp_fp:
            json.dump(output_json, tmp_fp, indent=4, sort_keys=True)

        shutil.move(tmp_fn, orig_fn)

    else:
        print(output_json)


if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Report which dates have full report')
    arguments.add_argument('-d', '--duplicates',   help = 'fasta to filter duplicates', required = True, type = str)
    arguments.add_argument('-m', '--map', help='output map to old sequence names', type=str, default=None)
    arguments.add_argument('-o', '--overwrite', help='overwrite duplicate file, otherwise write to stdout', action='store_true')
    args = arguments.parse_args()
    fix_duplicates(args.duplicates, args.map, args.overwrite)

