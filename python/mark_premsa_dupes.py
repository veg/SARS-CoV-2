import csv
import json
import sys
import argparse
import itertools
import shutil
import copy
import os
import multiprocessing
from multiprocessing import Pool
from datetime import date, timedelta
from operator import itemgetter
from Bio import SeqIO

from pymongo import MongoClient, InsertOne, DeleteMany, ReplaceOne, UpdateOne, UpdateMany


def update_dupe_records(dupes, reference, gene):

    # nuc_seq_str = str(nuc_seq.seq)
    # prot_seq_str = str(prot_seq.seq)

    # db.restaurant.updateMany(
    #    { violations: { $gt: 4 } },
    #    { $set: { "Review" : true } }
    # );

    to_update = {}

    nuc_seq_key = gene + '_premsa_nuc_seq'
    if(nuc_seq_key in reference.keys()):
        to_update[nuc_seq_key] = reference[gene + '_premsa_nuc_seq']

    prot_seq_key = gene + '_premsa_protein_seq'
    if(prot_seq_key in reference.keys()):
        to_update[prot_seq_key] = reference[gene + '_premsa_protein_seq']


    # print(epi_id)
    qc_key = ".".join(["qc", gene])
    if "qc" in reference.keys():
        if gene in reference["qc"].keys():
            qc_value = reference["qc"][gene]
        else:
            qc_value = {"passed": True }
    else:
        qc_value = {"passed": True }

    dupe_key = ".".join(["duplicate_of_by_gene", gene])

    to_update[qc_key] = qc_value
    to_update[dupe_key] = reference["id"]
    return UpdateMany({'id': {"$in" : dupes}}, {'$set': to_update})

def get_references(ids, db):
    query = {"id" : {"$in" : ids} }
    return list(db.gisaid.records.find(query))

def mark_premsa_dupes(dupe_input, gene):

    db = MongoClient(host='192.168.0.4')
    dupes = json.loads(open(dupe_input, 'r').read())

    # Shave meta
    shave_meta = lambda x: '_'.join(x.split('_')[:3])

    # Transform ids to only include epi_isl_xxxxxxx
    fmt_dupes = { shave_meta(ref): list(map(shave_meta, dupe.values()))[1:]  for ref, dupe in dupes.items() if len(dupe) > 1}

    # Get documents for all references
    ref_keys = list(fmt_dupes.keys())
    ref_docs = { ref["id"]: ref for ref in get_references(ref_keys, db) }

    # key references by id
    get_update_queries = [update_dupe_records(v, ref_docs[k], gene) for k,v in fmt_dupes.items()]

    results = db.gisaid.records.bulk_write(get_update_queries)
    print("Updated " +  str(results.modified_count) + " of " + str(sum([len(v) for k,v in fmt_dupes.items()])) + " items to update")


if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Mark duplicates in MongoDB')
    arguments.add_argument('-d', '--dupe-input',   help = 'fasta to update', required = True, type = str)
    arguments.add_argument('-t', '--type',   help = 'gene region', required = True, type = str)
    args = arguments.parse_args()
    mark_premsa_dupes(args.dupe_input, args.type)

