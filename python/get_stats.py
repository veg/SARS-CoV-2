import csv
import json
import sys
import argparse
import copy
import itertools
import multiprocessing
from datetime import date, timedelta
from operator import itemgetter
from multiprocessing import Pool

import pymongo
from pymongo import MongoClient

def get_clade_counts():

    db = MongoClient(host='192.168.0.4')
    # Aggregate is too slow, so get distinct lineage then counts

    # All lineages
    lineages = db.gisaid.records.distinct("pangolinLineage")

    counts = []
    for lineage in lineages:
        count = db.gisaid.records.find({"pangolinLineage": lineage }).count()
        counts.append((lineage, count))

    return counts

def get_unique_haplos(gene):
    db = MongoClient(host='192.168.0.4')
    reference_key = 'duplicate_of_by_gene.' + gene
    count = db.gisaid.records.find({reference_key:'reference'}).count()
    return count

def get_all_unique_haplos(genes):
    counts = map(get_unique_haplos, genes)
    return list(zip(genes, counts))

if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Report which dates have full report')
    # arguments.add_argument('-g', '--gene', help = 'gene to get stats with', required = True, type = str)
    arguments.add_argument('-o', '--output', help = 'write compressed fasta here', type = argparse.FileType('w'), default = sys.stdout)
    args = arguments.parse_args()
    genes = ["leader","nsp2","nsp3","nsp4","3C","nsp6","nsp7","nsp8","nsp9","nsp10","helicase","exonuclease","endornase","S","E","M","N","ORF3a","ORF6","ORF7a","ORF8","RdRp","methyltransferase"]
    counts = get_all_unique_haplos(genes)

    spamwriter = csv.writer(args.output, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    for x in counts:
        spamwriter.writerow(x)

