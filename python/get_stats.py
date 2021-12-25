import csv
import json
import sys
import argparse
import copy
import itertools
import multiprocessing
from datetime import date, timedelta
import datetime
from dateutil.relativedelta import *
from dateutil.rrule import *
from dateutil.parser import *
from operator import itemgetter
from multiprocessing import Pool

import pymongo
from pymongo import MongoClient

def get_sliding_window_counts(sliding_windows):

    db = MongoClient(host='129.32.209.134')

    # All sliding windows
    counts = []

    for sliding_window in sliding_windows:
        print("Processing")
        print(sliding_window)
        start_date = sliding_window[0]
        end_date = sliding_window[1]
        mongo_query = { "seq": {"$exists":True} }
        mongo_query["collected"] = { "$gt": datetime.datetime.strptime(start_date, "%Y-%m-%d"), "$lt": datetime.datetime.strptime(end_date, "%Y-%m-%d") }
        mongo_query["originalCollected"] = { "$regex": "[0-9]{4}-[0-9]{2}" }
        count = db.gisaid.records.find(mongo_query).count();
        counts.append((sliding_window, count))

    return counts


def get_clade_counts():

    db = MongoClient(host='129.32.209.134')
    # Aggregate is too slow, so get distinct lineage then counts

    # All lineages
    lineages = db.gisaid.records.distinct("pangolinLineage")

    counts = []
    for lineage in lineages:
        count = db.gisaid.records.find({"pangolinLineage": lineage }).count()
        counts.append((lineage, count))

    return counts

def get_unique_haplos(gene):
    db = MongoClient(host='129.32.209.134')
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
    # genes = ["leader","nsp2","nsp3","nsp4","3C","nsp6","nsp7","nsp8","nsp9","nsp10","helicase","exonuclease","endornase","S","E","M","N","ORF3a","ORF6","ORF7a","ORF8","RdRp","methyltransferase"]
    # counts = get_all_unique_haplos(genes)

    # Supplement with 3 month sliding windows since beginning of pandemic
    sliding_windows = []
    TODAY = datetime.date.today()
    LASTMONTH = TODAY-relativedelta(months=+1, day=31)
    THREEMONTHSAGO = TODAY-relativedelta(months=+3, day=1)

    starts = [dt.strftime('%Y-%m-%d') for dt in rrule(MONTHLY, interval=1,bymonthday=(1),dtstart=parse("20191201T000000"), until=THREEMONTHSAGO)]
    ends = [dt.strftime('%Y-%m-%d') for dt in rrule(MONTHLY, interval=1,bymonthday=(-1),dtstart=parse("20200228T000000"), until=LASTMONTH)]
    sliding_windows = set(list(zip(starts,ends)) + sliding_windows)
    counts = get_sliding_window_counts(list(sliding_windows)[:1])

    spamwriter = csv.writer(args.output, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    spamwriter.writerow(['begin', 'end', 'count'])
    for x in counts:
        spamwriter.writerow([x[0][0], x[0][1], x[1]])

