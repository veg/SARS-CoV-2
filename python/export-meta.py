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

import pymongo
from pymongo import MongoClient

arguments = argparse.ArgumentParser(description='Report which dates have full report')
arguments.add_argument('-o', '--output',   help = 'json output', type = argparse.FileType('w'))

args = arguments.parse_args()
db = MongoClient()

records = list(db.gisaid.records.find({}))
for rec in records:
    rec.pop("_id")
    if "seq" in rec:
        rec.pop("seq")

output_json = { row["id"] : row for row in records }

# Write attributes.csv file and MASTER-NO-JSON file
json.dump(output_json, args.output, indent=4, sort_keys=True)

