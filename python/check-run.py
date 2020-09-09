from os import path
from itertools import product
import json
import numpy as np
import collections
import csv
import multiprocessing
from multiprocessing import Pool
from datetime import date, timedelta

# gene list
genes=['leader','nsp2','nsp3','nsp4','3C','nsp6','nsp7','nsp8','nsp9','nsp10','helicase','exonuclease','endornase','S','E','M','N','ORF3a','ORF6','ORF7a','ORF8','RdRp','methyltransferase']

# Specify dates
dates = [date(2020, 9, 1)]
gene_fn = lambda x,y: path.join(basedir, x.strftime('%Y-%m-%d'), 'sequences.' + y + '.json')

# get directory listing
basedir = 'data/fasta/'

def collect_info(item):
    gene = gene_fn(*item)
    return (gene, path.exists(gene))

def main():
    cpus = multiprocessing.cpu_count()
    combos = list(product(dates, genes))
    exists = [collect_info(x) for x in combos]
    print(list(filter(lambda x: not x[1], exists)))

main()
