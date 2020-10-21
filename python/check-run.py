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
#genes=['leader','nsp2','nsp4','3C','nsp6','nsp7','nsp8','nsp9','nsp10','helicase','exonuclease','endornase','E','M','ORF3a','ORF6','ORF7a','ORF8','methyltransferase']

# Specify dates
dates = [date(2020, 9, 15), date(2020, 9, 22), date(2020, 10, 5), date(2020, 10, 12)]
#dates = [date(2020, 9, 22), date(2020, 9, 15)]
gene_fn = lambda x,y: path.join(basedir, x.strftime('%Y-%m-%d'), 'sequences.' + y + '.json')

# get directory listing
basedir = 'data/fasta/'

def collect_info(item):
    gene = gene_fn(*item)
    # path exists and is not size 0
    if path.exists(gene):
        return (gene, True, path.getsize(gene))
    else:
        return (gene, False, 0)

def main():
    cpus = multiprocessing.cpu_count()
    combos = list(product(dates, genes))
    exists = [collect_info(x) for x in combos]
    print(list(filter(lambda x: not x[1] or not x[2], exists)))

main()
