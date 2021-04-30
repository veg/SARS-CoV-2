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
# dates = [date(2021, 2, 22)]

# sliding_windows = [
#                     ("2020-08-01", "2020-10-31"),
#                     ("2021-02-01", "2021-03-31"),
#                     ("2019-12-01", "2020-02-28"),
#                     ("2020-01-01", "2020-03-31"),
#                     ("2020-02-01", "2020-04-30"),
#                     ("2020-03-01", "2020-05-31"),
#                     ("2020-04-01", "2020-06-30"),
#                     ("2020-05-01", "2020-07-31"),
#                     ("2020-06-01", "2020-08-31"),
#                     ("2020-07-01", "2020-09-30")
#                   ]

#dirs = ['_'.join(window) for window in sliding_windows]
dirs = ['2021-04-23']

gene_fn = lambda x,y: path.join(basedir, x, 'sequences.' + y + '.json')

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
    combos = list(product(dirs, genes))
    exists = [collect_info(x) for x in combos]
    print(len(list(filter(lambda x: not x[1] or not x[2], exists))))
    print(list(filter(lambda x: not x[1] or not x[2], exists)))

main()

