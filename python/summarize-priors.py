import os
from os import path
from itertools import product
import json
import numpy as np
import collections
import csv
import multiprocessing
from multiprocessing import Pool
from datetime import date, timedelta
from itertools import chain

'''
Data dictionary

1. analysis_date, Analysis date
2. gene, Gene/segment name
3. site, Number of sites in the analysis
'''
# gene list
genes = ['leader','nsp2','nsp3','nsp4','3C','nsp6','nsp7','nsp8','nsp9','nsp10','helicase','exonuclease','endornase','S','E','M','N','ORF3a','ORF6','ORF7a','ORF8','RdRp','methyltransferase']
offsets = {'leader':1, 'ORF1b' : 1, 'RdRp' : 1, 'helicase' : 1, 'endornase' : 1, 'exonuclease' : 1, 'ORF3a' : 1, '3C': 1, 'M': 1, 'S':1, 'nsp10':1, 'nsp2':1, 'nsp3':1, 'nsp7':1, 'nsp8':1, 'E':1, 'N':1};

# Specify dates
sdate = date(2020, 3, 30)
edate = date(2020, 9, 8)
#edate = date.today() - timedelta(days = 1)
delta = edate - sdate
dates = [(sdate + timedelta(days=i)).strftime('%Y-%m-%d') for i in range(delta.days + 1)]

gene_sum_fn = lambda x,y: path.join(basedir, x, 'sequences.' + y + '.json')

# get directory listing
basedir = 'data/fasta/'

def collect_info(item):

    pval = 0.1

    sum_fn = gene_sum_fn(*item)
    entries = []

    try:
        with open(sum_fn) as sum_fh:
            sums = json.load(sum_fh)
    except:
        print(f'No sum results for : {item}')
        return entries

    offset = 0
    if item[1] in offsets.keys():
        offset = offsets[item[1]]

    # Get sites with selection
    sites = list(sums['selection'].keys())
    sites = [sums['map'][int(x)]+1 for x in sites]
    entries = [(item[0],item[1],s) for s in sites]

    return entries

def main():

    cpus = multiprocessing.cpu_count()
    combos = list(product(dates, genes))

    #combos = [('2020-09-01', 'leader')]
    #row_items = [collect_info(combos[0])]
    #print(row_items)

    row_items = [collect_info(combo) for combo in combos]

    # with Pool(cpus) as p:
    #     row_items = p.map(collect_info, combos)

    row_items = filter(lambda x : len(x) > 0, row_items)
    row_items = list(chain.from_iterable(row_items))

    with open('priors-report.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        headers = ["date", "region", "site"]
        writer.writerow(headers)
        for row_item in row_items:
            writer.writerow(row_item)


main()
