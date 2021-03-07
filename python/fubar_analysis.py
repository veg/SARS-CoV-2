import os
from os import path
import json
import numpy as np
import collections
import csv
import multiprocessing

from scipy.spatial import distance
from multiprocessing import Pool
from datetime import date, timedelta
from itertools import chain, tee, combinations, product
from Bio import SeqIO
from collections import Counter
import itertools

# Can you write a  script that scrapes FUBAR results (just the weight vectors, each should be 400 numbers) and compute https://en.wikipedia.org/wiki/Cosine_similarity between them all?

genes = ['leader','nsp2','nsp3','nsp4','3C','nsp6','nsp7','nsp8','nsp9','nsp10','helicase','exonuclease','endornase','S','E','M','N','ORF3a','ORF6','ORF7a','ORF8','RdRp','methyltransferase']
#genes = ['leader']
dates = ['2020-09-01']

gene_fubar_fn = lambda x,y: path.join(basedir, x, 'sequences.' + y + '.compressed.fas.FUBAR.json')

# get directory listing
basedir = 'data/fasta/'

def get_cosine_similarity(item):

    # get mapping
    summ_fn = gene_summ_fn(*item)
    with open(summ_fn) as summ_fh:
        summ = json.load(summ_fh)

    fas_fn = gene_msa_fn(*item)
    msa = np.array([list(seq.seq) for seq in SeqIO.parse(fas_fn, 'fasta')]).transpose()
    cnts = [[item[1], summ['map'][int(i)]+1, *get_minor_cnt(Counter(msa[i]))] for i in range(len(msa))]
    # for each site, get variant composition
    return cnts

def collect_info(item):

    fubar_fn = gene_fubar_fn(*item)

    try:
        with open(fubar_fn) as fubar_fh:
            fubars = json.load(fubar_fh)
    except:
        return []

    weights = [x[2] for x in fubars['grid']]
    return (item[1], weights)

def main():

    cpus = multiprocessing.cpu_count()
    combos = list(product(dates, genes))

    ## Linear instead of pool
    row_items = [collect_info(combo) for combo in combos]
    pairwise = list(combinations(row_items, 2))

    # compute cosign similarity
    distances = [(pair[0][0], pair[1][0], distance.cosine(pair[0][1], pair[1][1])) for pair in pairwise]


    with open('fubars-report.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        headers = ["gene_one", "gene_two", "cosine similarity"]
        writer.writerow(headers)
        for row_item in distances:
            writer.writerow(row_item)


main()
