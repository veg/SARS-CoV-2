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
from Bio import SeqIO
from collections import Counter
import itertools

#genes = ['leader','nsp2','nsp3','nsp4','3C','nsp6','nsp7','nsp8','nsp9','nsp10','helicase','exonuclease','endornase','S','E','M','N','ORF3a','ORF6','ORF7a','ORF8','RdRp','methyltransferase']
genes = ['S']
dates = ['2020-09-01-bk']

gene_meme_fn = lambda x,y: path.join(basedir, x, 'sequences.' + y + '.MEME.json')
gene_msa_fn = lambda x,y: path.join(basedir, x, 'sequences.' + y + '.msa')
gene_summ_fn = lambda x,y: path.join(basedir, x, 'sequences.' + y + '.json')

# get directory listing
basedir = 'data/fasta/'

def get_minor_cnt(cnt):
    # use counter, remove largest, remove gap, then count
    # remove largest
    cnt.pop(cnt.most_common()[0][0])
    # remove gaps, if any
    if '-' in cnt.keys() : cnt.pop('-')
    return (len(cnt.keys()),sum(cnt.values()))


def get_variance_at_site(item):

    # get mapping
    summ_fn = gene_summ_fn(*item)
    with open(summ_fn) as summ_fh:
        summ = json.load(summ_fh)

    fas_fn = gene_msa_fn(*item)
    msa = np.array([list(seq.seq) for seq in SeqIO.parse(fas_fn, 'fasta')]).transpose()
    cnts = [[item[1], summ['map'][int(i)]+1, *get_minor_cnt(Counter(msa[i]))] for i in range(len(msa))]
    # for each site, get variant composition
    # cnt_dist = [Counter(x) for x in msa]
    return cnts

def collect_info(item):

    meme_fn = gene_meme_fn(*item)

    try:
        with open(meme_fn) as meme_fh:
            memes = json.load(meme_fh)
    except:
        return []

    return memes

def main():

    cpus = multiprocessing.cpu_count()
    combos = list(product(dates, genes))

    ## Linear instead of pool
    cnts = [get_variance_at_site(combo) for combo in combos]
    row_items = [collect_info(combo) for combo in combos]

    # with Pool(cpus) as p:
    #     cnts = p.map(get_variance_at_site, combos)

    # with Pool(cpus) as p:
    #     row_items = p.map(collect_info, combos)


    pvals = []
    for i in range(len(row_items)):
        seq_pvals = [c[6] for c in row_items[i]["MLE"]["content"]["0"]]
        pvals.append(list(zip(seq_pvals, cnts[i])))

    all_pvals = list(itertools.chain(*pvals))

    with open('memes-report.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        headers = ["pval", "gene", "index", "num_minor_variants", "total_num_variants"]
        writer.writerow(headers)
        for row_item in all_pvals:
            writer.writerow([row_item[0], *row_item[1]])


main()
