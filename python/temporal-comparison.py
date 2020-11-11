from os import path
from itertools import product
from datetime import date, timedelta
from operator import itemgetter
import csv
import json
import sys
import argparse
import itertools
import multiprocessing


basedir = 'data/fasta/'

genes=['leader','nsp2','nsp3','nsp4','3C','nsp6','nsp7','nsp8','nsp9','nsp10','helicase','exonuclease','endornase','S','E','M','N','ORF3a','ORF6','ORF7a','ORF8','RdRp','methyltransferase']

sdate = date(2020, 4, 1)
edate = date(2020, 9, 1)
delta = edate - sdate
dates = [(sdate + timedelta(days=i)).strftime('%Y-%m-%d') for i in range(delta.days + 1)]

gene_summary_fn = lambda x,y: path.join(basedir, x, 'sequences.' + y + '.json')

# get gene from filename
# add filenames to CSV

# Only pick up sites that are considered under positive selection
# Add p-value

# inv_map_filtered = {file_filtered['map'][i] : i for i in range(len(file_filtered['map']))}
# inv_map_unfiltered = {file_unfiltered['map'][i] : i for i in range(len(file_unfiltered['map']))}

def get_pval(results, site, method):
    try:
        inv_map = {results['map'][i] : str(i-1) for i in range(len(results['map']))}
        site = inv_map[site]
        if site in results['selection'].keys():
            return results['selection'][site][method]
        else:
            return None
    except:
        print('error for ' + str(site) + ' ' + method + ' for ' + fn_filtered)

def collect_info(item):

    summary_fn = gene_summary_fn(*item)

    # get inverse map for both
    with open(summary_fn) as summary_fh:
        summary = json.load(summary_fh)

    # Get selection sites from first dataset (mapped)
    sites = set([k for k in summary['selection'].keys() if 'kind' in summary['selection'][k].keys() and summary['selection'][k]['kind'] == 'positive'])
    sites = set([summary['map'][int(site)]+1 for site in sites])

    # to_return = [{ "gene": gene, "fel_filtered": get_pval(file_filtered, x, 'fel'), "fel_unfiltered": get_pval(file_unfiltered, x, 'fel'), "meme_filtered": get_pval(file_filtered, x, 'meme'), "meme_unfiltered": get_pval(file_unfiltered, x, 'meme'), "fn_filtered": fn_filtered, "fn_unfiltered": fn_unfiltered,"site" : x, "in_filtered": x in sites_filtered , "in_unfiltered": x in sites_unfiltered , "in_both": x in sites_filtered and x in sites_unfiltered } for x in all_keys]
    return { "date" : item[0], "gene" : item[1], "sites": list(sites)}


def main():
    cpus = multiprocessing.cpu_count()
    combos = list(product(dates, genes))

    # filter combos based on available dates
    # only use dates for files that exists
    valid_combos = list(filter(lambda x: path.exists(gene_summary_fn(*x)), combos))
    row_items = [collect_info(x) for x in valid_combos]

    output = sys.stdout

    # with Pool(cpus) as p:
    #     row_items = p.map(collect_info, combos)

    # Print intersection by creating dict for all keys and which set they are in
    # fieldnames = ['gene', 'site', **dates]
    fieldnames = ['date', 'gene', 'sites']
    writer = csv.DictWriter(output, fieldnames=fieldnames)
    writer.writeheader()
    for row in row_items:
        writer.writerow(row)

main()
