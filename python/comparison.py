import csv
import json
import sys
import argparse
import itertools
from datetime import date, timedelta
from operator import itemgetter

# genes=['leader','nsp2','nsp3','nsp4','3C','nsp6','nsp7','nsp8','nsp9','nsp10','helicase','exonuclease','endornase','S','E','M','N','ORF3a','ORF6','ORF7a','ORF8','RdRp','methyltransferase']
# genes=['nsp7']

arguments = argparse.ArgumentParser(description='Compares sites of significance between two runs')
arguments.add_argument('-i', '--input',   help = 'summary results file', required = True, type = argparse.FileType('r'))
arguments.add_argument('-j', '--input-unfiltered',   help = 'summary results file', required = True, type = argparse.FileType('r'))
arguments.add_argument('-o', '--output', help = 'Write results here', type = argparse.FileType('w'), default = sys.stdout)

args = arguments.parse_args()

gene_meme_fn = lambda x,y: path.join(basedir, x, 'sequences.' + y + '.MEME.json')

# get gene from filename
# add filenames to CSV

file_filtered = json.loads(args.input.read())
file_unfiltered = json.loads(args.input_unfiltered.read())

fn_filtered = args.input.name
fn_unfiltered = args.input_unfiltered.name
gene = fn_filtered.split('.')[1]

# Only pick up sites that are considered under positive selection
# Add p-value

inv_map_filtered = {file_filtered['map'][i] : i for i in range(len(file_filtered['map']))}
inv_map_unfiltered = {file_unfiltered['map'][i] : i for i in range(len(file_unfiltered['map']))}

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

# get inverse map for both

# Get selection sites from first dataset (mapped)
sites_filtered = set([k for k in file_filtered['selection'].keys() if 'kind' in file_filtered['selection'][k].keys() and file_filtered['selection'][k]['kind'] == 'positive'])
sites_filtered = set([file_filtered['map'][int(site)]+1 for site in sites_filtered])

# Get selection sites from second dataset (mapped)
sites_unfiltered = set([k for k in file_unfiltered['selection'].keys() if 'kind' in file_unfiltered['selection'][k].keys() and file_unfiltered['selection'][k]['kind'] == 'positive'])
sites_unfiltered = set([file_unfiltered['map'][int(site)]+1 for site in sites_unfiltered])

all_keys = sites_filtered.union(sites_unfiltered)

to_return = [{ "gene": gene, "fel_filtered": get_pval(file_filtered, x, 'fel'), "fel_unfiltered": get_pval(file_unfiltered, x, 'fel'), "meme_filtered": get_pval(file_filtered, x, 'meme'), "meme_unfiltered": get_pval(file_unfiltered, x, 'meme'), "fn_filtered": fn_filtered, "fn_unfiltered": fn_unfiltered,"site" : x, "in_filtered": x in sites_filtered , "in_unfiltered": x in sites_unfiltered , "in_both": x in sites_filtered and x in sites_unfiltered } for x in all_keys]

# Print intersection by creating dict for all keys and which set they are in
fieldnames = ['gene', 'fn_filtered', 'fn_unfiltered', 'site', 'in_filtered', 'in_unfiltered', 'in_both', 'fel_filtered', 'fel_unfiltered', 'meme_filtered', 'meme_unfiltered']
writer = csv.DictWriter(args.output, fieldnames=fieldnames)
writer.writeheader()
for row in to_return:
    writer.writerow(row)


