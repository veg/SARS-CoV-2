import glob
import shutil
import datetime
from os import path
from itertools import product
import json
import numpy as np
import itertools
import collections
import csv
import warnings
from datetime import date, timedelta

warnings.filterwarnings('ignore')

'''
Data dictionary

1. analysis_date, Analysis date
2. gene, Gene/segment name
3. num_seqs, Number of sequences in the analysis
4. num_unique_haps, Number of unique haplotypes in the analysis
5. VPH, Median # of variants per haplotypes (VPH)
6. VPH_Q_1, 1st quartile of VPH
7. VPH_Q_3, 3rd quartile of VPH
8. tree_len_internal, Tree length internal branches
9. tree_len_terminal, Tree length terminal branches
10. mean_dnds_leaf, Mean dN/dS on leaf branches
11. mean_dnds_internal, Mean dN/dS on internal branches
12. num_sites, Number of sites in the alignments
13. codon_num_var_sites, Number of variable sites in the alignment (codon-level)
14. codon_num_var_sites_minor, Number of sites with 2 or more sequences that have the same minor variant (codon-level)
15. prot_num_var_sites, Number of variable sites in the alignment (protein-level)
16. prot_num_var_sites_minor, Number of sites with 2 or more sequences that have the same minor variant (protein-level)
17. num_sites_pos_fel, Number of sites with significant FEL results (positively selected)
18. num_sites_neg_fel, Number of sites with significant FEL results (negatively selected)
19. num_sites_meme, Number of sites with significant MEME results
20. median_branches_meme, Median # of branches where there were significnt MEME results
'''

gene_dup_fn = lambda x,y: path.join(x, 'sequences.' + y + '.duplicates.json')
gene_meme_fn = lambda x,y: path.join(x, 'sequences.' + y + '.MEME.json')
gene_meme_full_fn = lambda x,y: path.join(x, 'sequences.' + y + '.FULL.MEME.json')
gene_slac_fn = lambda x,y: path.join(x, 'sequences.' + y + '.SLAC.json')
gene_fel_fn = lambda x,y: path.join(x, 'sequences.' + y + '.FEL.json')

gene_bgm_fn = lambda x,y: path.join(x, 'sequences.' + y + '.combined.fas.BGM.json')
gene_prime_fn = lambda x,y: path.join(x, 'sequences.' + y + '.PRIME.json')
gene_full_meme_fn = lambda x,y: path.join(x, 'sequences.' + y + '.FULL.MEME.json')
gene_relax_fn = lambda x,y: path.join(x, 'sequences.' + y + '.RELAX.json')
gene_busted_fn = lambda x,y: path.join(x, 'sequences.' + y + '.BUSTED.json')
gene_absrel_fn = lambda x,y: path.join(x, 'sequences.' + y + '.ABSREL.json')
gene_fade_fn = lambda x,y: path.join(x, 'sequences.' + y + '.FADE.json')
gene_cfel_fn = lambda x,y: path.join(x, 'sequences.' + y + '.CFEL.json')


'''
analysis_date, gene, num_seqs, num_unique_haps, VPH, VPH_Q_1, VPH_Q_3,
tree_len_internal, tree_len_terminal, mean_dnds_leaf,
mean_dnds_internal,codon_num_var_sites, codon_num_var_sites_minor,
prot_num_var_sites, prot_num_var_sites_minor, num_sites_pos_fel,
num_sites_neg_fel, num_sites_meme, median_branches_meme
'''

def get_variant_count(s, min_count=1):
    variants = collections.Counter(s)
    variants = { x: v for x,v in variants.items() if x != '---' }
    to_rtn = { x: v for x,v in variants.items() if v >= min_count }
    return len(to_rtn.keys()) - 1

def collect_info(item):

    pval = 0.1

    meme_fn = gene_meme_fn(*item)
    meme_full_fn = gene_meme_full_fn(*item)
    prime_fn = gene_prime_fn(*item)
    busted_fn = gene_busted_fn(*item)
    relax_fn = gene_relax_fn(*item)
    slac_fn = gene_slac_fn(*item)
    fel_fn = gene_fel_fn(*item)
    cfel_fn = gene_cfel_fn(*item)
    fade_fn = gene_fade_fn(*item)
    bgm_fn = gene_bgm_fn(*item)
    absrel_fn = gene_absrel_fn(*item)

    to_return = {
        "analysis_date": item[0],
        "gene": item[1],
        "num_seqs": None,
        "num_uniq_haps": None,
        "VPH": None,
        "VPH_Q_1": None,
        "VPH_Q_3": None,
        "tree_len_internal": None,
        "tree_len_terminal": None,
        "mean_dnds_leaf": None,
        "mean_dnds_internal": None,
        "num_sites": None,
        "codon_num_var_sites": None,
        "codon_num_var_sites_minor": None,
        "prot_num_var_sites": None,
        "prot_num_var_sites_minor": None,
        "num_sites_pos_fel": None,
        "num_sites_neg_fel": None,
        "num_sites_pos_cfel": None,
        "num_sites_neg_cfel": None,
        "num_sites_directional": None,
        "num_sites_coevolving": None,
        "num_sites_meme": None,
        "num_sites_meme_full": None,
        "num_sites_pos_prime": None,
        "num_sites_neg_prime": None,
        "is_sig_relax": None,
        "is_sig_busted": None,
        "num_branches_absrel": None,
        "median_branches_meme": None
    }

    # If aBSREL
    try:
        with open(absrel_fn) as absrel_fh:
            absrel = json.load(absrel_fh)

        pvals = [ p['Corrected P-value'] for p in absrel["branch attributes"]["0"].values() if p['Corrected P-value'] is not None and p['Corrected P-value'] < pval ]
        to_return["num_branches_absrel"] = len(pvals)

    except:
        print(f'No aBSREL results for : {item}')


    # If RELAX
    try:
        with open(relax_fn) as relax_fh:
            relax = json.load(relax_fh)

        to_return["is_sig_relax"] = relax["test results"]["p-value"] < pval

    except:
        print(f'No RELAX results for : {item}')

    # If BUSTED
    try:
        with open(busted_fn) as busted_fh:
            busted = json.load(busted_fh)

        to_return["is_sig_busted"] = busted["test results"]["p-value"] < pval

    except:
        print(f'No BUSTED results for : {item}')


    # If PRIME
    try:
        with open(prime_fn) as prime_fh:
            prime = json.load(prime_fh)

        sig_sites = [x[1] > x[0] for x in prime["MLE"]["content"]["0"] if x[4] <= pval]
        to_return["num_sites_pos_prime"] = len([x for x in sig_sites if x])
        to_return["num_sites_neg_prime"] = len(sig_sites) - to_return["num_sites_pos_prime"]

    except:
        print(f'No PRIME results for : {item}')


    # If Contrast-FEL
    try:
        with open(cfel_fn) as cfel_fh:
            cfel = json.load(cfel_fh)

        qval = 0.2

        header_key = 'Q-value (overall)'
        qval_index = [ header[0] for header in cfel["MLE"]["headers"]].index('Q-value (overall)')
        sig_sites = [x[2] > x[1] for x in cfel["MLE"]["content"]["0"] if x[qval_index] <= qval]

        to_return["num_sites_pos_cfel"] = len([x for x in sig_sites if x])
        to_return["num_sites_neg_cfel"] = len(sig_sites) - to_return["num_sites_pos_cfel"]

    except:
        print(f'No Contrast-FEL results for : {item}')

    # If FADE
    try:
        with open(fade_fn) as fade_fh:
            fade = json.load(fade_fh)

        # Get BayesFactor Index
        bayes_factor_index = [header[0] for header in fade['MLE']['headers']].index('BayesFactor[bias>0]')
        bayes_threshold = 100
        # flatten content
        items = itertools.chain(*[val['0'] for val in fade['MLE']['content'].values()])
        to_return["num_sites_directional"] = sum([ c[bayes_factor_index] > bayes_threshold for c in items])

    except:
        print(f'No FADE results for : {item}')

    # If BGM
    try:
        with open(bgm_fn) as bgm_fh:
            bgm = json.load(bgm_fh)

        bgm_prob_threshold = 0.5

        # Check if there is MLE content
        if 'MLE' not in bgm.keys():
            to_return["num_sites_coevolving"] = 0
        else:
            # Get all indices that are probabilities, any one of them over 0.5
            items = [ header[0].startswith('P ') for header in bgm["MLE"]["headers"] ]
            count = sum([ any(prob > bgm_prob_threshold for prob in itertools.compress(content, items)) for content in bgm["MLE"]["content"] ])
            to_return["num_sites_coevolving"] = count

    except:
        print(f'No BGM results for : {item}')

    # If MEME
    try:
        with open(meme_fn) as meme_fh:
            meme = json.load(meme_fh)

        # Get branch lengths and node name
        bls = [(k, v['Global MG94xREV']) for k,v in meme['branch attributes']['0'].items()]

        to_return["tree_len_internal"] = sum(map(lambda x: x[1], filter(lambda x: x[0].startswith('Node'), bls)))
        to_return["tree_len_terminal"] = sum(map(lambda x: x[1], filter(lambda x: not x[0].startswith('Node'), bls)))

        to_return["mean_dnds_leaf"] = meme['fits']['Global MG94xREV']['Rate Distributions']['non-synonymous/synonymous rate ratio for *background*'][0][0]
        to_return["mean_dnds_internal"] = meme['fits']['Global MG94xREV']['Rate Distributions']['non-synonymous/synonymous rate ratio for *test*'][0][0]

        to_return["num_sites_meme"] = len([row for row in meme["MLE"]["content"]["0"] if row[6] <= pval])
        to_return["median_branches_meme"] = np.median([row[7] for row in meme["MLE"]["content"]["0"] if row[6] <= pval])

    except:
        print(f'No MEME results for : {item}')

    # If MEME FULL
    try:
        with open(meme_full_fn) as meme_full_fh:
            meme_full = json.load(meme_full_fh)

        # Get branch lengths and node name
        to_return["num_sites_meme_full"] = len([row for row in meme_full["MLE"]["content"]["0"] if row[6] <= pval])

    except:
        print(f'No Full MEME results for : {item}')

    # If SLAC
    try:
        with open(slac_fn) as slac_fh:
            slac = json.load(slac_fh)

        to_return["num_sites"] = slac["input"]["number of sites"]
        to_return["num_seqs"] = slac["input"]["number of sequences"]

        # Construct codon matrix, transpose, and count
        by_site = np.transpose([v["codon"][0] for k,v in slac["branch attributes"]["0"].items()])
        to_return["codon_num_var_sites"] = len([x for x in map(get_variant_count, by_site) if x > 0])
        to_return["codon_num_var_sites_minor"] = len([x for x in map(lambda x: get_variant_count(x,2), by_site) if x > 0])

        by_prot_site = np.transpose([v["amino-acid"][0] for k,v in slac["branch attributes"]["0"].items()])
        to_return["prot_num_var_sites"] = len([x for x in map(get_variant_count, by_prot_site) if x > 0])
        to_return["prot_num_var_sites_minor"] = len([x for x in map(lambda x: get_variant_count(x,2), by_prot_site) if x > 0])

    except:
        print(f'No SLAC results for : {item}')

    # If FEL
    try:
        with open(fel_fn) as fel_fh:
            fel = json.load(fel_fh)

        sig_sites = [x[1] > x[0] for x in fel["MLE"]["content"]["0"] if x[4] <= pval]
        to_return["num_sites_pos_fel"] = len([x for x in sig_sites if x])
        to_return["num_sites_neg_fel"] = len(sig_sites) - to_return["num_sites_pos_fel"]

    except:
        print(f'No FEL results for : {item}')

    # Return dictionary with all items
    return to_return

def rascl_summary_report(input_dir, output_fn):

    # gene list
    genes=['leader','nsp2','nsp3','nsp4','3C','nsp6','nsp7','nsp8','nsp9','nsp10','helicase','exonuclease','endornase','S','E','M','N','ORF3a','ORF6','ORF7a','ORF8','RdRp','methyltransferase']

    # get directory listing
    combos = list(product([input_dir], genes))
    row_items = [ collect_info(combo) for combo in combos ]

    with open(output_fn, 'w', newline='') as csvfile:
        fieldnames = row_items[0].keys()
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row_item in row_items:
            writer.writerow(row_item)

def copy_reports(input_fn):

    basepath = "/data/shares/web/web/covid-19/selection-analyses/rascl/"
    # datetime.datetime.fromtimestamp(1629590400)
    # Get timestamps from directories
    timestamp = int(input_fn.split('/')[-2])
    gene = input_fn.split('/')[-3]
    dt = datetime.datetime.fromtimestamp(timestamp)
    exec_date = datetime.datetime(*dt.timetuple()[:3]) + datetime.timedelta(days=1)
    exec_date_str = exec_date.isoformat() + "+00:00"

    full_path = basepath + gene + '/scheduled__' + exec_date_str + '/report.csv'
    alt_full_path = basepath + gene + '/backfill__' + exec_date_str + '/report.csv'

    try:
        shutil.copyfile(input_fn, full_path, follow_symlinks=True)
    except OSError as err:
        print(err)
        try:
            shutil.copyfile(input_fn, alt_full_path, follow_symlinks=True)
        except OSError as err:
            print(err)
            print("no destination directory present! " + full_path + " or " + alt_full_path)

def main():

    # basedir = '/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/data/rascl/AY.1/1630195200/'
    # output_fn = 'report.csv'
    # rascl_summary_report(basedir, output_fn)

    dirs = glob.glob("/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/rascl/**/*")

    for dir in dirs:
        output_fn = dir + "/report.csv"
        print(dir)
        print(output_fn)
        rascl_summary_report(dir, output_fn)
        copy_reports(output_fn)

#main()
