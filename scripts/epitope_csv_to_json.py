import csv
import json
import collections

with open('./data/ctl/predictions.csv', newline='') as csvfile:
    reader = csv.DictReader(csvfile, delimiter=',', quotechar='"')
    rows = [row for row in reader]

# Read in Nelde
with open('./data/ctl/nelde_epitopes.csv', newline='') as csvfile:
    reader = csv.DictReader(csvfile, delimiter=',', quotechar='"')
    nelde_rows = [row for row in reader]

# Only get peptides of length 9
epitopes = list(filter(lambda x: len(x["Peptide"]) == 9, rows))
nelde_epitopes = list(filter(lambda x: len(x["Sequence"]) >= 9, nelde_rows))
epis = {epitope["Peptide"]: {} for epitope in epitopes}

nelde_epis = {nepitope["Sequence"]: {} for nepitope in nelde_epitopes}

for epi in epitopes:
    key = "HLA-" + epi["HLA Allele"]
    gene = epi["hlagene"]
    item = { 'type' : 'predicted', 'source' : 'boni', 'gene' : gene }
    score = { key : float(format(float(epi["Best Score (nM)"]), '.3f')) }
    if "score" in epis[epi["Peptide"]].keys():
        score_dict = epis[epi["Peptide"]]["score"]
    else :
        score_dict = {}
    score_dict.update(score)
    item["score"] = score_dict
    epis[epi["Peptide"]].update(item)

for epi in nelde_epitopes:
    key = epi["hlafam"]
    gene = epi["hlagene"]
    item = { 'gene' : gene, 'type' : 'experimentally validated', 'source': 'nelde'}
    score = { key : None }
    item["score"] = score
    nelde_epis[epi["Sequence"]].update(item)

epis.update(nelde_epis)

with open('./data/ctl/epitopes.json', 'w') as outfile:
    json.dump(epis, outfile, indent=4, separators=(',', ': '))
