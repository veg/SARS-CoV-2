import csv
import json
import collections

with open('./data/ctl/predictions.csv', newline='') as csvfile:
    reader = csv.DictReader(csvfile, delimiter=',', quotechar='"')
    rows = [row for row in reader]

# Only get peptides of length 9
epitopes = list(filter(lambda x: len(x["Peptide"]) == 9, rows))

epis = {epitope["Peptide"]: {} for epitope in epitopes}

for epi in epitopes:
    key = "HLA-" + epi["HLA Allele"]
    item = {key : float(format(float(epi["Best Score (nM)"]), '.3f'))}
    epis[epi["Peptide"]].update(item)


with open('./data/ctl/epitopes.json', 'w') as outfile:
    json.dump(epis, outfile, indent=4, separators=(',', ': '))
