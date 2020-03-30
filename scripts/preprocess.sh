#!/bin/bash

fdate=$(date +"%Y-%m-%d")
python3 python/extract-attributes.py -j /data/shares/gisaid/$fdate.master.json > /data/shares/gisaid/$fdate.attributes.csv
mkdir -p data/fasta/$fdate/
cp /data/shares/gisaid/$fdate.attributes.csv data/fasta/$fdate/attributes.csv
cp /data/shares/gisaid/$fdate.master.nofasta.json data/fasta/$fdate/master-no-fasta.json
cp /data/shares/gisaid/$fdate.master.json data/fasta/$fdate/master.json
python3 python/extract-sequences.py -j /data/shares/gisaid/$fdate.master.json -f "host" "re" "[hH]uman" -f "sequence" ">" 28000 > data/fasta/$fdate/sequences
