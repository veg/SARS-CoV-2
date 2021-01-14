#!/bin/bash 

DATE=$(date +"%Y-%m-%d")
mkdir data/fasta/$DATE
python3 python/export-sequences.py -o data/fasta/$DATE/sequences
python3 python/export-meta.py -o data/fasta/$DATE/master-no-fasta.json
./submit_jobs.sh $DATE
