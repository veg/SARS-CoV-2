#!/usr/bin/zsh

source $HOME/.zshrc

fdate=$(date +"%Y-%m-%d")
fdate="2020-04-27"
#echo "Getting metadata"
#node /data/shares/veg/SARS-CoV-2/gisaid-downloader/get_metadata.js
#mv /data/shares/veg/SARS-CoV-2/gisaid-downloader/downloads/* /data/shares/gisaid/
#echo "Getting sequences"
#node /data/shares/veg/SARS-CoV-2/gisaid-downloader/get_seqs.js
#mv /data/shares/veg/SARS-CoV-2/gisaid-downloader/downloads/gisaid_cov2020_sequences.fasta  /data/shares/gisaid/$fdate.fasta
#echo "Translating to pipeline format"
#node /data/shares/veg/SARS-CoV-2/gisaid-downloader/translate-to-master.js /data/shares/gisaid/$fdate-metadata.json -o /data/shares/gisaid/$fdate.master.json -f /data/shares/gisaid/$fdate.fasta
#echo "Generating file with no seqs"
#jq 'del(.[].sequence)' /data/shares/gisaid/$fdate.master.json > /data/shares/gisaid/$fdate.master.nofasta.json
/data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/preprocess.sh
cp /data/shares/gisaid/$fdate.master.json /data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate/master.json
#cp /data/shares/gisaid/$fdate.master.nofasta.json /data/shares/web/web/covid-19/$fdate.master.nofasta.json
#ln -sf /data/shares/web/web/covid-19/$fdate.master.nofasta.json  /data/shares/web/web/covid-19/latest.json 
