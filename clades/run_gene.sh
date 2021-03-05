#!/usr/bin/bash
#PBS -N run_gene
#PBS -o /home/aglucaci/SARS-CoV-2/clades/STDOUT
#PBS -e /home/aglucaci/SARS-CoV-2/clades/STDOUT

module load aocc


# parameter passed from submit_jobs.sh
GENE=$1

# Working-directory
BASEDIR="/home/aglucaci/SARS-CoV-2/clades"

# Input fasta
#INPUT=$BASEDIR"/TEST_DATA/gisaid_romania.fasta"

# South Africa Clade B.1.351
#INPUT=$BASEDIR"/B-1-351/gisaid_hcov-19_2021_02_07_17.fasta.fa"

# Bronx clade
#INPUT=$BASEDIR"/Bronx/gisaid_hcov-19_2021_02_16_20_Bronx.fasta.fa"

# New York clade
INPUT=$BASEDIR"/NYC/gisaid_hcov-19_2021_03_04_19.fasta.fa"


# Reference data
#REF_ALN=$BASEDIR"/REFERENCE_DATA/CODON/sequences.${GENE}.compressed.fas"
REF_ALN=$BASEDIR"/REFERENCE_DATA/ViPR_ReferenceSet/sequences.${GENE}.compressed.fas"
REF_SEQ=$BASEDIR"/REFERENCE_DATA/reference_genes/${GENE}.fas"

# Custom label for query sequences
#LABEL="ROM"
#LABEL="SouthAfrica_B1_351"
#LABEL="BRONX"
LABEL="NYC"

# Set where the config is
CONFIG=$BASEDIR"/config.json"

# Output directory
#OUTPUT_DIR=$BASEDIR"/TEST_ROM"
#OUTPUT_DIR=$BASEDIR"/SouthAfrica_B1_351"
#OUTPUT_DIR=$BASEDIR"/Bronx"
OUTPUT_DIR=$BASEDIR"/NYC"

# Create output directory
mkdir -p $OUTPUT_DIR

#python3 run.py -A --max_reference 200 --max_query 500 --threshold_query 0.0001 -i input/B.1.351.fasta  -O gisaid/sequences.${GENE}.compressed.fas -r ref/${GENE}.fa -o output_B1351 -L "B.1.351" -P config_mp.json  -g $GENE

#echo 
python3 $BASEDIR"/run.py" -A --max_reference 400 --max_query 1000 --threshold_query 0.0001 -i $INPUT  -O $REF_ALN -r $REF_SEQ -o $OUTPUT_DIR -L $LABEL -P $CONFIG  -g $GENE


# End of file
