#!/bin/bash
clear

genes=(leader nsp2 nsp3 nsp4 3C nsp6 nsp7 nsp8 nsp9 nsp10 helicase exonuclease endornase S E M N ORF3a ORF6 ORF7a ORF8 RdRp methyltransferase) 

BASEDIR="/home/aglucaci/SARS-CoV-2/clades"

# TEST_ROM
#FASTA=$BASEDIR"/TEST_DATA/gisaid_romania.fasta"

# South-Africa clade B.1.351
#FASTA=$BASEDIR"/B-1-351/gisaid_hcov-19_2021_02_07_17.fasta"

# Bronx sequences
FASTA=$BASEDIR"/Bronx/gisaid_hcov-19_2021_02_16_20_Bronx.fasta"


# Format FASTA headers
echo "# Pre-processing FASTAs"

awk '{ if ($0 ~ "^[^>]") {a = gensub(/[^acgt]/, "", "g"); print a;} else print;}' $FASTA > $FASTA".fa" 
#awk '{ if ($0 ~ "^>") {sub(" ", "_"); print ;} else print;}' => replace ' ' with '_' in FASTA
sed 's, ,_,g' -i $FASTA
awk '{ if ($0 ~ "^>") {b=gensub(/>(.+)\|(EPI_ISL_|epi_isl_)([0-9]+)\|(.+)/, ">epi_isl_\\3/\\1","g"); print b;} else print;}' $FASTA > $FASTA".fa"
#exit 0

# Submit Jobs
echo "# Clade specific analyses"
#echo ""
#genes=(S)

for i in "${genes[@]}"; 
do 
    GENE=${i} 
    echo "Processing: "$GENE 
    echo qsub -l nodes=1:ppn=8,walltime=72:00:00 -q epyc -d `pwd` -F "$GENE" run_gene.sh 
    qsub -l nodes=1:ppn=8,walltime=72:00:00 -q epyc -d `pwd` -F "$GENE" run_gene.sh 

done


# End of file
