#!/bin/bash
#fdate=$(date +"%Y-%m-%d")
#fdate=$1
#fdate="2020-11-18"
fdate="2020-10-12"

FQUEUE='epyc'
QUEUE='epyc'
OUTPUT_DIR=`pwd`/logs/$fdate/
#BASE_DIR="/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel"
BASE_DIR=`pwd`
mkdir $OUTPUT_DIR

# Gene
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $FQUEUE -F "$BASE_DIR/data/fasta/$fdate S 16" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate M 8" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate N 8" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate E 8" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate ORF3a 8" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate ORF6 8" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate ORF7a 8" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate ORF7b 8" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate ORF8 8" $BASE_DIR/scripts/extract_genes.sh 
#qsub -l nodes=2:ppn=64 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate ORF1a 128" $BASE_DIR/scripts/extract_genes.sh
#qsub -l nodes=2:ppn=64 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIRkdata/fasta/$fdate ORF1b 128" $BASE_DIR/scripts/extract_genes.sh

# Products
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate leader 8" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate nsp2 8" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $FQUEUE -F "$BASE_DIR/data/fasta/$fdate nsp3 16" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate nsp4 8" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate 3C 8" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate nsp6 8" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate nsp7 8" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate nsp8 8" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate nsp9 8" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate nsp10 8" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate RdRp 8" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate helicase 8" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate exonuclease 8" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate endornase 8" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate methyltransferase 8" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate ORF10 8" $BASE_DIR/scripts/extract_genes.sh

