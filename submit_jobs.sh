#!/bin/bash
fdate=$(date +"%Y-%m-%d")
#fdate=$1
#fdate="2021-02-23"

FQUEUE='epyc2'
QUEUE='epyc'
OUTPUT_DIR=`pwd`/logs/$fdate/
#BASE_DIR="/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel"
BASE_DIR=`pwd`
mkdir -p $OUTPUT_DIR

SMALLPPN=32
LARGEPPN=64

# Gene
qsub -l nodes=1:ppn=$LARGEPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $FQUEUE -F "$BASE_DIR/data/fasta/$fdate S $LARGEPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate M $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate N $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate E $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate ORF3a $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate ORF6 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate ORF7a $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate ORF7b $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate ORF8 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh 
#qsub -l nodes=2:ppn=64 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate ORF1a 128" $BASE_DIR/scripts/extract_genes.sh
#qsub -l nodes=2:ppn=64 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIRkdata/fasta/$fdate ORF1b 128" $BASE_DIR/scripts/extract_genes.sh

# Products
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate leader $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate nsp2 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$LARGEPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $FQUEUE -F "$BASE_DIR/data/fasta/$fdate nsp3 $LARGEPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate nsp4 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate 3C $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate nsp6 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate nsp7 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate nsp8 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate nsp9 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate nsp10 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate RdRp $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate helicase $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate exonuclease $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate endornase $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate methyltransferase $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/fasta/$fdate ORF10 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
