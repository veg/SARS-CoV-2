#!/bin/bash
#fdate=$(date +"%Y-%m-%d")
#fdate=$1
#fdate="2020-08-01_2020-10-31"
#fdate="2021-02-01_2021-03-31"
#fdate="2019-12-01_2020-02-28"
#fdate="2020-01-01_2020-03-31"
fdate="2020-02-01_2020-04-30"
#fdate="2020-03-01_2020-05-31"
#fdate="2020-04-01_2020-06-30"
#fdate="2020-05-01_2020-07-31"
#fdate="2020-06-01_2020-08-31"
#fdate="2020-07-01_2020-09-30"


FQUEUE='priority'
QUEUE='epyc2'
OUTPUT_DIR=`pwd`/logs/$fdate/
#BASE_DIR="/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel"
BASE_DIR=`pwd`
mkdir -p $OUTPUT_DIR

SMALLPPN=8
LARGEPPN=16

# Gene
qsub -l nodes=1:ppn=$LARGEPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $FQUEUE -F "$BASE_DIR/data/sliding-windows/$fdate S $LARGEPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/sliding-windows/$fdate M $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/sliding-windows/$fdate N $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/sliding-windows/$fdate E $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/sliding-windows/$fdate ORF3a $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/sliding-windows/$fdate ORF6 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/sliding-windows/$fdate ORF7a $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/sliding-windows/$fdate ORF7b $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/sliding-windows/$fdate ORF8 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh 
#qsub -l nodes=2:ppn=64 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/sliding-windows/$fdate ORF1a 128" $BASE_DIR/scripts/extract_genes.sh
#qsub -l nodes=2:ppn=64 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIRkdata/sliding-windows/$fdate ORF1b 128" $BASE_DIR/scripts/extract_genes.sh

# Products
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/sliding-windows/$fdate leader $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/sliding-windows/$fdate nsp2 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$LARGEPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $FQUEUE -F "$BASE_DIR/data/sliding-windows/$fdate nsp3 $LARGEPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/sliding-windows/$fdate nsp4 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/sliding-windows/$fdate 3C $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/sliding-windows/$fdate nsp6 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/sliding-windows/$fdate nsp7 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/sliding-windows/$fdate nsp8 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/sliding-windows/$fdate nsp9 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/sliding-windows/$fdate nsp10 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/sliding-windows/$fdate RdRp $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/sliding-windows/$fdate helicase $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/sliding-windows/$fdate exonuclease $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/sliding-windows/$fdate endornase $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/sliding-windows/$fdate methyltransferase $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/sliding-windows/$fdate ORF10 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh
