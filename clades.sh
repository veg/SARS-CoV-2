#!/bin/bash
#fdate=$(date +"%Y-%m-%d")
#fdate=$1
#fdate="2020-08-01_2020-10-31"
#fdate="2021-02-01_2021-03-31"
#fdate="2019-12-01_2020-02-28"
#fdate="2020-01-01_2020-03-31"
#fdate="2020-02-01_2020-04-30"
#fdate="2020-03-01_2020-05-31"
#fdate="2020-04-01_2020-06-30"
#fdate="2020-05-01_2020-07-31"
#fdate="2020-06-01_2020-08-31"
#fdate="2020-07-01_2020-09-30"


FQUEUE='epyc2'
QUEUE='epyc'
OUTPUT_DIR=`pwd`/logs/$fdate/
#BASE_DIR="/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel"
BASE_DIR=`pwd`
mkdir -p $OUTPUT_DIR

SMALLPPN=16
LARGEPPN=32

clades=(B.1.2 B.1.596 B.1 B.1.1.519 B.1.243 B.1.234 B.1.526.1 B.1.1 B.1.526.2 B.1.575 R.1 B.1.1.7 B.1.429 B.1.427 B.1.351 P.1 B.1.526 P.2 B.1.525)


for fdate in ${clades[@]};
  do
    echo $fdate;
    # Gene
    qsub -l nodes=1:ppn=$LARGEPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $FQUEUE -F "$BASE_DIR/data/clades/$fdate S $LARGEPPN" $BASE_DIR/scripts/extract_genes.sh;
    qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/clades/$fdate M $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh;
    qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/clades/$fdate N $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh;
    qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/clades/$fdate E $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh;
    qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/clades/$fdate ORF3a $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh;
    qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/clades/$fdate ORF6 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh;
    qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/clades/$fdate ORF7a $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh;
    qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/clades/$fdate ORF7b $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh;
    qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/clades/$fdate ORF8 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh ;
    #qsub -l nodes=2:ppn=64 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/clades/$fdate ORF1a 128" $BASE_DIR/scripts/extract_genes.sh
    #qsub -l nodes=2:ppn=64 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIRkdata/clades/$fdate ORF1b 128" $BASE_DIR/scripts/extract_genes.sh

    # Products
    qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/clades/$fdate leader $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh;
    qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/clades/$fdate nsp2 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh;
    qsub -l nodes=1:ppn=$LARGEPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $FQUEUE -F "$BASE_DIR/data/clades/$fdate nsp3 $LARGEPPN" $BASE_DIR/scripts/extract_genes.sh;
    qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/clades/$fdate nsp4 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh;
    qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/clades/$fdate 3C $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh;
    qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/clades/$fdate nsp6 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh;
    qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/clades/$fdate nsp7 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh;
    qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/clades/$fdate nsp8 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh;
    qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/clades/$fdate nsp9 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh;
    qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/clades/$fdate nsp10 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh;
    qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/clades/$fdate RdRp $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh;
    qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/clades/$fdate helicase $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh;
    qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/clades/$fdate exonuclease $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh;
    qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/clades/$fdate endornase $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh;
    qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/clades/$fdate methyltransferase $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh;
    qsub -l nodes=1:ppn=$SMALLPPN -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "$BASE_DIR/data/clades/$fdate ORF10 $SMALLPPN" $BASE_DIR/scripts/extract_genes.sh;
  done;
