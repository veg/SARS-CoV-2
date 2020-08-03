#!/bin/bash
#PBS -l walltime=99:0:0:0

export PATH=/usr/local/bin:$PATH
source /etc/profile.d/modules.sh
module load aocc/1.3.0
module load openmpi/gnu/3.0.2
export LC_ALL="en_US.UTF-8"
export LC_CTYPE="en_US.UTF-8"

DIRECTORY=$1
FILE=$1/sequences
ATTRIBUTES=$1/attributes.csv
LOG=$1/log.txt
MASTER=$1/master.json
MASTERNOFASTA=$1/master-no-fasta.json
GENE=$2
NP=$3
HYPHY=/data/shares/veg/SARS-CoV-2/hyphy/hyphy
HYPHYMPI=/data/shares/veg/SARS-CoV-2/hyphy/HYPHYMPI
HYPHYLIBPATH=/data/shares/veg/SARS-CoV-2/hyphy/res
MAFFT=/usr/local/bin/mafft
RAXML=/usr/local/bin/raxml-ng-mpi
TN93=/usr/local/bin/tn93
#PREMSA=/Users/sergei/Development/hyphy-analyses/codon-msa/pre-msa.bf
#POSTMSA=/Users/sergei/Development/hyphy-analyses/codon-msa/post-msa.bf
PREMSA=/data/shares/veg/SARS-CoV-2/hyphy-analyses/codon-msa/pre-msa.bf
POSTMSA=/data/shares/veg/SARS-CoV-2/hyphy-analyses/codon-msa/post-msa.bf
WORKING_DIR=/data/shares/veg/SARS-CoV-2/SARS-CoV-2/
PYTHON=/data/shares/veg/SARS-CoV-2/SARS-CoV-2/env/bin/python3

ZERO_LENGTHS_FLAGS='--kill-zero-lengths No ENV="_DO_TREE_REBALANCE_=1"'


function run_a_gene {


GENE=$1
REFERENCE_SEQUENCE=$2
TRIM_FROM=$3
TRIM_TO=$4
N_FRAC=$5

if [ -s ${FILE}.${GENE}.json ]
then 
   echo "$GENE alignment already processed"
else
    genes=(leader nsp2 nsp3 nsp4 3C nsp6 nsp7 nsp8 nsp9 nsp10 helicase exonuclease endornase  S E M N ORF3a ORF6 ORF7a ORF8 RdRp methyltransferase)
    offsets=(265 805 2719 8554 10054 10972 11842 12091 12685 13024 16236 18039 19620 21562 26244 26522 28273 25392 27201 27393 27893 13440 20658)
    fragments=(ORF1a ORF1a ORF1a ORF1a ORF1a ORF1a ORF1a ORF1a ORF1a ORF1a ORF1b ORF1b ORF1b S E M N ORF3a ORF6 ORF7a ORF8 ORF1b ORF1b)
    shifts=(0        180   818   2763  3263  3569  3859  3942  4140  4253  922   1523  2050  0 0 0 0 0     0    0     0    -10   2396)
    add_one=(0       0      0     0     0     0     0     0     0     0    1     1     1     0 0 0 0 0     0    0     0    1    1)

		GENE_INDEX=99

		for i in "${!genes[@]}"; do
			 if [[ "${genes[$i]}" = "${GENE}" ]]; then
					 GENE_INDEX=${i}
			 fi
		done

    FRAGMENT=${fragments[$GENE_INDEX]}
    ADDSHIFT=${add_one[$GENE_INDEX]}
    SHIFT=${shifts[$GENE_INDEX]}
    OFFSET=${offsets[$GENE_INDEX]}
    ANNOTATION=${FILE}.annotation.json
    cp data/comparative-annotation.json ${ANNOTATION}

    echo "$PYTHON $WORKING_DIR/python/summarize-gene.py -T data/ctl/epitopes.json -D $MASTERNOFASTA -d ${FILE}.${GENE}.duplicates.json -s ${FILE}.${GENE}.SLAC.json -f ${FILE}.${GENE}.FEL.json -m ${FILE}.${GENE}.MEME.json -P 0.1 --output  ${FILE}.${GENE}.json -c ${FILE}.${GENE}.compressed.fas -E data/evo_annotation.json -A data/mafs.csv -V data/evo_freqs.csv -F $FRAGMENT --frame_shift ${ADDSHIFT} --fragment_shift $SHIFT -S $OFFSET -O $ANNOTATION"
    $PYTHON $WORKING_DIR/python/summarize-gene.py -T data/ctl/epitopes.json -D $MASTERNOFASTA -d ${FILE}.${GENE}.duplicates.json -s ${FILE}.${GENE}.SLAC.json -f ${FILE}.${GENE}.FEL.json -m ${FILE}.${GENE}.MEME.json -P 0.1 --output  ${FILE}.${GENE}.json -c ${FILE}.${GENE}.compressed.fas -E data/evo_annotation.json -A data/mafs.csv -V data/evo_freqs.csv -F $FRAGMENT --frame_shift ${ADDSHIFT} --fragment_shift $SHIFT -S $OFFSET -O $ANNOTATION

fi

}

case $GENE in

  S)
    echo "Analyzing S Gene"
    run_a_gene "S" "$WORKING_DIR/reference_genes/S.fas" "20000" "27000" 0.005
    ;;

  M)
    echo "Analyzing M gene"
    run_a_gene "M" "$WORKING_DIR/reference_genes/M.fas" "25000" "30000" 0.01
    ;;

  N)
    echo "Analyzing N gene"
    run_a_gene "N" "$WORKING_DIR/reference_genes/N.fas" "26000" "35000" 0.01
    ;;

  E)
    echo "Analyzing E gene"
    run_a_gene "E" "$WORKING_DIR/reference_genes/E.fas" "25500" "27000" 0.01
    ;;

  ORF3a)
    echo "Analyzing ORF3a gene"
    run_a_gene "ORF3a" "$WORKING_DIR/reference_genes/ORF3a.fas" "24000" "27000" 0.01
    ;;

  ORF6)
    echo "Analyzing ORF6 gene"
    run_a_gene "ORF6" "$WORKING_DIR/reference_genes/ORF6.fas" "26000" "30000" 0.01
    ;;

  ORF7a)
    echo "Analyzing ORF7a gene"
    run_a_gene "ORF7a" "$WORKING_DIR/reference_genes/ORF7a.fas" "26000" "35000" 0.01
    ;;

  ORF7b)
    echo "Analyzing ORF7b gene"
    run_a_gene "ORF7b" "$WORKING_DIR/reference_genes/ORF7b.fas" "27500" "28000" 0.01
    ;;

  ORF8)
    echo "Analyzing ORF8 gene"
    run_a_gene "ORF8" "$WORKING_DIR/reference_genes/ORF8.fas" "26000" "35000" 0.01
    ;;

  ORF10)
    echo "Analyzing ORF10 gene"
    run_a_gene "ORF10" "$WORKING_DIR/reference_genes/ORF10.fas" "29500" "29800" 0.01
    ;;

  ORF1a)
    echo "Analyzing ORF1a gene"
    run_a_gene "ORF1a" "$WORKING_DIR/reference_genes/ORF1a.fas" "1" "15000" 0.001
    ;;

  ORF1b)
    echo "Analyzing ORF1b gene"
    run_a_gene "ORF1b" "$WORKING_DIR/reference_genes/ORF1b.fas" "12000" "24000" 0.001
    ;;

  leader)
    echo "Analyzing leader product"
    run_a_gene "leader" "$WORKING_DIR/reference_genes/leader.fas" "1" "1000" 0.001
    ;;

  nsp2)
    echo "Analyzing nsp2 product"
    run_a_gene "nsp2" "$WORKING_DIR/reference_genes/nsp2.fas" "500" "3000" 0.001
    ;;

  nsp3)
    echo "Analyzing nsp3 product"
    run_a_gene "nsp3" "$WORKING_DIR/reference_genes/nsp3.fas" "2000" "10000" 0.001
    ;;

  nsp4)
    echo "Analyzing nsp4 product"
    run_a_gene "nsp4" "$WORKING_DIR/reference_genes/nsp4.fas" "8000" "11000" 0.001
    ;;

  3C)
    echo "Analyzing 3C product"
    run_a_gene "3C" "$WORKING_DIR/reference_genes/3C.fas" "9000" "12000" 0.001
    ;;

  nsp6)
    echo "Analyzing nsp6 product"
    run_a_gene "nsp6" "$WORKING_DIR/reference_genes/nsp6.fas" "10500" "12500" 0.001
    ;;

  nsp7)
    echo "Analyzing nsp7 product"
    run_a_gene "nsp7" "$WORKING_DIR/reference_genes/nsp7.fas" "11500" "12500" 0.001
    ;;

  nsp8)
    echo "Analyzing nsp8 product"
    run_a_gene "nsp8" "$WORKING_DIR/reference_genes/nsp8.fas" "11500" "13000" 0.001
    ;;

  nsp9)
    echo "Analyzing nsp9 product"
    run_a_gene "nsp9" "$WORKING_DIR/reference_genes/nsp9.fas" "12000" "13500" 0.001
    ;;

  nsp10)
    echo "Analyzing nsp10 product"
    run_a_gene "nsp10" "$WORKING_DIR/reference_genes/nsp10.fas" "12500" "14000" 0.001
    ;;

  RdRp)
    echo "Analyzing RdRp product"
    run_a_gene "RdRp" "$WORKING_DIR/reference_genes/RdRp.fas" "13000" "17000" 0.001
    ;;

  helicase)
    echo "Analyzing helicase product"
    run_a_gene "helicase" "$WORKING_DIR/reference_genes/helicase.fas" "15500" "18500" 0.001
    ;;

  exonuclease)
    echo "Analyzing exonuclease product"
    run_a_gene "exonuclease" "$WORKING_DIR/reference_genes/exonuclease.fas" "17500" "20000" 0.001
    ;;

  endornase)
    echo "Analyzing endornase product"
    run_a_gene "endornase" "$WORKING_DIR/reference_genes/endornase.fas" "19000" "21000" 0.001
    ;;

  methyltransferase)
    echo "Analyzing methyltransferase product"
    run_a_gene "methyltransferase" "$WORKING_DIR/reference_genes/methyltransferase.fas" "20000" "22000" 0.001
    ;;


  *)
    echo -n "Unknown Gene, exiting"
    ;;
esac
