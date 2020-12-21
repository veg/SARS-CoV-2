#!/bin/bash
#PBS -l walltime=72:0:0:0

export PATH=/usr/local/bin:$PATH
source /etc/profile.d/modules.sh
module load aocc/1.3.0
module load openmpi/gnu/3.0.2
export LC_ALL=en_US
export LC_CTYPE=en_US

DIRECTORY=$1
FILE=$1/sequences
#ATTRIBUTES=$1/attributes.csv
LOG=$1/log.txt
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

WORKING_DIR=/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel
PYTHON=/data/shares/veg/SARS-CoV-2/SARS-CoV-2-devel/env/bin/python3

COMPRESSOR=$WORKING_DIR/scripts/compressor.bf
COMPRESSOR2=$WORKING_DIR/scripts/compressor-2.bf

ZERO_LENGTHS_FLAGS='--kill-zero-lengths Constrain ENV="_DO_TREE_REBALANCE_=1"'

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
    echo "EXTRACTING $GENE"
    if [ -s ${FILE}.${GENE}_protein.fas ] 
    then
        echo "Already extracted"
    else
        TMP_FILE=${FILE}.${GENE}.tmp
        cp ${FILE} ${TMP_FILE}
        echo "mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH $PREMSA --input $TMP_FILE --reference $REFERENCE_SEQUENCE --trim-from $TRIM_FROM --trim-to $TRIM_TO --E 0.01 --N-fraction $N_FRAC --remove-stop-codons Yes"
        mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH $PREMSA --input $TMP_FILE --reference $REFERENCE_SEQUENCE --trim-from $TRIM_FROM --trim-to $TRIM_TO --E 0.01 --N-fraction $N_FRAC --remove-stop-codons Yes
        mv ${TMP_FILE}_protein.fas ${FILE}.${GENE}_protein.fas
        mv ${TMP_FILE}_nuc.fas ${FILE}.${GENE}_nuc.fas
    fi

    echo "INITIAL COMPRESSION $GENE"
    if [ -s ${FILE}.${GENE}_protein.compressed.fas ] && [ -s ${FILE}.${GENE}_nuc.compressed.fas ]
    then
        echo "Already extracted"
    else
        echo "$PYTHON python/get-raw-duplicates.py -i ${FILE}.${GENE}_protein.fas -m ${FILE}.${GENE}_nuc.fas -o ${FILE}.${GENE}_protein.compressed.fas -n ${FILE}.${GENE}_nuc.compressed.fas -d ${FILE}.${GENE}_protein.duplicates.json -g ${FILE}.${GENE}_raw_nucleotide.duplicates.json"
        $PYTHON python/get-raw-duplicates.py -i ${FILE}.${GENE}_protein.fas -m ${FILE}.${GENE}_nuc.fas -o ${FILE}.${GENE}_protein.compressed.fas -n ${FILE}.${GENE}_nuc.compressed.fas -d ${FILE}.${GENE}_protein.duplicates.json -g ${FILE}.${GENE}_raw_nucleotide.duplicates.json
    fi

    echo "ALIGNING PROTEIN DATA"
    if [ -s ${FILE}.${GENE}.msa ] 
    then
        echo "Already aligned"
    else
        echo "$MAFFT --auto --thread -1 --addfragments ${FILE}.${GENE}_protein.compressed.fas reference_genes/reference.${GENE}_protein.fas >| ${FILE}.${GENE}.tmp.msa"
        $MAFFT --auto --thread -1 --addfragments ${FILE}.${GENE}_protein.compressed.fas reference_genes/reference.${GENE}_protein.fas >| ${FILE}.${GENE}.tmp.msa 
        # Remove reference from output
        $PYTHON python/remove-seq.py -i ${FILE}.${GENE}.tmp.msa -r reference_genes/reference.${GENE}_protein.fas -o ${FILE}.${GENE}.msa
        rm ${FILE}.${GENE}.tmp.msa 
    fi
        
    if [ -s ${FILE}.${GENE}.compressed.fas ] 
    then
        echo "Already reverse translated"
    else
        echo "$HYPHY LIBPATH=$HYPHYLIBPATH $POSTMSA --protein-msa ${FILE}.${GENE}.msa --nucleotide-sequences ${FILE}.${GENE}_nuc.fas --output ${FILE}.${GENE}.compressed.fas --duplicates ${FILE}.${GENE}_nucleotide.duplicates.json"
        $HYPHY LIBPATH=$HYPHYLIBPATH $POSTMSA --protein-msa ${FILE}.${GENE}.msa --nucleotide-sequences ${FILE}.${GENE}_nuc.fas --output ${FILE}.${GENE}.compressed.fas --duplicates ${FILE}.${GENE}_nucleotide.duplicates.json
        #Replace all unknown characters with N
        sed -i '/^>/! s/[^ACTG-]/N/g' ${FILE}.${GENE}.compressed.fas
    fi

    echo "MERGING DUPLICATES FROM POST-MSA and PRIOR"
    if [ -s ${FILE}.${GENE}.duplicates.json ] 
    then
        echo "Already MERGED"
    else
        echo "$PYTHON python/merge-duplicates.py -p ${FILE}.${GENE}_raw_nucleotide.duplicates.json -n ${FILE}.${GENE}_nucleotide.duplicates.json -o ${FILE}.${GENE}.duplicates.json"
        $PYTHON python/merge-duplicates.py -p ${FILE}.${GENE}_raw_nucleotide.duplicates.json -n ${FILE}.${GENE}_nucleotide.duplicates.json -o ${FILE}.${GENE}.duplicates.json
        # Fix duplicates 
        $PYTHON python/fix-duplicates.py -d ${FILE}.${GENE}.duplicates.json -m ${FILE}.${GENE}.map.json -o
        # Fix header files
        $PYTHON python/update-fasta-duplicates.py -f ${FILE}.${GENE}.compressed.fas -m ${FILE}.${GENE}.map.json
    fi

    echo "FILTERING ALIGNMENT"
    if [ -s ${FILE}.${GENE}.compressed.filtered.fas ] 
    then
        echo "ALIGNMENT ALREADY FILTERED"
    else

        echo "$HYPHY LIBPATH=$HYPHYLIBPATH $COMPRESSOR --msa ${FILE}.${GENE}.compressed.fas --duplicates ${FILE}.${GENE}.duplicates.json --output ${FILE}.${GENE}.variants.csv --json ${FILE}.${GENE}.variants.json"
        $HYPHY LIBPATH=$HYPHYLIBPATH $COMPRESSOR --msa ${FILE}.${GENE}.compressed.fas --duplicates ${FILE}.${GENE}.duplicates.json --output ${FILE}.${GENE}.variants.csv --json ${FILE}.${GENE}.variants.json

        echo "$HYPHY LIBPATH=$HYPHYLIBPATH $COMPRESSOR2 --msa ${FILE}.${GENE}.compressed.fas --duplicates ${FILE}.${GENE}.duplicates.json --csv ${FILE}.${GENE}.variants.csv --byseq ${FILE}.${GENE}.variants.json --p 0.9 --output ${FILE}.${GENE}.compressed.filtered.fas --json ${FILE}.${GENE}.filtered.json"
        $HYPHY LIBPATH=$HYPHYLIBPATH $COMPRESSOR2 --msa ${FILE}.${GENE}.compressed.fas --duplicates ${FILE}.${GENE}.duplicates.json --csv ${FILE}.${GENE}.variants.csv --byseq ${FILE}.${GENE}.variants.json --p 0.9 --output ${FILE}.${GENE}.compressed.filtered.fas --json ${FILE}.${GENE}.filtered.json
    fi

    SEQCOUNT=$(grep -c ">" ${FILE}.${GENE}.compressed.fas)

    #if [ "$SEQCOUNT" -lt "1000" ]
    #then
    #  ZERO_LENGTHS_FLAGS=""
    #fi

    #if [ -s ${FILE}.${GENE}.tn93 ] 
    #then
    #    echo "Already computed TN93"
    #else
    #    $TN93 -q -t 0.05 ${FILE}.${GENE}.compressed.filtered.fas > ${FILE}.${GENE}.tn93 2> ${FILE}.${GENE}.tn93.json
    #    echo $PYTHON $WORKING_DIR/python/tabulate-diversity-divergence.py -j $MASTER -t ${FILE}.${GENE}.tn93 > $DIRECTORY/evolution.${GENE}.csv
    #    $PYTHON $WORKING_DIR/python/tabulate-diversity-divergence.py -j $MASTER -t ${FILE}.${GENE}.tn93 > $DIRECTORY/evolution.${GENE}.csv
    #fi

    if [ -s ${FILE}.${GENE}.compressed.filtered.fas.rapidnj.bestTree ] 
    then
        echo "Already has tree"
    else
        echo "seqmagick convert ${FILE}.${GENE}.compressed.filtered.fas ${FILE}.${GENE}.compressed.filtered.sto"
        echo "rapidnj ${FILE}.${GENE}.compressed.filtered.sto -i sth > ${FILE}.${GENE}.compressed.filtered.fas.rapidnj.bestTree"
        seqmagick convert ${FILE}.${GENE}.compressed.filtered.fas ${FILE}.${GENE}.compressed.filtered.sto
        rapidnj ${FILE}.${GENE}.compressed.filtered.sto -i sth > ${FILE}.${GENE}.compressed.filtered.fas.rapidnj.bestTree
        sed -i "s/'//g" ${FILE}.${GENE}.compressed.filtered.fas.rapidnj.bestTree
        #$RAXML --tree pars{5} --msa ${FILE}.${GENE}.compressed.filtered.fas --threads 4 --model GTR+G --force
    fi


    #if [ -s ${FILE}.${GENE}.compressed.filtered.fas.raxml.bestTree ] 
    #then
    #    echo "Already has tree"
    #else
    #    $RAXML --tree pars{5} --msa ${FILE}.${GENE}.compressed.filtered.fas --threads 4 --model GTR+G --force
    #fi
    
    if [ -s ${FILE}.${GENE}.SLAC.json ] 
    then
        echo "Already has SLAC results"
    else
        echo mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH slac $ZERO_LENGTHS_FLAGS --alignment ${FILE}.${GENE}.compressed.filtered.fas --tree ${FILE}.${GENE}.compressed.filtered.fas.rapidnj.bestTree --branches All --samples 0 --output ${FILE}.${GENE}.SLAC.json
        mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH slac $ZERO_LENGTHS_FLAGS --alignment ${FILE}.${GENE}.compressed.filtered.fas --tree ${FILE}.${GENE}.compressed.filtered.fas.rapidnj.bestTree --branches All --samples 0 --output ${FILE}.${GENE}.SLAC.json
    fi

    if [ -s ${FILE}.${GENE}.FEL.json ] 
    then
        echo "Already has FEL results"
    else
        echo mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH fel $ZERO_LENGTHS_FLAGS --alignment ${FILE}.${GENE}.compressed.filtered.fas --tree ${FILE}.${GENE}.compressed.filtered.fas.rapidnj.bestTree --branches Internal --output ${FILE}.${GENE}.FEL.json
        mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH fel $ZERO_LENGTHS_FLAGS --alignment ${FILE}.${GENE}.compressed.filtered.fas --tree ${FILE}.${GENE}.compressed.filtered.fas.rapidnj.bestTree --branches Internal --output ${FILE}.${GENE}.FEL.json
    fi

    if [ -s ${FILE}.${GENE}.MEME.json ] 
    then
        echo "Already has MEME results"
    else
        echo mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH meme $ZERO_LENGTHS_FLAGS --alignment ${FILE}.${GENE}.compressed.filtered.fas --tree ${FILE}.${GENE}.compressed.filtered.fas.rapidnj.bestTree --branches Internal --output ${FILE}.${GENE}.MEME.json
        mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH meme $ZERO_LENGTHS_FLAGS --alignment ${FILE}.${GENE}.compressed.filtered.fas --tree ${FILE}.${GENE}.compressed.filtered.fas.rapidnj.bestTree --branches Internal --output ${FILE}.${GENE}.MEME.json
    fi

    #if [ -s ${FILE}.${GENE}.FUBAR.json ] 
    #then
    #    echo "Already has FUBAR results"
    #else
    #   echo "$HYPHY LIBPATH=$HYPHYLIBPATH  fubar $ZERO_LENGTHS_FLAGS --grid 40 --alignment ${FILE}.${GENE}.compressed.filtered.fas --tree ${FILE}.${GENE}.compressed.filtered.fas.rapidnj.bestTree --output ${FILE}.${GENE}.FUBAR.json"
    #   $HYPHY LIBPATH=$HYPHYLIBPATH  fubar $ZERO_LENGTHS_FLAGS --grid 40 --alignment ${FILE}.${GENE}.compressed.filtered.fas --tree ${FILE}.${GENE}.compressed.filtered.fas.rapidnj.bestTree --output ${FILE}.${GENE}.FUBAR.json
    #fi

    #if [ -s ${FILE}.${GENE}.PRIME.json ] 
    #then
    #    echo "Already has PRIME results"
    #else
    #    echo mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH prime --alignment ${FILE}.${GENE}.compressed.filtered.fas --tree ${FILE}.${GENE}.compressed.filtered.fas.rapidnj.bestTree --branches Internal --output ${FILE}.${GENE}.PRIME.json
    #    mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH prime --alignment ${FILE}.${GENE}.compressed.filtered.fas --tree ${FILE}.${GENE}.compressed.filtered.fas.rapidnj.bestTree --branches Internal --output ${FILE}.${GENE}.PRIME.json
    #fi

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

    echo "$PYTHON $WORKING_DIR/python/summarize-gene.py -T data/ctl/epitopes.json -B data/single_mut_effects.csv -D $MASTERNOFASTA -d ${FILE}.${GENE}.duplicates.json -s ${FILE}.${GENE}.SLAC.json -f ${FILE}.${GENE}.FEL.json -m ${FILE}.${GENE}.MEME.json -P 0.1 --output  ${FILE}.${GENE}.json -c ${FILE}.${GENE}.compressed.filtered.fas -E data/evo_annotation.json -A data/mafs.csv -V data/evo_freqs.csv -F $FRAGMENT --frame_shift ${ADDSHIFT} --fragment_shift $SHIFT -S $OFFSET -O $ANNOTATION"
    $PYTHON $WORKING_DIR/python/summarize-gene.py -T data/ctl/epitopes.json -B data/single_mut_effects.csv -D $MASTERNOFASTA -d ${FILE}.${GENE}.duplicates.json -s ${FILE}.${GENE}.SLAC.json -f ${FILE}.${GENE}.FEL.json -m ${FILE}.${GENE}.MEME.json -P 0.1 --output  ${FILE}.${GENE}.json -c ${FILE}.${GENE}.compressed.filtered.fas -E data/evo_annotation.json -A data/mafs.csv -V data/evo_freqs.csv -F $FRAGMENT --frame_shift ${ADDSHIFT} --fragment_shift $SHIFT -S $OFFSET -O $ANNOTATION

    #if [ -s ${FILE}.${GENE}.BGM.json ] 
    #then
    #    echo "Already has BGM results"
    #else
    #    echo mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH bgm --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.rapidnj.bestTree --branches Internal --min-subs 2 --type codon
    #    mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH bgm --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.rapidnj.bestTree --branches Internal --min-subs 2 --type codon
    #fi

    #if [ -s ${FILE}.${GENE}.BUSTED.json ] 
    #then
    #    echo "Already has BUSTED results"
    #else
    #    echo $HYPHY LIBPATH=$HYPHYLIBPATH busted --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.rapidnj.bestTree --branches Internal --output ${FILE}.${GENE}.BUSTED.json --rates 2 --syn-rates 2 --starting-points 10
    #    $HYPHY LIBPATH=$HYPHYLIBPATH busted --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.rapidnj.bestTree --branches Internal --output ${FILE}.${GENE}.BUSTED.json --rates 2 --syn-rates 2 --starting-points 10
    #fi

    #if [ -s ${FILE}.${GENE}.ABSREL.json ] 
    #then
    #    echo "Already has ABSREL results"
    #else
    #    mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH absrel --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.rapidnj.bestTree --branches Internal --output ${FILE}.${GENE}.ABSREL.json
    #fi
    
    #if [ -s ${FILE}.${GENE}.FADE.json ] 
    #then
    #    echo "Already has FADE results"
    #else
    #    echo $HYPHY LIBPATH=$HYPHYLIBPATH scripts/reroot-on-oldest.bf --tree ${FILE}.${GENE}.compressed.fas.rapidnj.bestTree --csv $ATTRIBUTES --output ${FILE}.${GENE}.compressed.fas.rooted
    #    $HYPHY LIBPATH=$HYPHYLIBPATH scripts/reroot-on-oldest.bf --tree ${FILE}.${GENE}.compressed.fas.rapidnj.bestTree --csv $ATTRIBUTES --output ${FILE}.${GENE}.compressed.fas.rooted
    #    echo $HYPHY LIBPATH=$HYPHYLIBPATH conv Universal "Keep Deletions" ${FILE}.${GENE}.compressed.fas  ${FILE}.${GENE}.compressed.fas.prot
    #    $HYPHY LIBPATH=$HYPHYLIBPATH conv Universal "Keep Deletions" ${FILE}.${GENE}.compressed.fas  ${FILE}.${GENE}.compressed.fas.prot
    #    echo $HYPHY LIBPATH=$HYPHYLIBPATH fade --alignment ${FILE}.${GENE}.compressed.fas.prot --tree ${FILE}.${GENE}.compressed.fas.rooted --branches Internal
    #    $HYPHY LIBPATH=$HYPHYLIBPATH fade --alignment ${FILE}.${GENE}.compressed.fas.prot --tree ${FILE}.${GENE}.compressed.fas.rooted --branches Internal
    #fi

    
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
