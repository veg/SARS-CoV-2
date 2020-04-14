#!/bin/bash
#PBS -l nodes=1:ppn=16
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
NP=16
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


function run_a_gene {


GENE=$1
REFERENCE_SEQUENCE=$2
TRIM_FROM=$3
TRIM_TO=$4
N_FRAC=$5

if [ -s ${FILE}.${GENE}.1.json ]
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
        echo "mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH $PREMSA --input $TMP_FILE --reference $REFERENCE_SEQUENCE --trim-from $TRIM_FROM --trim-to $TRIM_TO --E 0.01 --N-fraction $N_FRAC"
        mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH $PREMSA --input $TMP_FILE --reference $REFERENCE_SEQUENCE --trim-from $TRIM_FROM --trim-to $TRIM_TO --E 0.01 --N-fraction $N_FRAC
        mv ${TMP_FILE}_protein.fas ${FILE}.${GENE}_protein.fas
        mv ${TMP_FILE}_nuc.fas ${FILE}.${GENE}_nuc.fas
    fi
    
    echo "ALIGNING PROTEIN DATA"
    if [ -s ${FILE}.${GENE}.msa ] 
    then
        echo "Already aligned"
    else
        echo $MAFFT ${FILE}.${GENE}_protein.fas > ${FILE}.${GENE}.msa
        $MAFFT ${FILE}.${GENE}_protein.fas > ${FILE}.${GENE}.msa 2> mafft.error.log
    fi
        
    if [ -s ${FILE}.${GENE}.compressed.fas ] 
    then
        echo "Already reverse translated"
    else
        $HYPHY LIBPATH=$HYPHYLIBPATH $POSTMSA --protein-msa ${FILE}.${GENE}.msa --nucleotide-sequences ${FILE}.${GENE}_nuc.fas --output ${FILE}.${GENE}.compressed.fas --duplicates ${FILE}.${GENE}.duplicates.json
        $HYPHY LIBPATH=$HYPHYLIBPATH $POSTMSA --protein-msa ${FILE}.${GENE}.msa --nucleotide-sequences ${FILE}.${GENE}_nuc.fas --compress No --output ${FILE}.${GENE}.all.fas    
    fi
    
    if [ -s ${FILE}.${GENE}.withref.fas ]
    then 
        echo "Already has alignment with reference"
    else
        echo $MAFFT --add $REFERENCE_SEQUENCE --reorder ${FILE}.${GENE}.all.fas > ${FILE}.${GENE}.withref.fas
        $MAFFT --add $REFERENCE_SEQUENCE --reorder ${FILE}.${GENE}.all.fas > ${FILE}.${GENE}.withref.fas 2> mafft.error.log
    fi 

    if [ -s ${FILE}.${GENE}.tn93 ] 
    then
        echo "Already computed TN93"
    else
        $TN93 -q -t 0.05 ${FILE}.${GENE}.withref.fas > ${FILE}.${GENE}.tn93 2> ${FILE}.${GENE}.tn93.json
        echo python3 $WORKING_DIR/python/tabulate-diversity-divergence.py -j $MASTER -t ${FILE}.${GENE}.tn93 > $DIRECTORY/evolution.${GENE}.csv
        python3 $WORKING_DIR/python/tabulate-diversity-divergence.py -j $MASTER -t ${FILE}.${GENE}.tn93 > $DIRECTORY/evolution.${GENE}.csv
    fi

    if [ -s ${FILE}.${GENE}.compressed.fas.raxml.bestTree ] 
    then
        echo "Already has tree"
    else
        $RAXML --msa ${FILE}.${GENE}.compressed.fas --threads 4 --model GTR+G --force
    fi
    
    if [ -s $WORKING_DIR/${FILE}.${GENE}.SLAC.json ] 
    then
        echo "Already has SLAC results"
    else
        echo mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH slac --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.SLAC.json
        mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH slac --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.SLAC.json
    fi

    if [ -s ${FILE}.${GENE}.FEL.json ] 
    then
        echo "Already has FEL results"
    else
        echo mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH fel --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.FEL.json
        mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH fel --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.FEL.json
    fi

    if [ -s ${FILE}.${GENE}.MEME.json ] 
    then
        echo "Already has MEME results"
    else
        echo mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH meme --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.MEME.json
        mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH meme --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.MEME.json
    fi

    if [ -s ${FILE}.${GENE}.PRIME.json ] 
    then
        echo "Already has PRIME results"
    else
        echo mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH prime --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.PRIME.json
        mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH prime --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.PRIME.json
    fi

    echo python3 $WORKING_DIR/python/summarize-gene.py -D $MASTERNOFASTA -d ${FILE}.${GENE}.duplicates.json -s ${FILE}.${GENE}.SLAC.json -f ${FILE}.${GENE}.FEL.json -m ${FILE}.${GENE}.MEME.json -P 0.1 -p ${FILE}.${GENE}.PRIME.json --output  ${FILE}.${GENE}.json -E ${FILE}.${GENE}.evo_annotation.json -c ${FILE}.${GENE}.withref.fas
    python3 $WORKING_DIR/python/summarize-gene.py -D $MASTERNOFASTA -d ${FILE}.${GENE}.duplicates.json -s ${FILE}.${GENE}.SLAC.json -f ${FILE}.${GENE}.FEL.json -m ${FILE}.${GENE}.MEME.json -P 0.1 -p ${FILE}.${GENE}.PRIME.json --output  ${FILE}.${GENE}.json -E ${FILE}.${GENE}.evo_annotation.json -c ${FILE}.${GENE}.withref.fas


    #if [ -s ${FILE}.${GENE}.BGM.json ] 
    #then
    #    echo "Already has BGM results"
    #else
    #    echo mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH bgm --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --branches Internal --min-subs 2 --type codon
    #    mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH bgm --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --branches Internal --min-subs 2 --type codon
    #fi

    #if [ -s ${FILE}.${GENE}.BUSTED.json ] 
    #then
    #    echo "Already has BUSTED results"
    #else
    #    echo $HYPHY LIBPATH=$HYPHYLIBPATH busted --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.BUSTED.json --rates 2 --syn-rates 2 --starting-points 10
    #    $HYPHY LIBPATH=$HYPHYLIBPATH busted --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.BUSTED.json --rates 2 --syn-rates 2 --starting-points 10
    #fi

    #if [ -s ${FILE}.${GENE}.ABSREL.json ] 
    #then
    #    echo "Already has ABSREL results"
    #else
    #    mpirun -np $NP $HYPHYMPI LIBPATH=$HYPHYLIBPATH absrel --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.ABSREL.json
    #fi
    
    #if [ -s ${FILE}.${GENE}.FADE.json ] 
    #then
    #    echo "Already has FADE results"
    #else
    #    echo $HYPHY LIBPATH=$HYPHYLIBPATH scripts/reroot-on-oldest.bf --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --csv $ATTRIBUTES --output ${FILE}.${GENE}.compressed.fas.rooted
    #    $HYPHY LIBPATH=$HYPHYLIBPATH scripts/reroot-on-oldest.bf --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --csv $ATTRIBUTES --output ${FILE}.${GENE}.compressed.fas.rooted
    #    echo $HYPHY LIBPATH=$HYPHYLIBPATH conv Universal "Keep Deletions" ${FILE}.${GENE}.compressed.fas  ${FILE}.${GENE}.compressed.fas.prot
    #    $HYPHY LIBPATH=$HYPHYLIBPATH conv Universal "Keep Deletions" ${FILE}.${GENE}.compressed.fas  ${FILE}.${GENE}.compressed.fas.prot
    #    echo $HYPHY LIBPATH=$HYPHYLIBPATH fade --alignment ${FILE}.${GENE}.compressed.fas.prot --tree ${FILE}.${GENE}.compressed.fas.rooted --branches Internal
    #    $HYPHY LIBPATH=$HYPHYLIBPATH fade --alignment ${FILE}.${GENE}.compressed.fas.prot --tree ${FILE}.${GENE}.compressed.fas.rooted --branches Internal
    #fi

    
fi

}

case $GENE in

  S)
    echo -n "Analyzing S Gene"
    run_a_gene "S" "$WORKING_DIR/data/reference_genes/S.fas" "20000" "27000" 0.005
    ;;

  M)
    echo -n "Analyzing M gene"
    run_a_gene "M" "$WORKING_DIR/data/reference_genes/M.fas" "25000" "30000" 0.01
    ;;

  N)
    echo -n "Analyzing N gene"
    run_a_gene "N" "$WORKING_DIR/data/reference_genes/N.fas" "26000" "35000" 0.01
    ;;

  ORF3a)
    echo -n "Analyzing ORF3a gene"
    run_a_gene "ORF3a" "$WORKING_DIR/data/reference_genes/ORF3a.fas" "24000" "27000" 0.01
    ;;

  ORF7a)
    echo -n "Analyzing ORF7a gene"
    run_a_gene "ORF7a" "$WORKING_DIR/data/reference_genes/ORF7a.fas" "26000" "35000" 0.01
    ;;

  ORF8)
    echo -n "Analyzing ORF8 gene"
    run_a_gene "ORF8" "$WORKING_DIR/data/reference_genes/ORF8.fas" "26000" "35000" 0.01
    ;;

  ORF1a)
    echo -n "Analyzing ORF1a gene"
    run_a_gene "ORF1a" "$WORKING_DIR/data/reference_genes/ORF1a.fas" "1" "15000" 0.001
    ;;

  ORF1b)
    echo -n "Analyzing ORF1b gene"
    run_a_gene "ORF1b" "$WORKING_DIR/data/reference_genes/ORF1b.fas" "12000" "24000" 0.001
    ;;

  *)
    echo -n "Unknown Gene, exiting"
    ;;
esac
