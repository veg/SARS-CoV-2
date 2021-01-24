#!/bin/bash
#PBS -l walltime=72:0:0:0

################################################################################
# Help                                                                         #
################################################################################
Help()
{
   # Display Help
   echo "QA, align, compress, build trees, perform selection analyses on segments"
   echo
   echo "Syntax: ./submit_jobs [-h] data_directory working_directory gene number_cpu MP_flag"
   echo "options:"
   echo "h     Print this Help."
   echo "positional arguments:"   
   echo "data_directory     path to the directory with the data"
   echo "base_driectory     path to the root directory of the pipeline "
   echo "gene               a SARS-CoV-2 gene/ORF/product (one of S,M,N,E,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF10,ORF1a,ORF1b,leader,nsp2,nsp3,nsp4,3C,nsp6,nsp7,nsp8,nsp9,nsp10,RdRp,helicase,exonuclease,endornase,methyltransferase)"
   echo "number_cpu         allocate this many CPUs to a jobs"
   echo "[opt] MP_flag      if set to MP, run the jobs WITHOUT MPI (otherwise run them with MPI)"
}

################################################################################
# ARGUMENTS                                                                    #
################################################################################


check_directory () {
    if [ ! -d "$1" ] 
    then
        echo "$2 ($1) does not exist"
        exit 1
    fi
}

check_file () {
    if [ ! -f "$1" ] 
    then
        echo "$2 ($1) does not exist"
        exit 1
    fi
}

check_exec () {
    if [ ! -f "$1" ] 
    then
        echo "Executable $2 ($1) does not exist"
        exit 1
    fi
}


while getopts "h" option; do
   case $option in
      h) # display Help
         Help
         exit;;
   esac
done

shift $((OPTIND -1))

DIRECTORY=$1
check_directory $DIRECTORY "Data directory"

FILE=$1/sequences
check_file $FILE "Source 'sequences' file"

WORKING_DIR=$2
check_directory $DIRECTORY "Base directory"

LOG=$1/log.txt
GENE=$3
NP=$4
MP_OR_MPI=${5-"MPI"}

################################################################################
# CONFIGURATION                                                                #
################################################################################


# The following section is needed if you submit this as an MPI job via qsub
# 

#export PATH=/usr/local/bin:$PATH
#source /etc/profile.d/modules.sh
#module load aocc/1.3.0
#module load openmpi/gnu/3.0.2
#export LC_ALL=en_US
#export LC_CTYPE=en_US

HYPHY=/usr/local/bin/hyphy
check_exec $HYPHY "HyPhy binary"
## HYPHY binary

## you may comment out if you don't use MPI version of HyPhy
if [ $MP_OR_MPI == "MPI" ]
then
    HYPHYMPI=/usr/local/bin/HYPHYMPI
    check_exec $HYPHYMPI "HyPhy MPI binary"
    RUN_HYPHY="mpirun -np $4 $HYPHYMPI"
else
    RUN_HYPHY="$HYPHY CPU=$NP"
fi

HYPHYLIBPATH=/usr/local/share/hyphy/   
check_directory $DIRECTORY "HyPhy library directory"

MAFFT=/usr/local/bin/mafft
check_exec $MAFFT "MAFFT binary"

RAPIDNJ=/usr/local/bin/rapidnj
check_exec $RAPIDNJ "RAPIDNJ binary"

SEQMAGICK=/usr/local/bin/seqmagick
check_exec $RAPIDNJ "SEQMAGICK script binary"

### The following need to be installed from github.com/veg/hyphy-analyses
PREMSA=/usr/local/share/hyphy-analyses/codon-msa/pre-msa.bf
check_file $PREMSA "pre-msa.bf file"

POSTMSA=/usr/local/share/hyphy-analyses/codon-msa/post-msa.bf
check_file $PREMSA "post-msa.bf file"

PYTHON=/usr/local/bin/python3
check_exec $PYTHON "PYTHON interpreter"

COMPRESSOR=$WORKING_DIR/scripts/compressor.bf
check_file $PREMSA "compressor.bf file"

COMPRESSOR2=$WORKING_DIR/scripts/compressor-2.bf
check_file $PREMSA "compressor.bf file"

ZERO_LENGTHS_FLAGS='--kill-zero-lengths Constrain ENV="_DO_TREE_REBALANCE_=1"'

################################################################################
# FUNCTIONS                                                                    #
################################################################################



################################################################################
# MAIN                                                                         #
################################################################################


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
    echo "================"
    echo "EXTRACTING $GENE"
    echo "================"
    if [ -s ${FILE}.${GENE}_protein.fas ] 
    then
        echo "Already extracted"
    else
        TMP_FILE=${FILE}.${GENE}.tmp
        cp ${FILE} ${TMP_FILE}
        echo "$RUN_HYPHY LIBPATH=$HYPHYLIBPATH $PREMSA --input $TMP_FILE --reference $REFERENCE_SEQUENCE --trim-from $TRIM_FROM --trim-to $TRIM_TO --E 0.01 --N-fraction $N_FRAC --remove-stop-codons Yes"
        $RUN_HYPHY LIBPATH=$HYPHYLIBPATH $PREMSA --input $TMP_FILE --reference $REFERENCE_SEQUENCE --trim-from $TRIM_FROM --trim-to $TRIM_TO --E 0.01 --N-fraction $N_FRAC --remove-stop-codons Yes || exit 1
        mv ${TMP_FILE}_protein.fas ${FILE}.${GENE}_protein.fas
        mv ${TMP_FILE}_nuc.fas ${FILE}.${GENE}_nuc.fas
    fi

    echo "========================="
    echo "INITIAL COMPRESSION $GENE"
    echo "========================="
    if [ -s ${FILE}.${GENE}_protein.compressed.fas ] && [ -s ${FILE}.${GENE}_nuc.compressed.fas ]
    then
        echo "Already extracted"
    else
        echo "$PYTHON python/get-raw-duplicates.py -i ${FILE}.${GENE}_protein.fas -m ${FILE}.${GENE}_nuc.fas -o ${FILE}.${GENE}_protein.compressed.fas -n ${FILE}.${GENE}_nuc.compressed.fas -d ${FILE}.${GENE}_protein.duplicates.json -g ${FILE}.${GENE}_raw_nucleotide.duplicates.json"
        $PYTHON python/get-raw-duplicates.py -i ${FILE}.${GENE}_protein.fas -m ${FILE}.${GENE}_nuc.fas -o ${FILE}.${GENE}_protein.compressed.fas -n ${FILE}.${GENE}_nuc.compressed.fas -d ${FILE}.${GENE}_protein.duplicates.json -g ${FILE}.${GENE}_raw_nucleotide.duplicates.json || exit 1
    fi

    echo "====================="
    echo "ALIGNING PROTEIN DATA"
    echo "====================="
    if [ -s ${FILE}.${GENE}.msa ] 
    then
        echo "Already aligned"
    else
        echo "$MAFFT --auto --thread -1 --addfragments ${FILE}.${GENE}_protein.compressed.fas reference_genes/reference.${GENE}_protein.fas >| ${FILE}.${GENE}.tmp.msa"
        $MAFFT --auto --thread -1 --addfragments ${FILE}.${GENE}_protein.compressed.fas reference_genes/reference.${GENE}_protein.fas >| ${FILE}.${GENE}.tmp.msa || exit 1
        # Remove reference from output
        $PYTHON python/remove-seq.py -i ${FILE}.${GENE}.tmp.msa -r reference_genes/reference.${GENE}_protein.fas -o ${FILE}.${GENE}.msa || exit 1
        rm ${FILE}.${GENE}.tmp.msa 
    fi
        
    if [ -s ${FILE}.${GENE}.compressed.fas ] 
    then
        echo "Already reverse translated"
    else
        echo "$HYPHY LIBPATH=$HYPHYLIBPATH $POSTMSA --protein-msa ${FILE}.${GENE}.msa --nucleotide-sequences ${FILE}.${GENE}_nuc.fas --output ${FILE}.${GENE}.compressed.fas --duplicates ${FILE}.${GENE}_nucleotide.duplicates.json"
        $HYPHY LIBPATH=$HYPHYLIBPATH $POSTMSA --protein-msa ${FILE}.${GENE}.msa --nucleotide-sequences ${FILE}.${GENE}_nuc.fas --output ${FILE}.${GENE}.compressed.fas --duplicates ${FILE}.${GENE}_nucleotide.duplicates.json || exit 1
        #Replace all unknown characters with N
        sed -i '' '/^>/! s/[^ACTG-]/N/g' ${FILE}.${GENE}.compressed.fas || exit 1
    fi

    echo "=========================================="
    echo "MERGING DUPLICATES FROM POST-MSA and PRIOR"
    echo "=========================================="
    if [ -s ${FILE}.${GENE}.duplicates.json ] 
    then
        echo "Already MERGED"
    else
        echo "$PYTHON python/merge-duplicates.py -p ${FILE}.${GENE}_raw_nucleotide.duplicates.json -n ${FILE}.${GENE}_nucleotide.duplicates.json -o ${FILE}.${GENE}.duplicates.json"
        $PYTHON python/merge-duplicates.py -p ${FILE}.${GENE}_raw_nucleotide.duplicates.json -n ${FILE}.${GENE}_nucleotide.duplicates.json -o ${FILE}.${GENE}.duplicates.json || exit 1
        # Fix duplicates 
        $PYTHON python/fix-duplicates.py -d ${FILE}.${GENE}.duplicates.json -m ${FILE}.${GENE}.map.json -o || exit 1
        # Fix header files
        $PYTHON python/update-fasta-duplicates.py -f ${FILE}.${GENE}.compressed.fas -m ${FILE}.${GENE}.map.json || exit 1
    fi

    echo "==================="
    echo "FILTERING ALIGNMENT"
    echo "==================="
    if [ -s ${FILE}.${GENE}.compressed.filtered.fas ] 
    then
        echo "ALIGNMENT ALREADY FILTERED"
    else

        echo "$HYPHY LIBPATH=$HYPHYLIBPATH $COMPRESSOR --msa ${FILE}.${GENE}.compressed.fas --duplicates ${FILE}.${GENE}.duplicates.json --output ${FILE}.${GENE}.variants.csv --json ${FILE}.${GENE}.variants.json --duplicate-out ${FILE}.${GENE}.duplicates-merged.json --regexp '(epi_isl_[0-9A-Za-z]+).+'"
        $HYPHY LIBPATH=$HYPHYLIBPATH $COMPRESSOR --msa ${FILE}.${GENE}.compressed.fas --duplicates ${FILE}.${GENE}.duplicates.json --output ${FILE}.${GENE}.variants.csv --json ${FILE}.${GENE}.variants.json --duplicate-out ${FILE}.${GENE}.duplicates-merged.json --regexp '(epi_isl_[0-9A-Za-z]+).+'

        echo "$HYPHY LIBPATH=$HYPHYLIBPATH $COMPRESSOR2 --msa ${FILE}.${GENE}.compressed.fas --duplicates ${FILE}.${GENE}.duplicates-merged.json --csv ${FILE}.${GENE}.variants.csv --byseq ${FILE}.${GENE}.variants.json --p 0.9 --output ${FILE}.${GENE}.compressed.filtered.fas --output-edits ${FILE}.${GENE}.edits.json --json ${FILE}.${GENE}.filtered.json"
        $HYPHY LIBPATH=$HYPHYLIBPATH $COMPRESSOR2 --msa ${FILE}.${GENE}.compressed.fas --duplicates ${FILE}.${GENE}.duplicates-merged.json --csv ${FILE}.${GENE}.variants.csv --byseq ${FILE}.${GENE}.variants.json --p 0.9 --output ${FILE}.${GENE}.compressed.filtered.fas --json ${FILE}.${GENE}.filtered.json --output-edits ${FILE}.${GENE}.edits.json 
    fi

 
    echo "================"
    echo "MAKING A NJ TREE "
    echo "================"
    if [ -s ${FILE}.${GENE}.compressed.filtered.fas.rapidnj.bestTree ] 
    then
        echo "Already has tree"
    else
        echo "seqmagick convert ${FILE}.${GENE}.compressed.filtered.fas ${FILE}.${GENE}.compressed.filtered.sto"
        echo "rapidnj ${FILE}.${GENE}.compressed.filtered.sto -i sth > ${FILE}.${GENE}.compressed.filtered.fas.rapidnj.bestTree"
        seqmagick convert ${FILE}.${GENE}.compressed.filtered.fas ${FILE}.${GENE}.compressed.filtered.sto || exit 1
        $RAPIDNJ ${FILE}.${GENE}.compressed.filtered.sto -i sth > ${FILE}.${GENE}.compressed.filtered.fas.rapidnj.bestTree || exit 1
        sed -i '' "s/'//g" ${FILE}.${GENE}.compressed.filtered.fas.rapidnj.bestTree || exit 1
        #$RAXML --tree pars{5} --msa ${FILE}.${GENE}.compressed.filtered.fas --threads 4 --model GTR+G --force
    fi


    echo "========================="
    echo "RUNNING THE SLAC ANALYSIS"
    echo "========================="
     if [ -s ${FILE}.${GENE}.SLAC.json.gz ] 
    then
        echo "Already has SLAC results"
    else
        echo $RUN_HYPHY LIBPATH=$HYPHYLIBPATH slac $ZERO_LENGTHS_FLAGS --alignment ${FILE}.${GENE}.compressed.filtered.fas --tree ${FILE}.${GENE}.compressed.filtered.fas.rapidnj.bestTree --branches All --samples 0 --output ${FILE}.${GENE}.SLAC.json --intermediate-fits ${FILE}.${GENE}.cache
        $RUN_HYPHY LIBPATH=$HYPHYLIBPATH slac $ZERO_LENGTHS_FLAGS --alignment ${FILE}.${GENE}.compressed.filtered.fas --tree ${FILE}.${GENE}.compressed.filtered.fas.rapidnj.bestTree --branches All --samples 0 --output ${FILE}.${GENE}.SLAC.json --intermediate-fits ${FILE}.${GENE}.cache  || exit 1
        gzip ${FILE}.${GENE}.SLAC.json 
    fi

    echo "========================="
    echo "RUNNING THE FEL  ANALYSIS"
    echo "========================="
    if [ -s ${FILE}.${GENE}.FEL.json.gz ] 
    then
        echo "Already has FEL results"
    else
        echo $RUN_HYPHY LIBPATH=$HYPHYLIBPATH fel $ZERO_LENGTHS_FLAGS --alignment ${FILE}.${GENE}.compressed.filtered.fas --tree ${FILE}.${GENE}.compressed.filtered.fas.rapidnj.bestTree --branches Internal --output ${FILE}.${GENE}.FEL.json --intermediate-fits ${FILE}.${GENE}.cache
        $RUN_HYPHY LIBPATH=$HYPHYLIBPATH fel $ZERO_LENGTHS_FLAGS --alignment ${FILE}.${GENE}.compressed.filtered.fas --tree ${FILE}.${GENE}.compressed.filtered.fas.rapidnj.bestTree --branches Internal --output ${FILE}.${GENE}.FEL.json --intermediate-fits ${FILE}.${GENE}.cache || exit 1
        gzip ${FILE}.${GENE}.FEL.json 
    fi

    echo "=========================="
    echo "RUNNING THE MEME  ANALYSIS"
    echo "=========================="
    if [ -s ${FILE}.${GENE}.MEME.json.gz ] 
    then
        echo "Already has MEME results"
    else
        echo $RUN_HYPHY LIBPATH=$HYPHYLIBPATH meme $ZERO_LENGTHS_FLAGS --alignment ${FILE}.${GENE}.compressed.filtered.fas --tree ${FILE}.${GENE}.compressed.filtered.fas.rapidnj.bestTree --branches Internal --output ${FILE}.${GENE}.MEME.json --intermediate-fits ${FILE}.${GENE}.cache
        $RUN_HYPHY LIBPATH=$HYPHYLIBPATH meme $ZERO_LENGTHS_FLAGS --alignment ${FILE}.${GENE}.compressed.filtered.fas --tree ${FILE}.${GENE}.compressed.filtered.fas.rapidnj.bestTree --branches Internal --output ${FILE}.${GENE}.MEME.json --intermediate-fits ${FILE}.${GENE}.cache || exit 1
        gzip ${FILE}.${GENE}.MEME.json 
    fi
    
    if [ -s {FILE}.${GENE}.duplicates.json.gz ]
    then
        echo ;
    else
        gzip {FILE}.${GENE}.duplicates.json
    fi

    TMP_FILE=${FILE}.${GENE}.tmp
    rm -f $TMP_FILE
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
