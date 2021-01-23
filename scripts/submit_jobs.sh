#!/bin/bash

################################################################################
# Help                                                                         #
################################################################################
Help()
{
   # Display Help
   echo "Top level job dispatcher for the SARS-CoV-2 selection analysis pipeline"
   echo
   echo "Syntax: ./submit_jobs [-hm] [-b dir] [-l dir] [-s number] [-L number] [-q queue] analysis_directory "
   echo "options:"
   echo "h     Print this Help."
   echo "m     Enable MPI mode (requires qsub; otherwise run serially)."
   echo "b     path to the root directory of the pipeline (default = `pwd`)"
   echo "l     path to the directory for MPI job logs (default `pwd`/logs/analysis_dir_name)"
   echo "s     allocate this many CPUs to 'small' jobs (default 16)"
   echo "L     allocate this many CPUs to 'large' jobs (default 48)"
   echo "q     MPI-queue send jobs to this MPI queue (default : no queue specification)"
   echo "positional arguments:"
   echo "[req] analysis_directory path to the directory with the sequence file"
}

################################################################################
################################################################################
# Main program                                                                 #
################################################################################
################################################################################


DO_MPI="0"
BASE_DIR=`pwd`

while getopts "hmb:l:s:L:q:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      m)
        DO_MPI="1"
        ;;
      b)
        BASE_DIR=$OPTARG
        ;;
      l)
        LOG_DIR=$OPTARG
        ;;
      s)
        SMALLPPN=$OPTARG
        ;;
      L)
        LARGEPPN=$OPTARG
        ;;
      q)
        QUEUE=$OPTARG
        ;;
   esac
done

shift $((OPTIND -1))

if [ $# -lt 1 ]
  then
    echo "Analysis path is required"
    exit 1
fi

DATA_DIR=${1}
JOB_NAME=$(basename $DATA_DIR)

echo "...WILL RUN ANALYSES ON $DATA_DIR"
echo "...BASE DIRECTORY $BASE_DIR"

if [ -z $LOG_DIR ] 
then
    LOG_DIR=`pwd`/logs/$JOB_NAME/
fi
echo "...LOG DIRECTORY $LOG_DIR"

SMALLPPN=${SMALLPPN:-"16"}
echo "...$SMALLPPN processors per small job"
LARGEPPN=${LARGEPPN:-"48"}
echo "...$LARGEPPN processors per large job"
QUEUE=${QUEUE:-""}

if [ $DO_MPI == 1 ]
then
    echo "...$QUEUE for MPI submissions"
fi

submit_a_job () {
   if [ $DO_MPI = "0" ] 
   then
         bash $2/scripts/extract_genes.sh $1 $2 $3 $4 "MP" 
   else
        if [ -z $6 ]
        then
            qsub -l nodes=1:ppn=$4 -d `pwd` -o $5-e $5 -F "$1 $2 $3 $4" $2/scripts/extract_genes.sh        
        else
            qsub -l nodes=1:ppn=$4 -d `pwd` -o $5-e $5 -q $6 -F "$1 $2 $3 $4" $2/scripts/extract_genes.sh
        fi
   fi
} 


large=(S nsp3)

for i in ${!large[@]}; do
    GENE=${large[i]}
    submit_a_job $DATA_DIR $BASE_DIR $GENE $LARGEPPN $LOG_DIR $QUEUE 
done

genes=(M N E ORF3a ORF6 ORF7a ORF7b ORF8)

for i in ${!genes[@]}; do
    GENE=${genes[i]}
    submit_a_job $DATA_DIR $BASE_DIR $GENE $SMALLPPN $LOG_DIR $QUEUE 
done

products=(leader nsp2 nsp4 3C nsp6 nsp7 nsp8 nsp9 nsp10 RdRp helicase exonuclease endornase methyltransferase ORF10)

for i in ${!products[@]}; do
    GENE=${products[i]}
    submit_a_job $DATA_DIR $BASE_DIR $GENE $LARGEPPN $LOG_DIR $QUEUE 
done

exit 0
