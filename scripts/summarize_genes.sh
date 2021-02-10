#!/bin/bash

P3=python3
JQ=/usr/local/bin/jq


################################################################################
# Help                                                                         #
################################################################################
Help()
{
   # Display Help
   echo "Process results of the SARS-CoV-2 pipeline and generate summary JSON and TSV files"
   echo
   echo "Syntax: ./submit_jobs [-hvs] [-d dir] [-b base_dir] analysis_file"
   echo "options:"
   echo "h     Print this Help."
   echo "v     Generate a variants file."
   echo "b     path to the root directory of the pipeline (default = `pwd`)"
   echo "d     path to the data directory where meta information is located (default same as where the analysis file is)"
    echo "positional arguments:"
   echo "[req] analysis_file path to 'sequences' file"
}

check_file () {
    if [ ! -f "$1" ] 
    then
        echo "$2 ($1) does not exist"
        exit 1
    fi
}

BASE_DIR=`pwd`
DO_VARIANTS="0"

while getopts "hvb:d:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      v)
        DO_VARIANTS="1"
        ;;
      b)
        BASE_DIR=$OPTARG
        ;;
      d)
        DATA_DIR=$OPTARG
        ;;
   esac
done

shift $((OPTIND -1))

if [ $# -lt 1 ]
  then
    echo "Analysis path is required"
    exit 1
fi


FILE=$1
check_file $FILE "Analysis sequence file"

echo "...WILL RUN ANALYSES ON $FILE"
if [ -z $DATA_DIR ]; 
then
    DATA_DIR=$(dirname  $FILE)
fi

echo "...WILL USE META INFORMATION FROM $DATA_DIR"
echo "...COMPRESSING SEQUENCE ANNOTATIONS"

$P3 ${BASE_DIR}/python/obfuscate-master.py -i ${DATA_DIR}/annotation.json -o ${DATA_DIR}/sequence-info.json -m  ${DATA_DIR}/map.json



if [ $DO_VARIANTS == "1" ]; 
then
	$P3 ${BASE_DIR}/python/stitch_fasta.py -d $FILE -o ${FILE}.variants.json -r ${BASE_DIR}/reference_genes/sc2.mmi -m ${DATA_DIR}/map.json -c 0.00001
fi



genes=(leader nsp2 nsp3 nsp4 3C nsp6 nsp7 nsp8 nsp9 nsp10 helicase exonuclease endornase  S E M N ORF3a ORF6 ORF7a ORF8 RdRp methyltransferase)
offsets=(265 805 2719 8554 10054 10972 11842 12091 12685 13024 16236 18039 19620 21562 26244 26522 28273 25392 27201 27393 27893 13440 20658)
fragments=(ORF1a ORF1a ORF1a ORF1a ORF1a ORF1a ORF1a ORF1a ORF1a ORF1a ORF1b ORF1b ORF1b S E M N ORF3a ORF6 ORF7a ORF8 ORF1b ORF1b)
shifts=(0        180   818   2763  3263  3569  3859  3942  4140  4253  922   1523  2050  0 0 0 0 0     0    0     0    -10   2396)
add_one=(0       0      0     0     0     0     0     0     0     0    1     1     1     0 0 0 0 0     0    0     0    1    1)


#genes=(3C)
#offsets=(1004) 
#fragments=(ORF1a)
#shifts=(3263)
#add_one=(0)

#genes=(S)
#offsets=(21562) 
#fragments=(S)
#shifts=(0)
#add_one=(0)


OMNIBUS_FILE=${FILE}.report.json
FINAL_RESULT=$(dirname "${OMNIBUS_FILE})")"/report.json"
ANNOTATION=$(dirname "${OMNIBUS_FILE})")"/comparative-annotation.json"
echo '{}' > ${OMNIBUS_FILE}

cp data/comparative-annotation-between.json ${ANNOTATION}
rm ${FINAL_RESULT}

for i in ${!genes[@]}; do
    GENE=${genes[i]}
    OFFSET=${offsets[i]}
    FRAGMENT=${fragments[i]}
    SHIFT=${shifts[i]}
    ADDSHIFT=${add_one[i]}

    echo ">>>>>>>>>>"
    echo $GENE
    if [ -s ${FILE}.${GENE}.json ]; then
        echo  "Already done"
    else
        if $P3 ${BASE_DIR}/python/summarize-gene.py --map ${DATA_DIR}/map.json -T ${BASE_DIR}/data/ctl/epitopes.json -D ${DATA_DIR}/annotation.json --frame_shift ${ADDSHIFT} -d ${FILE}.${GENE}.duplicates.json.gz -s ${FILE}.${GENE}.SLAC.json.gz -f ${FILE}.${GENE}.FEL.json.gz -m ${FILE}.${GENE}.MEME.json.gz -P 0.05 --output  ${FILE}.${GENE}.json -c ${FILE}.${GENE}.compressed.fas -E ${BASE_DIR}/data/evo_annotation.json -F $FRAGMENT --fragment_shift $SHIFT -S $OFFSET -O $ANNOTATION > ${FILE}.${GENE}.json; then
            echo "Complete"
            echo "<<<<<<<<<<"
            $JQ -c -s ".[0] + {\"$GENE\" : .[1]}" $OMNIBUS_FILE ${FILE}.${GENE}.json > ${OMNIBUS_FILE}.2
            mv ${OMNIBUS_FILE}.2 ${OMNIBUS_FILE}
        else
            echo "Error"
            echo "<<<<<<<<<<"
            #exit 1
        fi
    fi
    
done;

mv ${OMNIBUS_FILE} ${FINAL_RESULT}

TRAJECTORY=$(dirname "${OMNIBUS_FILE})")
TRAJECTORY=$(dirname "${TRAJECTORY})")
$P3 ${BASE_DIR}/python/temporal-summary-paper.py -d $TRAJECTORY > ${TRAJECTORY}/temporal-gene-properties.csv
$P3 ${BASE_DIR}/python/temporal-summary-paper-sites.py -d $TRAJECTORY > ${TRAJECTORY}/temporal-gene-sites.csv
  

#if (( ${2:-0} == 2 || ${2:-0} == 3)); then
#	$P3 python/export-sites-to-tsv.py -f data/comparative-annotation.json > data/comparative-annotation.tsv
	#bzip2 --best ${OMNIBUS_FILE}
#fi

