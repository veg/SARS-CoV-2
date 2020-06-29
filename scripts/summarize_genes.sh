FILE=$1

#python3 python/stitch_fasta.py -d $FILE -o ${FILE}.combined.fas
#python3 python/extract_variants.py -i ${FILE}.combined.fas  -o ${FILE}.variants.json -c 5

#rm data/mafs.csv
#rm data/evo_freqs.csv

#gene_coordinates = [[265,13482, 'ORF1a', 0],
#  [13467,21554, 'ORF1b',-1],
#  [21562,25383,'S',0],
#  [25392,26219,'ORF3a',0],
#  [26244,26471,'E',0],
#  [26522,27190,'M',0],
#  [27201,27386,'ORF6',0],
#  [27393,27758,'ORF7a',0],
#  [27755,27886,'ORF7b',0],
#  [28273,29532,'N',0],
#  [29557,29673,'ORF10',0],
#  [265,804,'leader',0],
#  [805,2718, 'nsp2',0],
#  [2719,8553, 'nsp3',0],
#  [8554,10053, 'nsp4',0],
#  [10054,10971, '3C',0],
#  [10972,11841, 'nsp6',0],
#  [11842,12090, 'nsp7',0],
#  [12091,12684, 'nsp8',0],
#  [12685,13023, 'nsp9',0],
#  [13024,13441, 'nsp10',0],
#  [16236,18038, 'helicase',0],
#  [18039,19619, 'exonuclease',0],
#  [19620,20657, 'endornase',0],
#  [20658,21551, 'methyltransferase',0],
#  [26244,26471, 'E',0],
#  [29557,29673, 'ORF10',0],
#]

genes=(leader nsp2 nsp3 nsp4 3C nsp6 nsp7 nsp8 nsp9 nsp10 helicase exonuclease endornase  S E M N ORF3a ORF6 ORF7a ORF8 RdRp methyltransferase)
offsets=(265 805 2719 8554 10054 10972 11842 12091 12685 13024 16236 18039 19620 21562 26244 26522 28273 25392 27201 27393 27893 13440 20658)
fragments=(ORF1a ORF1a ORF1a ORF1a ORF1a ORF1a ORF1a ORF1a ORF1a ORF1a ORF1b ORF1b ORF1b S E M N ORF3a ORF6 ORF7a ORF8 ORF1b ORF1b)
shifts=(0        180   818   2763  3263  3569  3859  3942  4140  4253  922   1523  2050  0 0 0 0 0     0    0     0    -10   2396)
add_one=(0       0      0     0     0     0     0     0     0     0    1     1     1     0 0 0 0 0     0    0     0    1    1)

#genes=(leader)
#offsets=(265) 
#fragments=(ORF1a)
#shifts=(0)
#add_one=(0)

#for GENE in {S,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF8}; do

ANNOTATION=${FILE}.annotation.json

#cp data/comparative-annotation-between.json ${FILE}.annotation.json

for i in ${!genes[@]}; do
    GENE=${genes[i]}
    OFFSET=${offsets[i]}
    FRAGMENT=${fragments[i]}
    SHIFT=${shifts[i]}
    ADDSHIFT=${add_one[i]}
    
    echo $GENE
    #if [ -s ${FILE}.${GENE}.withref.fas ]
    #then 
    #    echo "Already has alignment with reference"
    #else
    #    echo ${FILE}.${GENE}.compressed.fas
    #    $MAFFT --add data/reference_genes/${GENE}.fas --reorder ${FILE}.${GENE}.compressed.fas > ${FILE}.${GENE}.withref.fas
    #    cp ${FILE}.${GENE}.withref.fas ${FILE}.${GENE}.bkup.withref.fas
    #fi 
    python3 python/summarize-gene.py -T data/ctl/epitopes.json -D data/db/master-no-fasta.json --frame_shift ${ADDSHIFT} -d ${FILE}.${GENE}.duplicates.json -u ${FILE}.${GENE}.compressed.fas.FUBAR.json -s ${FILE}.${GENE}.SLAC.json -f ${FILE}.${GENE}.FEL.json -m ${FILE}.${GENE}.MEME.json -P 0.1 --output  ${FILE}.${GENE}.json -c ${FILE}.${GENE}.compressed.fas -E data/evo_annotation.json -F $FRAGMENT -A data/mafs.csv -V data/evo_freqs.csv --fragment_shift $SHIFT -S $OFFSET -O $ANNOTATION > ${FILE}.${GENE}.json
done;

#cp $ANNOTATION data/comparative-annotation.json
#python3 python/export-sites-to-tsv.py -f data/comparative-annotation.json > data/comparative-annotation.tsv