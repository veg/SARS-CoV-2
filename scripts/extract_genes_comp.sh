FILE=$1
NP=10
HYPHY=/usr/local/bin/hyphy
HYPHYMPI=/usr/local/bin/HYPHYMPI
MAFFT=/usr/local/bin/mafft
RAXML=/usr/local/bin/raxml-ng
TN93=/usr/local/bin/tn93

function RunAGene {



GENE=$1
REFERENCE_SEQUENCE=$2
TRIM_FROM=$3
TRIM_TO=$4
N_FRAC=$5
N_FRAC_STRICT=$6

if [ -s ${FILE}.${GENE}.1.json ]
then 
   echo "$GENE alignment already processed"
else

    echo "EXTRACTING $GENE"
    if [ -s ${FILE}.${GENE}_protein.fas ] 
    then
        echo "Already extracted"
    else
        #echo         $HYPHY /Users/sergei/Development/hyphy-analyses/codon-msa/pre-msa.bf --input $FILE --reference $REFERENCE_SEQUENCE --E 0.000001 --trim-from $TRIM_FROM --trim-to $TRIM_TO --N-fraction $N_FRAC --keep-reference Yes
        mpirun -np $NP --use-hwthread-cpus $HYPHYMPI /Users/sergei/Development/hyphy-analyses/codon-msa/pre-msa.bf --input $FILE --reference $REFERENCE_SEQUENCE --E 0.1 --trim-from $TRIM_FROM --trim-to $TRIM_TO --N-fraction $N_FRAC --keep-reference Yes
        mv ${FILE}_protein.fas ${FILE}.${GENE}_protein.fas
        mv ${FILE}_nuc.fas ${FILE}.${GENE}_nuc.fas
        #$HYPHY scripts/filter_on_N.bf --nuc ${FILE}.${GENE}_nuc.fas --prot ${FILE}.${GENE}_protein.fas --N $N_FRAC_STRICT --output-nuc ${FILE}.${GENE}_nuc.strict.fas --output-prot ${FILE}.${GENE}_protein.strict.fas
    fi
    
    
    if [ -s ${FILE}.${GENE}.msa ] 
    then
        echo "Already aligned"
    else
        $MAFFT --auto ${FILE}.${GENE}_protein.fas > ${FILE}.${GENE}.msa
        cp ${FILE}.${GENE}.msa ${FILE}.${GENE}.bkup.msa
    fi

    #if [ -s ${FILE}.${GENE}.strict.msa ] 
    #then
    #    echo "Already aligned strict"
    #else
    #    $MAFFT --auto ${FILE}.${GENE}_protein.strict.fas > ${FILE}.${GENE}.strict.msa
    #fi
        
    if [ -s ${FILE}.${GENE}.duplicates.json ] 
    then
        echo "Already reverse translated"
    else
        $HYPHY /Users/sergei/Development/hyphy-analyses/codon-msa/post-msa.bf --protein-msa ${FILE}.${GENE}.msa --nucleotide-sequences ${FILE}.${GENE}_nuc.fas --output ${FILE}.${GENE}.compressed.fas --duplicates ${FILE}.${GENE}.duplicates.json
        $HYPHY /Users/sergei/Development/hyphy-analyses/codon-msa/post-msa.bf --protein-msa ${FILE}.${GENE}.msa --nucleotide-sequences ${FILE}.${GENE}_nuc.fas --compress No --output ${FILE}.${GENE}.all.fas    
    fi

    #if [ -s ${FILE}.${GENE}.duplicates.strict.json ] 
    #then
    #    echo "Already reverse translated strict"
    #else
    #    $HYPHY /Users/sergei/Development/hyphy-analyses/codon-msa/post-msa.bf --protein-msa ${FILE}.${GENE}.strict.msa --nucleotide-sequences ${FILE}.${GENE}_nuc.strict.fas --output ${FILE}.${GENE}.compressed.strict.fas --duplicates ${FILE}.${GENE}.duplicates.strict.json
    #fi
    
     
    #if [ -s ${FILE}.${GENE}.withref.fas ]
    #then 
    #    echo "Already has alignment with reference"
    #else
    #    $MAFFT --add $REFERENCE_SEQUENCE --reorder ${FILE}.${GENE}.all.fas > ${FILE}.${GENE}.withref.fas
    #    cp ${FILE}.${GENE}.withref.fas ${FILE}.${GENE}.bkup.withref.fas
    #fi 

    #if [ -s ${FILE}.${GENE}.withref.strict.fas ]
    #then 
    #    echo "Already has alignment with reference strict"
    #else
    #    $MAFFT --add $REFERENCE_SEQUENCE --reorder ${FILE}.${GENE}.compressed.strict.fas > ${FILE}.${GENE}.withref.strict.fas
    #fi 
    
 
    
    if [ -s ${FILE}.${GENE}.compressed.fas.raxml.bestTree ] 
    then
        echo "Already has tree"
    else
        $RAXML --tree pars{10} --msa ${FILE}.${GENE}.compressed.fas --model GTR+G --force --redo
    fi
    

    #if [ -s ${FILE}.${GENE}.compressed.strict.fas.raxml.bestTree ] 
    #then
    #    echo "Already has strict tree"
    #else
    #    $RAXML --msa ${FILE}.${GENE}.compressed.strict.fas --model GTR+G --force
    #fi

    if [ -s ${FILE}.${GENE}.SLAC.json ] 
    then
        echo "Already has SLAC results"
    else
        mpirun --use-hwthread-cpus -np $NP $HYPHYMPI slac --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.SLAC.json
    fi

    #if [ -s ${FILE}.${GENE}.SLAC.strict.json ] 
    #then
    #    echo "Already has strict SLAC results"
    #else
    #    mpirun --use-hwthread-cpus -np $NP $HYPHYMPI slac --alignment ${FILE}.${GENE}.compressed.strict.fas --tree ${FILE}.${GENE}.compressed.strict.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.SLAC.strict.json
    #fi

    #if [ -s ${FILE}.${GENE}.FEL.json ] 
    #then
    #    echo "Already has FEL results"
    #else
    #    mpirun --use-hwthread-cpus -np $NP $HYPHYMPI fel --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.FEL.json
    #fi

    #if [ -s ${FILE}.${GENE}.FEL.strict.json ] 
    #then
    #    echo "Already has strict FEL results"
    #else
    #    mpirun --use-hwthread-cpus -np $NP $HYPHYMPI fel --alignment ${FILE}.${GENE}.compressed.strict.fas --tree ${FILE}.${GENE}.compressed.strict.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.FEL.strict.json
    #fi


    #if [ -s ${FILE}.${GENE}.MEME.json ] 
    #then
    #    echo "Already has MEME results"
    #else
    #    mpirun --use-hwthread-cpus -np $NP $HYPHYMPI meme --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.MEME.json
    #fi
    

    #if [ -s ${FILE}.${GENE}.MEME.strict.json ] 
    #then
    #    echo "Already has strict MEME results"
    #else
    #    mpirun --use-hwthread-cpus -np $NP $HYPHYMPI meme --alignment ${FILE}.${GENE}.compressed.strict.fas --tree ${FILE}.${GENE}.compressed.strict.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.MEME.strict.json
    #fi

    #if [ -s ${FILE}.${GENE}.PRIME.json ] 
    #then
    #    echo "Already has PRIME results"
    #else
    #    mpirun --use-hwthread-cpus -np $NP $HYPHYMPI prime --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.PRIME.json
    #fi
    
    #if [ -s ${FILE}.${GENE}.PRIME.strict.json ] 
    #then
    #    echo "Already has strict PRIME results"
    #else
    #    mpirun --use-hwthread-cpus -np $NP $HYPHYMPI prime --alignment ${FILE}.${GENE}.compressed.strict.fas --tree ${FILE}.${GENE}.compressed.strict.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.PRIME.strict.json
    #fi

    #python3 python/summarize-gene.py -D data/db/master-no-fasta.json -d ${FILE}.${GENE}.duplicates.json -s ${FILE}.${GENE}.SLAC.json -f ${FILE}.${GENE}.FEL.json -m ${FILE}.${GENE}.MEME.json -P 0.1 -p ${FILE}.${GENE}.PRIME.json --output  ${FILE}.${GENE}.json -c ${FILE}.${GENE}.withref.fas
    #python3 python/summarize-gene.py -D data/db/master-no-fasta.json -d ${FILE}.${GENE}.duplicates.strict.json -s ${FILE}.${GENE}.SLAC.strict.json -f ${FILE}.${GENE}.FEL.strict.json -m ${FILE}.${GENE}.MEME.strict.json -P 0.1 -p ${FILE}.${GENE}.PRIME.strict.json --output  ${FILE}.${GENE}.json -c ${FILE}.${GENE}.withref.strict.fas

    #if [ -s ${FILE}.${GENE}.BGM.json ] 
    #then
    #    echo "Already has BGM results"
    #else
    #   $HYPHY bgm --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --min-subs 2 --type codon
    #   mv ${FILE}.${GENE}.compressed.fas.BGM.json ${FILE}.${GENE}.BGM.json
    #fi

    #if [ -s ${FILE}.${GENE}.FADE.json ] 
    #then
    #    echo "Already has FADE results"
    #else
    #    $HYPHY scripts/reroot-on-oldest.bf --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --csv data/attributes.csv --output ${FILE}.${GENE}.compressed.fas.rooted
    #    $HYPHY conv Universal "Keep Deletions" ${FILE}.${GENE}.compressed.fas  ${FILE}.${GENE}.compressed.fas.prot
    #    $HYPHY fade --alignment ${FILE}.${GENE}.compressed.fas.prot --tree ${FILE}.${GENE}.compressed.fas.rooted 
    #fi


    #if [ -s ${FILE}.${GENE}.BUSTED.json ] 
    #then
    #    echo "Already has BUSTED results"
    #else
    #    $HYPHY busted --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.BUSTED.json --rates 2 --syn-rates 2 --starting-points 10
    #fi


    #if [ -s ${FILE}.${GENE}.ABSREL.json ] 
    #then
    #    echo "Already has ABSREL results"
    #else
    #    mpirun -np $NP $HYPHYMPI absrel --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.ABSREL.json
    #fi
    
    
fi

}
RunAGene "nsp3" "data/reference_genes/nsp3.fas" "1000" "10000" 0.005 0.0005

RunAGene "RdRp" "data/reference_genes/RdRp.fas" "13000" "17000" 0.01 0.001
RunAGene "helicase" "data/reference_genes/helicase.fas" "15000" "19000" 0.01 0.001
RunAGene "exonuclease" "data/reference_genes/exonuclease.fas" "17000" "21000" 0.01 0.001
RunAGene "endornase" "data/reference_genes/endornase.fas" "19000" "22000" 0.01 0.001
RunAGene "methyltransferase" "data/reference_genes/methyltransferase.fas" "19000" "23000" 0.01 0.001
RunAGene "S" "data/reference_genes/S.fas" "20000" "27000" 0.01 0.001
RunAGene "leader" "data/reference_genes/leader.fas" "0" "2000" 0.005 0.0005
RunAGene "nsp2" "data/reference_genes/nsp2.fas" "0" "10000" 0.005 0.0005
RunAGene "nsp4" "data/reference_genes/nsp4.fas" "7000" "12000" 0.005 0.0005
RunAGene "3C" "data/reference_genes/3C.fas" "8000" "13000" 0.005 0.0005
RunAGene "nsp7" "data/reference_genes/nsp7.fas" "10000" "15000" 0.005 0.0005
RunAGene "nsp8" "data/reference_genes/nsp8.fas" "10000" "15000" 0.005 0.0005
RunAGene "nsp9" "data/reference_genes/nsp9.fas" "10000" "15000" 0.005 0.0005
RunAGene "nsp10" "data/reference_genes/nsp10.fas" "10000" "15000" 0.005 0.0005
RunAGene "M" "data/reference_genes/M.fas" "25000" "30000" 0.01 0.001
RunAGene "N" "data/reference_genes/N.fas" "26000" "35000" 0.01 0.001
RunAGene "ORF3a" "data/reference_genes/ORF3a.fas" "24000" "27000" 0.01 0.001
RunAGene "ORF6" "data/reference_genes/ORF6.fas" "26000" "30000" 0.01 0.001
RunAGene "ORF7a" "data/reference_genes/ORF7a.fas" "26000" "35000" 0.01 0.001
RunAGene "ORF8" "data/reference_genes/ORF8.fas" "26000" "35000" 0.01 0.001
# 
# 
