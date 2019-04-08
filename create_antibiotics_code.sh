#!/bin/bash  

grep < ARGannot_r3.fasta \> | sed 's/ /;/' | \
    tr ';' '\t' | awk -v OFS='\t' '{ print $5, $4}' | \
    sed 's/Agly/AGly/g' | tee >( sort -k1,2 > code.raw.tsv)\
    >( cut -f1 | sort | uniq > categories.txt)\
    > /dev/null

samples=(
         "AGly:Aminoglucósidos"
         "Bla:Betalactámicos"
         "Colistin:Colistina"
         "Fcd:Ácido_Fusídico"
         "Fcyn:Fosfomicina"
         "Flq:Fluoroquinolonas"
         "Gly:Glicopéptidos"
         "MLS:Macrólidos-Lincosamidas-Estreptograminas"
         "Ntmdz:Nitroimidazol"
         "Oxzln:Oxazolidinonas"
         "Phe:Fenicoles"
         "Rif:Rifampicina"
         "Sul:Sulfonamidas"
         "Tet:Tetraciclinas"
         "Tmt:Trimetroprima"
        )

for i in "${samples[@]}"
do
    key="${i%%:*}"
    val="${i##*:}"
    while read id_cat
    do
        if [ ${key} == ${id_cat} ]; then
            sed "s/$id_cat/$val/g" code.raw.tsv >> test.tsv
        else
            echo -e "$id_cat category not found in: \n${samples[@]}" 
        fi
    done < <(cat categories.txt | cut -f1 | uniq)
done
