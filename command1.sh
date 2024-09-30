#!/bin/bash

input=$1
output=$2

for line in $(cat $input)
do
genome1=$(cut -f1 -d'@' <<< "$line")
genome2=$(cut -f2 -d'@' <<< "$line")
res=$(orthoani -q $genome1 -r $genome2 -j 20)
echo $genome1 $genome2 $res
done > a.txt

sed '1i genome_1 genome_2 distance' a.txt | tr "" "\t" > ${output}.tsv ;
rm a.txt ; 
cat ${output}.tsv ;
ls -lh ;

