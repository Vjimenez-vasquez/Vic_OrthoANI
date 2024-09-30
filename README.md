# Vic_OrthoANI
A collection of codes for OrthoANI distances estimations for fasta files
# step 1 : generate a file containing pairwise combinantions
```r
#!/bin/bash
input=$1
ls *.fasta > sequences.txt
ruby -0777 -F'\n' -lane '$F.combination(2) { |c| puts c.join("@")}' sequences.txt > ${1}.txt
cat ${1}.txt
```
# step 2 : estimate OrthoANI distances 
```r
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
sed '1i genome_1 genome_2 distance' a.txt | tr " " "\t" > ${output}.tsv ;
rm a.txt ; 
cat ${output}.tsv ;
ls -lh ;
```
# usage: 
```r
./command0.sh input_name_prefix
./command1.sh input_name_prefix.txt output_name_prefix
```
