#!/bin/bash

input=$1
ls *.fasta > sequences.txt
ruby -0777 -F'\n' -lane '$F.combination(2) { |c| puts c.join("@")}' sequences.txt > ${1}.txt
cat ${1}.txt



