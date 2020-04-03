#!/usr/bin/env bash
BASE_FILENAME=$1

for i in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
do
# echo -I ${INPUT1}/${BASENAME_PREFIX}.${i}.g.vcf.gz
done > files_to_cat.txt
# cat files_to_cat.txt
gatk GatherVcfs $(cat files_to_cat.txt) -O $OUTPUT1
gatk IndexFeatureFile -F $OUTPUT1
