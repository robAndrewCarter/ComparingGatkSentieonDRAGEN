#!/usr/bin/env bash
echo $2 | sed 's/,/\n/g' | parallel --max-args=1 --jobs 8 /home/rcarter/rtg-tools/bin/rtg-tools-3.10.1-4d58eadb/rtg RTG_MEM=6G vcfeval -b $1 -c $1 --sample=NA12878_Bulk2,{} -o scratch/NA12878_Bulk2_vs_{}-${3} -t /home/rcarter/R_bioskryb_na12878_analysis/NA12878_comparison/SDF/ --squash-ploidy
