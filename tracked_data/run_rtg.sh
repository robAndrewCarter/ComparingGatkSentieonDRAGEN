# baseline vcf
# baseline sample_name
# comparison vcf
# joined comparison vcfs
#filder suffix
# extra options
#!/usr/bin/env bash
echo $4 | sed 's/,/\n/g' | parallel --max-args=1 --jobs 8 /home/rcarter/rtg-tools/bin/rtg-tools-3.10.1-4d58eadb/rtg RTG_MEM=6G vcfeval -b $1 -c $3 --sample=${2},{} -o scratch/${2}_vs_{}-${5} -e /home/rcarter/Homo_sapiens_assembly38_n25chr.bed -t /home/rcarter/R_bioskryb_na12878_analysis/NA12878_comparison/SDF/  --squash-ploidy $6
