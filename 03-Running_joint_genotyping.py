# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.4
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %load_ext autoreload
# %autoreload 2

# +
import re, json, os, logging, io, pprint, subprocess

import pandas as pd

from cromwell_tools.cromwell_api import CromwellAPI as cwt
from cromwell_tools import cromwell_auth
from google.cloud import storage
# -

os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = "/home/rcarter/.google/bioskryb-81ce35d92471.json"

# Here we run joint genotyping using the g.vcfs generated using halotypecaller with bwa aln, but swapping in the g.vcfs from 

# +
# %%writefile tracked_data/gather_vcfs.sh
#!/usr/bin/env bash
BASE_FILENAME=$1

for i in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
do
  echo -I ${INPUT1}/${BASENAME_PREFIX}.${i}.g.vcf.gz
done > files_to_cat.txt
# # cat files_to_cat.txt
gatk GatherVcfs $(cat files_to_cat.txt) -O $OUTPUT1
gatk IndexFeatureFile -F $OUTPUT1
# -

task_list = []
for _sample_name in ["JW-1","JW-10","JW-11","JW-12","JW-15","JW-16","JW-17","JW-18","JW-19","JW-2","JW-20","JW-22","JW-23","JW-25","JW-26","JW-27","JW-28","JW-29","JW-3","JW-30","JW-33","JW-34","JW-36","JW-37","JW-38","JW-5","JW-6","JW-7","JW-8","JW-9","NA12878_Bulk1"]:
    input_folder = "gs://cromwell_runs/vumc_subsampled_bams/{}".format(_sample_name)
    basename_prefix = "{}_merged_n450x10e6".format(_sample_name)
    output_vcf = "{}/{}.g.vcf.gz".format(input_folder,basename_prefix)
    output_vcf_index = "{}/{}.g.vcf.gz.tbi".format(input_folder,basename_prefix)
    task_list.append([_sample_name, basename_prefix, input_folder, output_vcf, output_vcf_index])
pd.DataFrame(task_list, columns = ['--env SAMPLE_NAME', '--env BASENAME_PREFIX',  '--input-recursive INPUT1', '--output OUTPUT1', '--output OUTPUT2']).to_csv("tracked_data/00-gather_vcfs.tsv", sep = "\t", header= True, index = False)

gathervcf_sb_out = subprocess.run(
'/home/rcarter/anacon3/bin/dsub \
--label cost-type=gathervcf \
--retries 3 \
--preemptible 3 \
--wait \
--project bioskryb \
--name gathervcf \
--zones "us-central1-a us-central1-b us-central1-c us-central1-d" \
--logging gs://scratch_space/logs \
--env TMP_DIR=/mnt/data/input/gs/ \
--min-ram 4.0 \
--min-cores 1 \
--disk-size 35 \
--boot-disk-size 15 \
--image gcr.io/bioskryb/gatk:4.1.3.0 \
--ssh \
--enable-stackdriver-monitoring \
--tasks tracked_data/00-gather_vcfs.tsv \
--script tracked_data/gather_vcfs.sh'.format(
        ).split(" "),
    shell = False,
    capture_output=True
)

gathervcf_sb_out

# The GATK's joint genotyping was run on the subsampled data by swapping the sentieon-based g.vcfs in place of the GATK samples for JW-11, JW-23, and JW-31


