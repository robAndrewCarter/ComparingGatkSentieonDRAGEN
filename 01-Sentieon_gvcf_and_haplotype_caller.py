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

# ## Run Sentieon on 3 samples

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

# +
## task_list = []
## for _bam_location in target_bam_locations:
##     output_prefix = "gs://bioskryb-illumina-share-3j5h2s/"
##     read1_out = "{}{}".format(output_prefix, os.path.basename(re.sub(".bam", "_R1.fastq.gz", _bam_location)))
##     read2_out = "{}{}".format(output_prefix, os.path.basename(re.sub(".bam", "_R2.fastq.gz", _bam_location)))
##     task_list.append([_bam_location, read1_out, read2_out])
## pd.DataFrame(task_list, columns = ['--input INPUT1', '--output OUTPUT1', '--output OUTPUT2']).to_csv("sam_to_fastq.tsv", sep = "\t", header= True, index = False)
# -

sentieon_input_dict = {
  "BAM": "",
  "REF": "gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.fasta",
  "OUTPUT_BUCKET": "gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/01-Sentieon_output/{}",
  "BQSR_SITES": "gs://bioskryb-dev-resources-4j2g6d/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz,gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.known_indels.vcf.gz,gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
  "DBSNP": "gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
  "PREEMPTIBLE_TRIES": "2",
  "NONPREEMPTIBLE_TRY": True,
  "STREAM_INPUT": "True",
  "ZONES": "us-central1-a,us-central1-b,us-central1-c,us-central1-f",
  "PROJECT_ID": "bioskryb",
  "EMAIL": "rob.carter@bioskryb.com"
}

target_bam_locations = [
    'gs://cromwell_runs/vumc_subsampled_bams/JW-31_merged_n450x10e6.bam',
    'gs://cromwell_runs/vumc_subsampled_bams/JW-23_merged_n450x10e6.bam',
    'gs://cromwell_runs/vumc_subsampled_bams/JW-11_merged_n450x10e6.bam'
]

# The sentieon runner script does not handle bai files of the form foo.bai and instead requires foo.bam.bai. Therefore, these will need to be created for the threee bams to test. They will be subsequently deleted.

for _bam in target_bam_locations:
    bai_filename = re.sub("\.bam", ".bai", _bam)
    new_bai_filename = _bam + ".bai"
#     !gsutil cp $bai_filename $new_bai_filename

from copy import deepcopy

# %%writefile tracked_data/time_and_run_sentieon.sh
#!/usr/bin/env bash
JSON_FILE=$1
LOG_FILE=${JSON_FILE}.run.log
# echo Start: $(date +'%m-%d-%Y-%H-%M-%S') > $LOG_FILE
# echo Running: ~/sentieon-google-genomics/runner/sentieon_runner.py $JSON_FILE
# echo $(~/sentieon-google-genomics/runner/sentieon_runner.py $JSON_FILE) >> $LOG_FILE
# echo End: $(date +'%m-%d-%Y-%H-%M-%S') >> $LOG_FILE

# !chmod 755 tracked_data/time_and_run_sentieon.sh

for _target_loc in target_bam_locations:
    new_input_dict = deepcopy(sentieon_input_dict)
    sample_name = re.sub(".bam", "", os.path.basename(_target_loc))
    new_input_dict['OUTPUT_BUCKET'] = sentieon_input_dict['OUTPUT_BUCKET'].format(sample_name)
    new_input_dict['BAM'] = _target_loc
    json_filename = "tracked_data/" + re.sub(".bam", ".sentieoninput.json", os.path.basename(_target_loc))
    with open(json_filename, 'w') as ofh:
        json.dump(new_input_dict, ofh)

json_files = !ls tracked_data/*.json
json_files

# !ls tracked_data/*.json | parallel -j 3 --max-args 1 "tracked_data/time_and_run_sentieon.sh {}"

# ## Run joint genotyping with Sentieon samples in place of GATK samples

# To accurately compare the Sentieon GATK replacement with GATK, we need to rerun the joint genotyping with the same samples as were run in GATK, but swapping out the three samples that were run SENTIEON. This will cause differences in variant calling to be causesd only by the presence of the three Sentieon-run cells


