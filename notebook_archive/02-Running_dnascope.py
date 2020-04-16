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

# ## Run DNAscope on individual samples

# To compares variant calling that is not affected by the inclusion of bulk data and other cells, we use DNAscope to assess variant calling on individual samples 

# We use a modified input json that we used when we ran DNAseq

sentieon_dnascope_input_dict = {
  "BAM": "",
  "REF": "gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.fasta",
  "OUTPUT_BUCKET": "gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/01-Sentieon_output/DNAScope/{}",
  "BQSR_SITES": "gs://bioskryb-dev-resources-4j2g6d/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz,gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.known_indels.vcf.gz,gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
  "DBSNP": "gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
  "PREEMPTIBLE_TRIES": "2",
  "NONPREEMPTIBLE_TRY": True,
  "STREAM_INPUT": "True",
  "ZONES": "us-central1-a,us-central1-b,us-central1-c,us-central1-f",
  "PROJECT_ID": "bioskryb",
  "EMAIL": "rob.carter@bioskryb.com",
  "CALLING_ALGO": "DNAscope"
}



target_bam_locations = [
    'gs://cromwell_runs/vumc_subsampled_bams/JW-31_merged_n450x10e6.bam',
    'gs://cromwell_runs/vumc_subsampled_bams/JW-23_merged_n450x10e6.bam',
    'gs://cromwell_runs/vumc_subsampled_bams/JW-11_merged_n450x10e6.bam'
]

from copy import deepcopy

for _target_loc in target_bam_locations:
    new_input_dict = deepcopy(sentieon_dnascope_input_dict)
    sample_name = re.sub(".bam", "", os.path.basename(_target_loc))
    new_input_dict['OUTPUT_BUCKET'] = sentieon_dnascope_input_dict['OUTPUT_BUCKET'].format(sample_name)
    new_input_dict['BAM'] = _target_loc
    json_filename = "tracked_data/" + re.sub(".bam", ".sentieoninput.json", os.path.basename(_target_loc))
    with open(json_filename, 'w') as ofh:
        json.dump(new_input_dict, ofh)

json_files = !ls tracked_data/*dnascope.json
json_files

# !ls tracked_data/*.dnascope.json | parallel -j 3 --max-args 1 "tracked_data/time_and_run_sentieon.sh {}"


