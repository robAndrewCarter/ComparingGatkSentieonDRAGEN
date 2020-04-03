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
import re, json, os, logging, io, pprint, subprocess, requests
from copy import deepcopy

import pandas as pd

from cromwell_tools.cromwell_api import CromwellAPI as cwt
from cromwell_tools import cromwell_auth
from google.cloud import storage

import cromwell_functions
# -

storage_client = storage.Client('bioskryb')

os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = "/home/rcarter/.google/bioskryb-81ce35d92471.json"

auth_obj = cromwell_auth.CromwellAuth(url = "http://34.69.50.44:8000", header= None, auth = None)
auth_obj

gvcf_updated_sentieon_filenames = !gsutil ls gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/03-joint_genotyping_input_gvcfs/sentieon**.g.vcf.gz
gvcf_updated_sentieon_filenames_n16 = gvcf_updated_sentieon_filenames + ["gs://bioskryb-cromwell-outputs-jd7fh2/temp-folder-for-deletion/NA12878_Bulk2_merged_n450x10e6.g.vcf.gz"]
gvcf_updated_sentieon_filenames_n16

options_dict = {
        "read_from_cache": False,
        "default_runtime_attributes": {
                "zones": "us-central1-a us-central1-b us-central1-c us-central1-f",
                "noAddress": False
        },
        "final_workflow_outputs_dir": "gs://scratch_space",
        "use_relative_output_paths": False,
        "final_call_logs_dir": "{}/call_logs".format('gs://scratch_space'),
        "google_labels": {
                "pipeline-name": "gatk4",
                "project-name": "delme"
        }
}

# +
default_inputs_dict = requests.get("https://raw.githubusercontent.com/gatk-workflows/generic-wdl-script/master/generic-workflow.input.json").json()
ref_sequence = "gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.fasta"
ref_sequence_index = "gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.fasta.fai"
dbsnp_filename = "gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
local_dbsnp_filename = os.path.basename(dbsnp_filename)
local_ref = os.path.basename(ref_sequence)
local_ref_index = os.path.basename(ref_sequence_index)

driver_options = "-t 64"
algo_options = "-d {}".format(local_dbsnp_filename)
gvcf_string = " -v ".join([os.path.basename(_path) for _path in gvcf_updated_sentieon_filenames_n16[:2]])

sentieon_jj_command = "printenv && sentieon driver -r {ref} {driver_options} --algo GVCFtyper -v {vcf_string} {algo_options} ./output/output.vcf.gz".format(ref = local_ref, driver_options = driver_options, vcf_string = gvcf_string, algo_options = algo_options)

wdl_iobytes = io.BytesIO(requests.get("https://raw.githubusercontent.com/gatk-workflows/generic-wdl-script/master/generic-workflow.wdl").text.encode())
new_sentieon_inputs_dict = deepcopy(default_inputs_dict)
new_sentieon_inputs_dict["genericworkflow.GenericTask.shell_command"] = sentieon_jj_command
new_sentieon_inputs_dict["genericworkflow.GenericTask.input_files"] = ["{}.tbi".format(i) for i in gvcf_updated_sentieon_filenames_n16[:2]] + gvcf_updated_sentieon_filenames_n16[:2]  + [ref_sequence, ref_sequence_index, dbsnp_filename, dbsnp_filename + ".tbi"]
new_sentieon_inputs_dict["genericworkflow.docker_override"] = "gcr.io/bioskryb/sentieon-201911:latest"
new_sentieon_inputs_dict["genericworkflow.GenericTask.machine_mem_gb"] = 32
new_sentieon_inputs_dict["genericworkflow.GenericTask.disk_space_gb"] = 100
new_sentieon_inputs_dict["genericworkflow.GenericTask.cpu"] = 32
wdl_sentieon_inputs_iobytes = io.BytesIO(json.dumps(new_sentieon_inputs_dict).encode())
wdl_sentieon__iobytes = io.BytesIO(json.dumps(new_sentieon_inputs_dict).encode())
temp_sentieon_resp = cwt.submit(auth = auth_obj, wdl_file=wdl_iobytes, inputs_files=wdl_sentieon_inputs_iobytes, options_file=io.BytesIO(json.dumps(options_dict).encode()))
# -

cwt.metadata(temp_sentieon_resp.json()['id'], auth_obj).json()



new_sentieon_inputs_dict


