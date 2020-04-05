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

# ## Testing Variant callers

# Here we test GATK, Sentieon, and Dragen variant callers to determine what the best mehod for variant calling is in terms of precision and sensitivity

# %load_ext autoreload
# %autoreload 2

# +
import re, json, os, logging, io, pprint, subprocess, requests

import pandas as pd

from cromwell_tools.cromwell_api import CromwellAPI as cwt
from cromwell_tools import cromwell_auth
from google.cloud import storage

from cromwell_functions import defaultCromwellRunner
# -

auth_obj = cromwell_auth.CromwellAuth(url = "http://34.69.50.44:8000", header= None, auth = None)
auth_obj

os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = "/home/rcarter/.google/bioskryb-81ce35d92471.json"

# ## GATK - with Ref

# Here we perform joint genotyping on 16 samples from lots 16 - 20. The samples are not subsampled, but are all similar in file size.
# We first run the GATK gVCFs through GATK's haplotypeCaller, but including the subsampled Bulk2 reference sequence.

non_subsampled_sample_names = list(as_and_wgs_and_dups_metrics_df[as_and_wgs_and_dups_metrics_df["sample_prop"]< 1.06].sort_values('sample_prop', ascending = False).head(16).loc[:, ['sample_name', "sample_prop"]].head(24)['sample_name'])

all_gvcf_files = !gsutil ls gs://bioskryb-cromwell-outputs-jd7fh2/VUMC_Nova183_H37TLDSXY_DUAL_PTA_4249-JW-49-64/**.g.vcf.gz
all_gvcf_files[:2]


gvcf_filenames = !gsutil ls gs://bioskryb-vumc-data/Full_PTA_exome_deepdive_with_other_techs/data/01-gatk_cromwell_output/**.g.vcf.gz
gvcf_filenames_n16 = [_file for _file in gvcf_filenames if re.sub(".g.vcf.gz", "", os.path.basename(_file)) in non_subsampled_sample_names] + ["gs://bioskryb-cromwell-outputs-jd7fh2/temp-folder-for-deletion/NA12878_Bulk2_merged_n450x10e6.g.vcf.gz"]

df_list = [[re.sub("(.g.vcf.gz)|(_merged.+)", "", os.path.basename(_blobname)), _blobname] for _blobname in gvcf_filenames_n16]
pd.DataFrame(df_list).to_csv('tracked_data/gatk_n16_nonsubsampled_sample_map.tsv', sep = "\t", header = False, index= False)
# !gsutil cp tracked_data/gatk_n16_nonsubsampled_sample_map.tsv gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/04-nosubsampling_n16_joint_genotyping/

jj_input_json = cromwell_functions.download_and_read_inputs_json("gs://bioskryb_dev_wdl_and_inputs/gatk-workflows/gatk4-germline-snps-indels/2.0.0/JointGenotyping.hg38.wgs.inputs.json", storage_client)

jj_input_json['JointGenotyping.sample_name_map'] = "gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/04-nosubsampling_n16_joint_genotyping/gatk_n16_nonsubsampled_sample_map.tsv"
jj_input_json['JointGenotyping.callset_name'] = "nosubsampling_n16_gatk_joint_genotyping_gatk"

# +
project_name = "ComparingGatkSentieonDRAGEN"
final_output_bucket = "gs://bioskryb-work-d8f6s9"
cromwell_runs_bucket = "gs://bioskryb-dev-cromwell-runs-8dhf5d"

    
output_base_location = "{}/{}/{}".format(final_output_bucket, project_name, "04-nosubsampling_n16_joint_genotyping")
jj_options_dict = {
        "read_from_cache": False,
        "default_runtime_attributes": {
                "zones": "us-central1-a us-central1-b us-central1-c us-central1-f"
        },
        "final_workflow_outputs_dir": output_base_location,
        "use_relative_output_paths": True,
        "final_call_logs_dir": "{}/call_logs".format(output_base_location),
        "jes_gcs_root": cromwell_runs_bucket,
        "google_labels": {
                "pipeline-name": "gatk4-germline-snps-indels",
                "project-name": "comparing-gatk-sentieon-dragen"
        }
}
# -

input_iobytes = io.BytesIO(json.dumps(jj_input_json).encode())
options_iobytes = io.BytesIO(json.dumps(jj_options_dict).encode())
jj_resp = cwt.submit(auth_obj, 
       wdl_file = cromwell_functions.get_wdl_iobytes("gs://bioskryb_dev_wdl_and_inputs/gatk-workflows/gatk4-germline-snps-indels/2.0.0/JointGenotyping.wdl", storage_client), 
       inputs_files = input_iobytes,
       options_file = options_iobytes
)

jj_resp.content

gvcf_filenames_n16

# ## GATK - no Ref

# Here we perform joint genotyping on 16 samples from lots 16 - 20. The samples are not subsampled, but are all similar in file size.
# We exclude the bulk data from this run

non_subsampled_sample_names = list(as_and_wgs_and_dups_metrics_df[as_and_wgs_and_dups_metrics_df["sample_prop"]< 1.06].sort_values('sample_prop', ascending = False).head(16).loc[:, ['sample_name', "sample_prop"]].head(24)['sample_name'])

all_gvcf_files = !gsutil ls gs://bioskryb-cromwell-outputs-jd7fh2/VUMC_Nova183_H37TLDSXY_DUAL_PTA_4249-JW-49-64/**.g.vcf.gz
all_gvcf_files[:2]


gvcf_filenames = !gsutil ls gs://bioskryb-vumc-data/Full_PTA_exome_deepdive_with_other_techs/data/01-gatk_cromwell_output/**.g.vcf.gz
gvcf_filenames_n16_no_ref = [_file for _file in gvcf_filenames if re.sub(".g.vcf.gz", "", os.path.basename(_file)) in non_subsampled_sample_names] 

df_list = [[re.sub("(.g.vcf.gz)|(_merged.+)", "", os.path.basename(_blobname)), _blobname] for _blobname in gvcf_filenames_n16]
pd.DataFrame(df_list).to_csv('tracked_data/gatk_n16_nonsubsampled_sample_map.tsv', sep = "\t", header = False, index= False)
# !gsutil cp tracked_data/gatk_n16_nonsubsampled_sample_map.tsv gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/04-nosubsampling_n16_joint_genotyping/

jj_input_json = cromwell_functions.download_and_read_inputs_json("gs://bioskryb_dev_wdl_and_inputs/gatk-workflows/gatk4-germline-snps-indels/2.0.0/JointGenotyping.hg38.wgs.inputs.json", storage_client)

jj_input_json['JointGenotyping.sample_name_map'] = "gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/04-nosubsampling_n16_joint_genotyping/gatk_n16_nonsubsampled_sample_map.tsv"
jj_input_json['JointGenotyping.callset_name'] = "nosubsampling_n16_gatk_joint_genotyping_gatk"

# +
project_name = "ComparingGatkSentieonDRAGEN"
final_output_bucket = "gs://bioskryb-work-d8f6s9"
cromwell_runs_bucket = "gs://bioskryb-dev-cromwell-runs-8dhf5d"

    
output_base_location = "{}/{}/{}".format(final_output_bucket, project_name, "04-nosubsampling_n16_joint_genotyping")
jj_options_dict = {
        "read_from_cache": False,
        "default_runtime_attributes": {
                "zones": "us-central1-a us-central1-b us-central1-c us-central1-f"
        },
        "final_workflow_outputs_dir": output_base_location,
        "use_relative_output_paths": True,
        "final_call_logs_dir": "{}/call_logs".format(output_base_location),
        "jes_gcs_root": cromwell_runs_bucket,
        "google_labels": {
                "pipeline-name": "gatk4-germline-snps-indels",
                "project-name": "comparing-gatk-sentieon-dragen"
        }
}
# -

input_iobytes = io.BytesIO(json.dumps(jj_input_json).encode())
options_iobytes = io.BytesIO(json.dumps(jj_options_dict).encode())
jj_resp = cwt.submit(auth_obj, 
       wdl_file = cromwell_functions.get_wdl_iobytes("gs://bioskryb_dev_wdl_and_inputs/gatk-workflows/gatk4-germline-snps-indels/2.0.0/JointGenotyping.wdl", storage_client), 
       inputs_files = input_iobytes,
       options_file = options_iobytes
)

jj_resp.content

gvcf_filenames_n16

# ## Sentieon

# The same 16 samples from above (plus the reference) are joint gneotyped using Sentieon

# First, generate gVCFs

r1_fastq_files = !gsutil ls gs://bioskryb-vumc-data/Nova181_H2VMKDSXY/*_R1_*fastq.gz
r2_fastq_files = !gsutil ls gs://bioskryb-vumc-data/Nova182_H2VGNDSXY/*_R1_*fastq.gz
all_r1_fastqs = r1_fastq_files + r2_fastq_files

all_r1_fastqs
valid_fastq_r1_filenames_n16 = [_file for _file in all_r1_fastqs if re.sub("_S[01-9].*.fastq.gz", "", os.path.basename(_file)) in non_subsampled_sample_names] + ["gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/NA12878_Bulk2_merged_n450x10e6_R1_001.fastq.gz"]
valid_fastq_r1_filenames_n16 = sorted(valid_fastq_r1_filenames_n16)

sample_to_fastq_dict = {}
for _file in valid_fastq_r1_filenames_n16:
    sample_name = re.sub("_S[01-9].*.fastq.gz", "", os.path.basename(_file))
    try:
        sample_to_fastq_dict[sample_name].append(_file)
    except:
        sample_to_fastq_dict[sample_name] = [_file]

sentieon_input_dict = {
  "FQ1": "",
  "FQ2": "",
  "REF": "gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.fasta",
  "READGROUP": "",
  "OUTPUT_BUCKET": "gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/03-joint_genotyping_input_gvcfs/sentieon/{}",
  "BQSR_SITES": "gs://bioskryb-dev-resources-4j2g6d/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz,gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.known_indels.vcf.gz,gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
  "DBSNP": "gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
  "GVCF_OUTPUT": "true",
  "PREEMPTIBLE_TRIES": "2",
  "NONPREEMPTIBLE_TRY": True,
  "STREAM_INPUT": "True",
  "ZONES": "us-central1-a,us-central1-b,us-central1-c,us-central1-f",
  "PROJECT_ID": "bioskryb",
  "EMAIL": "rob.carter@bioskryb.com"
}

for _sample_name, file_list in sample_to_fastq_dict.items():
    rg_list = []
    for _i in range(len(file_list)):
        rg_list.append("@RG\\tID:rgid-{}\\tSM:{}\\tPL:ILLUMINA".format(_i, _sample_name))
    rg = ",".join(rg_list)
    fq1 = ",".join(file_list)
    fq2 = re.sub("_R1_", "_R2_", ",".join(file_list))
    new_input_dict = deepcopy(sentieon_input_dict)
    new_input_dict['OUTPUT_BUCKET'] = sentieon_input_dict['OUTPUT_BUCKET'].format(_sample_name)
    new_input_dict['FQ1'] = fq1
    new_input_dict['FQ2'] = fq2
    new_input_dict['READGROUP'] = rg
    json_filename = "tracked_data/{}.n16.sentieoninput.json".format(_sample_name)
    with open(json_filename, 'w') as ofh:
        json.dump(new_input_dict, ofh)


# %%writefile tracked_data/time_and_run_sentieon.sh
#!/usr/bin/env bash
JSON_FILE=$1
LOG_FILE=${JSON_FILE}.run.log
# echo Start: $(date +'%m-%d-%Y-%H-%M-%S') > $LOG_FILE
# echo Running: ~/sentieon-google-genomics/runner/sentieon_runner.py $JSON_FILE
# echo $(~/sentieon-google-genomics/runner/sentieon_runner.py $JSON_FILE) >> $LOG_FILE
# echo End: $(date +'%m-%d-%Y-%H-%M-%S') >> $LOG_FILE

# !chmod 755 tracked_data/time_and_run_sentieon.sh

json_files = !ls tracked_data/*.n16.sentieoninput.json
json_files

# !ls tracked_data/*JW-32.n16.sentieoninput.json | parallel -j 16 --max-args 1 "tracked_data/time_and_run_sentieon.sh {}"

options_dict = {
        "default_runtime_attributes": {
                "zones": "us-central1-a us-central1-b us-central1-c us-central1-f",
                "noAddress": False,
                "bootDiskSizeGb": 60
        },
}

# +
default_inputs_dict = requests.get("https://raw.githubusercontent.com/gatk-workflows/generic-wdl-script/master/generic-workflow.input.json").json()
ref_sequence = "gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.fasta"
ref_sequence_index = "gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.fasta.fai"
dbsnp_filename = "gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
dbsnp_index = filename = "gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi"
local_dbsnp_filename = os.path.basename(dbsnp_filename)
local_ref = os.path.basename(ref_sequence)
local_ref_index = os.path.basename(ref_sequence_index)

bam = "gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/03-joint_genotyping_input_gvcfs/sentieon/4249-JW-32/aligned_reads/dedup.bam"
bai = "gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/03-joint_genotyping_input_gvcfs/sentieon/4249-JW-32/aligned_reads/dedup.bam.bai"
recal = "gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/03-joint_genotyping_input_gvcfs/sentieon/4249-JW-32/aligned_reads/recal_data.table"

haplotyper_command = "sentieon driver -t 64 -r {ref} -i {bam} -q {recal} --algo Haplotyper -d {dbsnp} --emit_mode gvcf 4249-JW-32_hc.g.vcf.gz && mv 4249-JW-32_hc.g.vcf.gz* ./output/".format(ref = os.path.basename(ref_sequence),  bam = os.path.basename(bam), recal = os.path.basename(recal), dbsnp = os.path.basename(dbsnp_filename))

wdl_iobytes = io.BytesIO(requests.get("https://raw.githubusercontent.com/gatk-workflows/generic-wdl-script/master/generic-workflow.wdl").text.encode())
new_sentieon_inputs_dict = deepcopy(default_inputs_dict)
new_sentieon_inputs_dict["genericworkflow.GenericTask.shell_command"] = haplotyper_command
new_sentieon_inputs_dict["genericworkflow.GenericTask.input_files"] =[ref_sequence, ref_sequence_index, dbsnp_filename, dbsnp_filename_index, bam, bai, recal]
new_sentieon_inputs_dict["genericworkflow.docker_override"] = "gcr.io/bioskryb/sentieon-201911:latest"
new_sentieon_inputs_dict["genericworkflow.GenericTask.machine_mem_gb"] = 84
new_sentieon_inputs_dict["genericworkflow.GenericTask.disk_space_gb"] = 400
new_sentieon_inputs_dict["genericworkflow.GenericTask.cpu"] = 64
wdl_sentieon_inputs_iobytes = io.BytesIO(json.dumps(new_sentieon_inputs_dict).encode())
temp_sentieon_resp = cwt.submit(auth = auth_obj, wdl_file=wdl_iobytes, inputs_files=wdl_sentieon_inputs_iobytes, options_file= io.BytesIO(json.dumps(options_dict).encode()))
# -

temp_sentieon_resp.content

# ### Joint-genotyping

gvcf_sentieon_filenames = !gsutil ls gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/03-joint_genotyping_input_gvcfs/sentieon/variants/**.g.vcf.gz
gvcf_sentieon_filenames_n16 = gvcf_sentieon_filenames
gvcf_sentieon_filenames_n16

len(gvcf_sentieon_filenames_n16)

# +
project_name = "ComparingGatkSentieonDRAGEN"
final_output_bucket = "gs://bioskryb-work-d8f6s9"
cromwell_runs_bucket = "gs://bioskryb-dev-cromwell-runs-8dhf5d"

    
output_base_location = "{}/{}/{}".format(final_output_bucket, project_name, "data/04-nosubsampling_n16_joint_genotyping/sentieon")
jj_sentieon_options_dict = {
        "read_from_cache": False,
        "default_runtime_attributes": {
                "zones": "us-central1-a us-central1-b us-central1-c us-central1-f",
                "noAddress": False
        },
        "final_workflow_outputs_dir": output_base_location,
        "use_relative_output_paths": True,
        "final_call_logs_dir": "{}/call_logs".format(output_base_location),
        "jes_gcs_root": cromwell_runs_bucket,
        "google_labels": {
                "pipeline-name": "sentieon-gvcftyper",
                "project-name": "comparing-gatk-sentieon-dragen"
        }
}

# +
vcf_df = pd.DataFrame([gvcf_sentieon_filenames_n16], columns= [re.sub("_hc.g.vcf.gz", "_vcf=input", os.path.basename(_filename)) for _filename in gvcf_sentieon_filenames_n16])

gvcf_sentieon_index_filenames_n16 = ["{}.tbi".format(_filename) for  _filename in gvcf_sentieon_filenames_n16]
vcf_index_df = pd.DataFrame([gvcf_sentieon_index_filenames_n16], columns= [re.sub("_hc.g.vcf.gz.tbi", "_vcf_index=input", os.path.basename(_filename)) for _filename in gvcf_sentieon_index_filenames_n16])

ref_sequence = "gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.fasta"
ref_sequence_index = "gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.fasta.fai"
dbsnp_filename = "gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
dbsnp_filename_index = "gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi"
other_input_files_df = pd.DataFrame([[ref_sequence, ref_sequence_index, dbsnp_filename, dbsnp_filename_index]], columns=["ref=input", "ref_index=input", "dbsnp=input", "dbsnp_index=input"])
outputs_df = pd.DataFrame([["nosubsampling_n16_sentieon_joint_genotyping_sentieon.vcf.gz", "nosubsampling_n16_sentieon_joint_genotyping_sentieon.vcf.gz.tbi"]], columns = ["out_vcf=output", "out_vcf_index=output"])

inputs_df = pd.concat([vcf_df, vcf_index_df, other_input_files_df, outputs_df], axis = 1)

variant_subcommand = "{" + "} -v {".join(_colname.split("=")[0] for _colname in  vcf_df.columns) + "}"
gvcftyper_command_template = "sentieon driver -r {ref} -t 64  --algo GVCFtyper -v " + variant_subcommand + " -d {dbsnp}"
gvcftyper_command_template

# +
options_dict = {
        "default_runtime_attributes": {
                "zones": "us-central1-a us-central1-b us-central1-c us-central1-f",
                "noAddress": False,
                "bootDiskSizeGb": 20
        },
}

inputs_dict = requests.get("https://raw.githubusercontent.com/gatk-workflows/generic-wdl-script/master/generic-workflow.input.json").json()
inputs_dict["genericworkflow.docker_override"] = "gcr.io/bioskryb/sentieon-201911:latest"
inputs_dict["genericworkflow.GenericTask.machine_mem_gb"] = 64
inputs_dict["genericworkflow.GenericTask.disk_space_gb"] = 100
inputs_dict["genericworkflow.GenericTask.cpu"] = 64

gctyper_dcr = defaultCromwellRunner(auth_obj, gvcftyper_command_template, inputs_df, options_dict, inputs_dict, wdl_text= requests.get("https://raw.githubusercontent.com/gatk-workflows/generic-wdl-script/master/generic-workflow.wdl").text)
# -

gctyper_dcr.submit_jobs()

gctyper_dcr.get_status()

# ## Testing Joint Genotyping methods

# ### Testing Sentieon's Joint genotyper on sentieon data

gvcf_sentieon_filenames_n16

# Sentieon's output files were renamed at the command line:
#
# `for i in $(gsutil ls gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/03-joint_genotyping_input_gvcfs/sentieon/**.vcf.gz*); do NEW_NAME=$(dirname $(dirname $i))_$(basename $i); gsutil mv $i $NEW_NAME; done`
#
# `for i in $(gsutil ls gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/03-joint_genotyping_input_gvcfs/sentieon/**.vcf.gz*); do NEW_NAME=$(dirname $i)/variants/$(basename $i); gsutil mv $i $NEW_NAME; done`

gvcf_updated_sentieon_filenames = !gsutil ls gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/03-joint_genotyping_input_gvcfs/sentieon**.g.vcf.gz
gvcf_updated_sentieon_filenames_n16 = gvcf_updated_sentieon_filenames + ["gs://bioskryb-cromwell-outputs-jd7fh2/temp-folder-for-deletion/NA12878_Bulk2_merged_n450x10e6.g.vcf.gz"]
gvcf_updated_sentieon_filenames_n16

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
new_sentieon_inputs_dict["genericworkflow.GenericTask.input_files"] = gvcf_updated_sentieon_filenames_n16[:2]  + [ref_sequence, ref_sequence_index, dbsnp_filename]
new_sentieon_inputs_dict["genericworkflow.docker_override"] = "gcr.io/bioskryb/sentieon-201911:latest"
new_sentieon_inputs_dict["genericworkflow.GenericTask.machine_mem_gb"] = 64
new_sentieon_inputs_dict["genericworkflow.GenericTask.disk_space_gb"] = 100
new_sentieon_inputs_dict["genericworkflow.GenericTask.cpu"] = 64
wdl_sentieon_inputs_iobytes = io.BytesIO(json.dumps(new_sentieon_inputs_dict).encode())
temp_sentieon_resp = cwt.submit(auth = auth_obj, wdl_file=wdl_iobytes, inputs_files=wdl_sentieon_inputs_iobytes)
# -

temp_sentieon_resp.content

cwt.abort(temp_sentieon_resp.json()['id'], auth_obj).json()


