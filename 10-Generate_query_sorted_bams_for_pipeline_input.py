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

# # Subsample bams

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

# ## Subsample bams

# ### Collect the WGS and AS metrics files from the full genome GATK run to determine subsampling level

lot16_20_dup_metrics_filenames = !gsutil ls gs://bioskryb-vumc-data/Full_PTA_exome_deepdive_with_other_techs/data/01-gatk_cromwell_output/**.duplicate_metrics | grep -vP "JW-[1-9]/" | grep -vP "JW-7[01-9]" | grep -v SRR | grep -v ".readgroup.alignment"
lot16_20_dup_metrics_filenames[:2]

# +
# #!gsutil ls gs://bioskryb-vumc-data/Full_PTA_exome_deepdive_with_other_techs/data/01-gatk_cromwell_output/**metrics*
# -

lot16_20_as_metrics_filenames = !gsutil ls gs://bioskryb-vumc-data/Full_PTA_exome_deepdive_with_other_techs/data/01-gatk_cromwell_output/**alignment_summary_metrics | grep -vP "JW-[1-9]/" | grep -vP "JW-7[01-9]" | grep -v SRR | grep -v ".readgroup.alignment"
lot16_20_as_metrics_filenames[:2]

lot16_20_wgs_metrics_filenames = !gsutil ls gs://bioskryb-vumc-data/Full_PTA_exome_deepdive_with_other_techs/data/01-gatk_cromwell_output/**.wgs_metrics | grep -vP "JW-[1-9]/" | grep -vP "JW-7[01-9]" | grep -v SRR
lot16_20_wgs_metrics_filenames[:2]

# + jupyter={"outputs_hidden": true}
df_list = []    
for _dup_metrics_filename in lot16_20_dup_metrics_filenames:
    temp_blob = cromwell_functions.get_blob_from_gcp_location(_dup_metrics_filename, storage_client)
    sample_name = re.sub(".duplicate_metrics", "", os.path.basename(_dup_metrics_filename))
    temp_df = pd.read_table(io.BytesIO(temp_blob.download_as_string()), sep = "\t", skiprows= 6, nrows=1)
    temp_df['sample_name'] = sample_name
    df_list.append(temp_df)
dup_metrics_df = pd.concat(df_list)
dup_metrics_df
# -

df_list = []    
for _wgs_metrics_filename in lot16_20_wgs_metrics_filenames:
    temp_blob = cromwell_functions.get_blob_from_gcp_location(_wgs_metrics_filename, storage_client)
    sample_name = re.sub(".wgs_metrics", "", os.path.basename(_wgs_metrics_filename))
    temp_df = pd.read_table(io.BytesIO(temp_blob.download_as_string()), sep = "\t", skiprows= 6, nrows=1)
    temp_df['sample_name'] = sample_name
    df_list.append(temp_df)
wgs_metrics_df = pd.concat(df_list)
wgs_metrics_df

# + jupyter={"outputs_hidden": true}
df_list = []    
for _as_metrics_filename in lot16_20_as_metrics_filenames:
    temp_blob = cromwell_functions.get_blob_from_gcp_location(_as_metrics_filename, storage_client)
    sample_name = re.sub(".alignment_summary_metrics", "", os.path.basename(_as_metrics_filename))
    temp_df = pd.read_table(io.BytesIO(temp_blob.download_as_string()), sep = "\t", skiprows= 6, nrows=3).tail(1)
    temp_df['sample_name'] = sample_name
    df_list.append(temp_df)
as_metrics_df = pd.concat(df_list)
# -

as_metrics_df.head()

# + jupyter={"outputs_hidden": true}
as_and_wgs_and_dups_metrics_df = pd.merge(pd.merge(wgs_metrics_df, as_metrics_df, on = "sample_name"), dup_metrics_df, on = 'sample_name')
as_and_wgs_and_dups_metrics_df['read_ids_examined'] = (as_and_wgs_and_dups_metrics_df['UNPAIRED_READS_EXAMINED']+ as_and_wgs_and_dups_metrics_df["READ_PAIRS_EXAMINED"])
as_and_wgs_and_dups_metrics_df['read_ids_with_dups'] = (as_and_wgs_and_dups_metrics_df['UNPAIRED_READ_DUPLICATES'] + as_and_wgs_and_dups_metrics_df['READ_PAIR_DUPLICATES'])
as_and_wgs_and_dups_metrics_df['read_ids_without_dups'] = as_and_wgs_and_dups_metrics_df['read_ids_examined'] - as_and_wgs_and_dups_metrics_df['read_ids_with_dups']
as_and_wgs_and_dups_metrics_df['sample_prop'] =   (4.25e8 / 2) / as_and_wgs_and_dups_metrics_df['read_ids_without_dups']
as_and_wgs_and_dups_metrics_df.loc[:, ['sample_name', 'read_ids_examined', 'read_ids_with_dups', 'read_ids_without_dups', 'sample_prop']]
# -

sample_and_sample_prop_df = as_and_wgs_and_dups_metrics_df[as_and_wgs_and_dups_metrics_df["sample_prop"]< 1.0].sort_values('sample_prop', ascending = False).head(24).loc[:, ['sample_name', "sample_prop"]]
sample_and_sample_prop_df

# ### Run Subsampling

cram_files = !gsutil ls gs://bioskryb-vumc-data/Full_PTA_exome_deepdive_with_other_techs/data/01-gatk_cromwell_output/**.cram | grep -vP "JW-[1-9]/" | grep -vP "JW-7[01-9]" | grep -v SRR
cram_files[:2]

cram_info_df = pd.DataFrame({"cram_location": cram_files, "sample_name": [re.sub(".cram", "", os.path.basename(_file)) for _file in cram_files]})

all_cromwell_info_df = pd.merge(sample_and_sample_prop_df, cram_info_df, how = "left", on = "sample_name")

all_cromwell_info_df

inputs_dict = requests.get("https://raw.githubusercontent.com/gatk-workflows/generic-wdl-script/master/generic-workflow.input.json").json()
ref_sequence = "gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.fasta"
ref_sequence_index = "gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.fasta.fai"
commands_list = []
resp_list = []
for _row_index, _row in all_cromwell_info_df.iterrows():
    wdl_iobytes = io.BytesIO(requests.get("https://raw.githubusercontent.com/gatk-workflows/generic-wdl-script/master/generic-workflow.wdl").text.encode())
    output_filename = "{}.query_sorted.n425e6.bam".format(_row['sample_name'])
    target_cram = _row['cram_location']
    target_cram_index = _row['cram_location'] + ".crai"
    local_cram = os.path.basename(target_cram)
    local_cram_index = os.path.basename(target_cram_index)
    local_ref = os.path.basename(ref_sequence)
    local_ref_index = os.path.basename(ref_sequence_index)
    subsampling_command = "samtools view -T {ref} -b  -s {sp} {tc} | samtools sort -n -o {ob} && mv {ob} ./output/".format(ref = local_ref, sp = _row["sample_prop"], tc = local_cram, ob = output_filename)
    commands_list.append(subsampling_command)
    new_inputs_dict = deepcopy(inputs_dict)
    new_inputs_dict["genericworkflow.GenericTask.shell_command"] = subsampling_command
    new_inputs_dict["genericworkflow.GenericTask.input_files"] =  [target_cram, target_cram_index, ref_sequence, ref_sequence_index]
    new_inputs_dict["genericworkflow.docker_override"] = "gcr.io/bioskryb/bwa_and_samtools:latest"
    new_inputs_dict["genericworkflow.GenericTask.machine_mem_gb"] = 6
    new_inputs_dict["genericworkflow.GenericTask.disk_space_gb"] = 100
    wdl_inputs_iobytes = io.BytesIO(json.dumps(new_inputs_dict).encode())
    temp_resp = cwt.submit(auth = auth_obj, wdl_file=wdl_iobytes, inputs_files=wdl_inputs_iobytes, options_file = io.BytesIO(json.dumps(options_dict).encode()))
    resp_list.append([_row["sample_name"], temp_resp.json()['id']])

# ### Grab the reference sequence

wdl_iobytes = io.BytesIO(requests.get("https://raw.githubusercontent.com/gatk-workflows/generic-wdl-script/master/generic-workflow.wdl").text.encode())
output_filename = "{}.query_sorted.n425e6.bam".format("NA12878_Bulk2")
target_bam = "gs://cromwell_runs/vumc_subsampled_bams/NA12878_Bulk2_merged_n450x10e6.bam"
target_bam_index = "gs://cromwell_runs/vumc_subsampled_bams/NA12878_Bulk2_merged_n450x10e6.bai"
local_bam = os.path.basename(target_bam)
local_bam_index = os.path.basename(target_bam_index)
local_ref = os.path.basename(ref_sequence)
local_ref_index = os.path.basename(ref_sequence_index)
subsampling_command = "samtools view -b {tc} | samtools sort -n -o {ob} && mv {ob} ./output/".format(ref = local_ref, tc = local_bam, ob = output_filename)
commands_list.append(subsampling_command)
new_inputs_dict = deepcopy(inputs_dict)
new_inputs_dict["genericworkflow.GenericTask.shell_command"] = subsampling_command
new_inputs_dict["genericworkflow.GenericTask.input_files"] =  [target_bam, target_bam_index, ref_sequence, ref_sequence_index]
new_inputs_dict["genericworkflow.docker_override"] = "gcr.io/bioskryb/bwa_and_samtools:latest"
new_inputs_dict["genericworkflow.GenericTask.machine_mem_gb"] = 6
new_inputs_dict["genericworkflow.GenericTask.disk_space_gb"] = 100
wdl_inputs_iobytes = io.BytesIO(json.dumps(new_inputs_dict).encode())
temp_resp = cwt.submit(auth = auth_obj, wdl_file=wdl_iobytes, inputs_files=wdl_inputs_iobytes, options_file = io.BytesIO(json.dumps(options_dict).encode()))
resp_list.append(['NA12878_Bulk2', temp_resp.json()['id']])

pd.DataFrame([[_a, _b, cwt.abort(_b, auth_obj).json()['status']] for _a, _b in resp_list])

temp_resp.content

inputs_dict

options_dict = {
    "google_labels": 
        {
                "run-name": "generic-wdl-practice"
        },
    "default_runtime_attributes": {
        "noAddress": False
    },
    "jes_gcs_root": "gs://bioskryb-dev-cromwell-runs-8dhf5d",
    "workflow_failure_mode": "ContinueWhilePossible",
    "use_relative_output_paths": False,
    "final_workflow_outputs_dir": "gs://scratch_space/delme/"
}

temp_resp = cwt.submit(auth = auth_obj, wdl_file=wdl_iobytes, inputs_files=wdl_inputs_iobytes, options_file = io.BytesIO(json.dumps(options_dict).encode()))

temp_resp.content

mdata = cwt.metadata(temp_resp.json()['id'], auth_obj)

# # **Non-subsampled Testing**

# ## Testing gCVF generation

# ### Joint Genotyping GATK gVCFs using GATK's joint genotyper

# Here we perform joint genotyping on 16 samples from lots 16 - 20. The samples are not subsampled,but are all similar in file size.
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

# ### Joint Genotyping DRAGEN gVCFs using GATK's joint genotyper

# We run the following from the commandline int the grab-illumina-data bucket to acquire the DRAGEN results for the above samples. We only acquire the gVCFs and indexes intitially:

# `gsutil cp $(find -L  /home/rob_carter/BaseSpace/Projects/4249/AppResults/*/Files/ -name "*gvcf.gz*"  |xargs -i readlink -f {} | grep -P "(4249-JW-43)|(4249-JW-38)|(4249-JW-18)|(4249-JW-19)|(4249-JW-40)|(4249-JW-36)|(4249-JW-42)|(4249-JW-16)|(4249-JW-17)|(4249-JW-23)|(4249-JW-10)|(4249-JW-27)|(4249-JW-32)|(4249-JW-25)|(4249-JW-30)|(4249-JW-14)") gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/03-joint_genotyping_input_gvcfs/dragen/`
#

dragen_gvcf_filenames_n16 = !gsutil ls gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/03-joint_genotyping_input_gvcfs/dragen/**vcf.gz
dragen_gvcf_filenames_n16 = dragen_gvcf_filenames_n16 + ["gs://bioskryb-cromwell-outputs-jd7fh2/temp-folder-for-deletion/NA12878_Bulk2_merged_n450x10e6.g.vcf.gz"]
dragen_gvcf_filenames_n16

df_list = [[re.sub("(_combined.+)|(_merged.+)", "", os.path.basename(_blobname)), _blobname] for _blobname in dragen_gvcf_filenames_n16]
pd.DataFrame(df_list).to_csv('tracked_data/dragen_n16_nonsubsampled_sample_map.tsv', sep = "\t", header = False, index= False)
# !gsutil cp tracked_data/dragen_n16_nonsubsampled_sample_map.tsv gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/04-nosubsampling_n16_joint_genotyping/

jj_dragen_input_json = cromwell_functions.download_and_read_inputs_json("gs://bioskryb_dev_wdl_and_inputs/gatk-workflows/gatk4-germline-snps-indels/2.0.0/JointGenotyping.hg38.wgs.inputs.json", storage_client)

jj_dragen_input_json['JointGenotyping.sample_name_map'] = "gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/04-nosubsampling_n16_joint_genotyping/dragen_n16_nonsubsampled_sample_map.tsv"
jj_dragen_input_json['JointGenotyping.callset_name'] = "nosubsampling_n16_dragen_joint_genotyping_gatk"

# +
project_name = "ComparingGatkSentieonDRAGEN"
final_output_bucket = "gs://bioskryb-work-d8f6s9"
cromwell_runs_bucket = "gs://bioskryb-dev-cromwell-runs-8dhf5d"

    
output_base_location = "{}/{}/{}".format(final_output_bucket, project_name, "04-nosubsampling_n16_joint_genotyping")
jj_dragen_options_dict = {
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

input_dragen_iobytes = io.BytesIO(json.dumps(jj_dragen_input_json).encode())
options_dragen_iobytes = io.BytesIO(json.dumps(jj_dragen_options_dict).encode())
jj_resp = cwt.submit(auth_obj, 
       wdl_file = cromwell_functions.get_wdl_iobytes("gs://bioskryb_dev_wdl_and_inputs/gatk-workflows/gatk4-germline-snps-indels/2.0.0/JointGenotyping.wdl", storage_client), 
       inputs_files = input_dragen_iobytes,
       options_file = options_dragen_iobytes
)

jj_resp.content

# ### Generating Sentieon gVCFs and joint genotyping using GATK

# First, generate gVCFs

r1_fastq_files = !gsutil ls gs://bioskryb-vumc-data/Nova181_H2VMKDSXY/*_R1_*fastq.gz
r2_fastq_files = !gsutil ls gs://bioskryb-vumc-data/Nova182_H2VGNDSXY/*_R1_*fastq.gz
all_r1_fastqs = r1_fastq_files + r2_fastq_files

all_r1_fastqs
valid_fastq_r1_filenames_n16 = [_file for _file in all_r1_fastqs if re.sub("_S[01-9].*.fastq.gz", "", os.path.basename(_file)) in non_subsampled_sample_names]
valid_fastq_r1_filenames_n16 = sorted(valid_fastq_r1_filenames_n16)

sample_to_fastq_dict = {}
for _file in valid_fastq_r1_filenames_n16:
    sample_name = re.sub("_S[01-9].*.fastq.gz", "", os.path.basename(_file))
    try:
        sample_to_fastq_dict[sample_name].append(_file)
    except:
        sample_to_fastq_dict[sample_name] = [_file]

# + jupyter={"outputs_hidden": true}
sample_to_fastq_dict
# -

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

# !ls tracked_data/*.n16.sentieoninput.json | parallel -j 16 --max-args 1 "tracked_data/time_and_run_sentieon.sh {}"

# #### Then, run joint genotyping using GATK's haplotypecaller

all_gvcf_files = !gsutil ls 


gvcf_sentieon_filenames = !gsutil ls gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/03-joint_genotyping_input_gvcfs/sentieon**.g.vcf.gz
gvcf_sentieon_filenames_n16 = gvcf_sentieon_filenames + ["gs://bioskryb-cromwell-outputs-jd7fh2/temp-folder-for-deletion/NA12878_Bulk2_merged_n450x10e6.g.vcf.gz"]
gvcf_sentieon_filenames_n16

df_list = [] #[[re.sub("(.g.vcf.gz)|(_merged.+)", "", os.path.basename(_blobname)), _blobname] for _blobname in gvcf_filenames_n16]
for _filename in gvcf_sentieon_filenames_n16:
    if re.search("4249", _filename):
        sample_name = _filename.split("/")[-3]
        df_list.append([sample_name, _filename])
    else:
        sample_name = re.sub("_merged.+", "", os.path.basename(_filename))
        df_list.append([sample_name, _filename])
pd.DataFrame(df_list).to_csv('tracked_data/sentieon_n16_nonsubsampled_sample_map.tsv', sep = "\t", header = False, index= False)
# !gsutil cp tracked_data/sentieon_n16_nonsubsampled_sample_map.tsv gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/04-nosubsampling_n16_joint_genotyping/

jj_sentieon_input_json = cromwell_functions.download_and_read_inputs_json("gs://bioskryb_dev_wdl_and_inputs/gatk-workflows/gatk4-germline-snps-indels/2.0.0/JointGenotyping.hg38.wgs.inputs.json", storage_client)

jj_sentieon_input_json['JointGenotyping.sample_name_map'] = "gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/04-nosubsampling_n16_joint_genotyping/sentieon_n16_nonsubsampled_sample_map.tsv"
jj_sentieon_input_json['JointGenotyping.callset_name'] = "nosubsampling_n16_sentieon_joint_genotyping_gatk"

# +
project_name = "ComparingGatkSentieonDRAGEN"
final_output_bucket = "gs://bioskryb-work-d8f6s9"
cromwell_runs_bucket = "gs://bioskryb-dev-cromwell-runs-8dhf5d"

    
output_base_location = "{}/{}/{}".format(final_output_bucket, project_name, "data/04-nosubsampling_n16_joint_genotyping")
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
                "pipeline-name": "gatk4-germline-snps-indels",
                "project-name": "comparing-gatk-sentieon-dragen"
        }
}
# -

input_sentieon_iobytes = io.BytesIO(json.dumps(jj_sentieon_input_json).encode())
options_sentieon_iobytes = io.BytesIO(json.dumps(jj_sentieon_options_dict).encode())
jj_sentieon_resp = cwt.submit(auth_obj, 
       wdl_file = cromwell_functions.get_wdl_iobytes("gs://bioskryb_dev_wdl_and_inputs/gatk-workflows/gatk4-germline-snps-indels/2.0.0/JointGenotyping.wdl", storage_client), 
       inputs_files = input_sentieon_iobytes,
       options_file = options_sentieon_iobytes
)

jj_sentieon_resp.content

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


