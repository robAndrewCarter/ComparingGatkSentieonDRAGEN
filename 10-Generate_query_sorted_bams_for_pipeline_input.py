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

import vcf

from cromwell_tools.cromwell_api import CromwellAPI as cwt
from cromwell_tools import cromwell_auth
from google.cloud import storage

import cromwell_functions
from cromwell_functions import defaultCromwellRunner
# -

auth_obj = cromwell_auth.CromwellAuth(url = "http://34.69.50.44:8000", header= None, auth = None)
auth_obj

os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = "/home/rcarter/.google/bioskryb-81ce35d92471.json"

storage_client = storage.Client('bioskryb')

# ## GATK - with Ref

# Here we perform joint genotyping on 16 samples from lots 16 - 20. The samples are not subsampled, but are all similar in file size.
# We first run the GATK gVCFs through GATK's haplotypeCaller, but including the subsampled Bulk2 reference sequence.

non_subsampled_sample_names = list(as_and_wgs_and_dups_metrics_df[as_and_wgs_and_dups_metrics_df["sample_prop"]< 1.06].sort_values('sample_prop', ascending = False).head(16).loc[:, ['sample_name', "sample_prop"]].head(24)['sample_name'])

all_gvcf_files = !gsutil ls gs://bioskryb-cromwell-outputs-jd7fh2/VUMC_Nova183_H37TLDSXY_DUAL_PTA_4249-JW-49-64/**.g.vcf.gz
all_gvcf_files[:2]


gvcf_filenames = !gsutil ls gs://bioskryb-vumc-data/Full_PTA_exome_deepdive_with_other_techs/data/01-gatk_cromwell_output/**.g.vcf.gz
gvcf_filenames_n16 = [_file for _file in gvcf_filenames if re.sub(".g.vcf.gz", "", os.path.basename(_file)) in non_subsampled_sample_names] + ["gs://bioskryb-cromwell-outputs-jd7fh2/temp-folder-for-deletion/NA12878_Bulk2_merged_n450x10e6.g.vcf.gz"]

df_list = [[re.sub("(.g.vcf.gz)|(_merged.+)", "", os.path.basename(_blobname)), _blobname] for _blobname in gvcf_filenames_n16_noref]
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

non_subsampled_sample_names = ["4249-JW-10","4249-JW-14","4249-JW-16","4249-JW-17","4249-JW-18","4249-JW-19","4249-JW-23","4249-JW-25","4249-JW-27","4249-JW-30","4249-JW-32","4249-JW-36","4249-JW-38","4249-JW-40","4249-JW-42","4249-JW-43"]

all_gvcf_files = !gsutil ls gs://bioskryb-cromwell-outputs-jd7fh2/VUMC_Nova183_H37TLDSXY_DUAL_PTA_4249-JW-49-64/**.g.vcf.gz
all_gvcf_files[:2]


gvcf_filenames = !gsutil ls gs://bioskryb-vumc-data/Full_PTA_exome_deepdive_with_other_techs/data/01-gatk_cromwell_output/**.g.vcf.gz
gvcf_filenames_n16_no_ref = [_file for _file in gvcf_filenames if re.sub(".g.vcf.gz", "", os.path.basename(_file)) in non_subsampled_sample_names] 

df_list = [[re.sub("(.g.vcf.gz)|(_merged.+)", "", os.path.basename(_blobname)), _blobname] for _blobname in gvcf_filenames_n16_no_ref]
pd.DataFrame(df_list).to_csv('tracked_data/gatk_n16_no_ref_nonsubsampled_sample_map.tsv', sep = "\t", header = False, index= False)
# !gsutil cp tracked_data/gatk_n16_no_ref_nonsubsampled_sample_map.tsv gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/04-nosubsampling_n16_joint_genotyping/

jj_input_json = cromwell_functions.download_and_read_inputs_json("gs://bioskryb_dev_wdl_and_inputs/gatk-workflows/gatk4-germline-snps-indels/2.0.0/JointGenotyping.hg38.wgs.inputs.json", storage_client)

jj_input_json['JointGenotyping.sample_name_map'] = "gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/04-nosubsampling_n16_joint_genotyping/gatk_n16_no_ref_nonsubsampled_sample_map.tsv"
jj_input_json['JointGenotyping.callset_name'] = "nosubsampling_n16_no_ref_gatk_joint_genotyping_gatk"

# +
project_name = "ComparingGatkSentieonDRAGEN"
final_output_bucket = "gs://bioskryb-work-d8f6s9"
cromwell_runs_bucket = "gs://bioskryb-dev-cromwell-runs-8dhf5d"

    
output_base_location = "{}/{}/{}".format(final_output_bucket, project_name, "data/04-nosubsampling_n16_no_ref_joint_genotyping")
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

# ## Sentieon - gVCF generation

# The same 16 samples from above (plus the reference) are joint genotyped using Sentieon

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

# ### Sentieon - Joint-genotyping

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
gvcftyper_command_template = "sentieon driver -r {ref} -t 32  --algo GVCFtyper -v " + variant_subcommand + " -d {dbsnp} {out_vcf}"
gvcftyper_command_template
# -

vcf_df

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
inputs_dict["genericworkflow.GenericTask.machine_mem_gb"] = 32
inputs_dict["genericworkflow.GenericTask.disk_space_gb"] = 200
inputs_dict["genericworkflow.GenericTask.cpu"] = 32

gctyper_dcr = defaultCromwellRunner(auth_obj, gvcftyper_command_template, inputs_df, jj_sentieon_options_dict, inputs_dict, wdl_text= requests.get("https://raw.githubusercontent.com/gatk-workflows/generic-wdl-script/master/generic-workflow.wdl").text)
# -

gctyper_dcr.submit_jobs()

gctyper_dcr.get_status()

# ### Sentieion - VarCal

# We next run Sentieon's SNP variant recalibration

# +
varcal_df = pd.DataFrame([
    ["gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/04-nosubsampling_n16_joint_genotyping/sentieon/f7791408-5ccb-415d-94c0-72c1802d5104/call-GenericTask/glob-eb940d6b5fb9e12ecfa084f32facbc84/nosubsampling_n16_sentieon_joint_genotyping_sentieon.vcf.gz", 
     "gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/04-nosubsampling_n16_joint_genotyping/sentieon/f7791408-5ccb-415d-94c0-72c1802d5104/call-GenericTask/glob-eb940d6b5fb9e12ecfa084f32facbc84/nosubsampling_n16_sentieon_joint_genotyping_sentieon.vcf.gz.tbi",
     "gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.fasta",
     "gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.fasta.fai",
     "gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz",
     "gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi",
     "gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz",
     "gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi",
     "gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
     "gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi",
     "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf",
     "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx",
     "snp_tranches_out.txt",
     "nosubsampling_n16_sentieon_joint_genotyping_sentieon.snp.recal.vcf.gz"
    ]
], columns= ['vcf=input', 'vcf_index=input', 'ref=input', 'ref_index=input', "hapmap=input", "hapmap_index=input", "omni=input", "omni_index=input", "onekg=input", "onekg_index=input","dbsnp=input","dbsnp_index=input", "snp_tranches_out=output", "snp_recal_vcf=output"])

varcal_command_template = "sentieon driver -r {ref} -t 6  --algo VarCal --annotation QD --annotation MQRankSum --annotation ReadPosRankSum --annotation FS --annotation MQ --annotation SOR --annotation DP -v {vcf} --var_type SNP --tranches_file {snp_tranches_out} --resource {hapmap} --resource_param hapmap,known=false,training=true,truth=true,prior=15  --resource {omni} --resource_param omni,known=false,training=true,truth=true,prior=12 --resource {onekg} --resource_param onekg,known=false,training=true,truth=false,prior=10 --resource {dbsnp} --resource_param dbsnp,known=true,training=false,truth=false,prior=7 {snp_recal_vcf}"

# +
project_name = "ComparingGatkSentieonDRAGEN"
final_output_bucket = "gs://bioskryb-work-d8f6s9"
cromwell_runs_bucket = "gs://bioskryb-dev-cromwell-runs-8dhf5d"

    
output_base_location = "{}/{}/{}".format(final_output_bucket, project_name, "data/04-nosubsampling_n16_joint_genotyping/sentieon/varcal")
jj_varcal_sentieon_options_dict = {
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
inputs_dict = requests.get("https://raw.githubusercontent.com/gatk-workflows/generic-wdl-script/master/generic-workflow.input.json").json()
inputs_dict["genericworkflow.docker_override"] = "gcr.io/bioskryb/sentieon-201911:latest"
inputs_dict["genericworkflow.GenericTask.machine_mem_gb"] = 6
inputs_dict["genericworkflow.GenericTask.disk_space_gb"] = 150
inputs_dict["genericworkflow.GenericTask.cpu"] = 6

snp_varcal_dcr = defaultCromwellRunner(auth_obj, varcal_command_template, varcal_df, jj_varcal_sentieon_options_dict, inputs_dict, wdl_text= requests.get("https://raw.githubusercontent.com/gatk-workflows/generic-wdl-script/master/generic-workflow.wdl").text)
# -

snp_varcal_dcr.submit_jobs()

snp_varcal_dcr.get_status()



# We next run Sentieon's indel variant recalibration

# +
indel_varcal_df = pd.DataFrame([
    ["gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/04-nosubsampling_n16_joint_genotyping/sentieon/f7791408-5ccb-415d-94c0-72c1802d5104/call-GenericTask/glob-eb940d6b5fb9e12ecfa084f32facbc84/nosubsampling_n16_sentieon_joint_genotyping_sentieon.vcf.gz", 
     "gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/04-nosubsampling_n16_joint_genotyping/sentieon/f7791408-5ccb-415d-94c0-72c1802d5104/call-GenericTask/glob-eb940d6b5fb9e12ecfa084f32facbc84/nosubsampling_n16_sentieon_joint_genotyping_sentieon.vcf.gz.tbi",
     "gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.fasta",
     "gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.fasta.fai",
     "gs://gcp-public-data--broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz", 
     "gs://gcp-public-data--broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi", 
     "gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
     "gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi",
     "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf",
     "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx",
     "indel_tranches_out.txt",
     "nosubsampling_n16_sentieon_joint_genotyping_sentieon.indel.recal.vcf.gz"
    ]
], columns= ['vcf=input', 'vcf_index=input', 'ref=input', 'ref_index=input', "axiom=input", "axiom_index=input", "mills=input", "mills_index=input", "dbsnp=input","dbsnp_index=input", "indel_tranches_out=output", "indel_recal_vcf=output"])

indel_varcal_command_template = "sentieon driver  -r {ref} -t 6  --algo VarCal --annotation FS --annotation ReadPosRankSum --annotation MQRankSum --annotation QD --annotation SOR --annotation DP -v {vcf} --var_type INDEL --tranches_file {indel_tranches_out} --resource {mills} --resource_param  mills,known=false,training=true,truth=true,prior=12 --resource {axiom} --resource_param axiom,known=false,training=true,truth=false,prior=10 --resource {dbsnp} --resource_param dbsnp,known=true,training=false,truth=false,prior=2 {indel_recal_vcf}"
indel_varcal_command_template

# +
project_name = "ComparingGatkSentieonDRAGEN"
final_output_bucket = "gs://bioskryb-work-d8f6s9"
cromwell_runs_bucket = "gs://bioskryb-dev-cromwell-runs-8dhf5d"

    
output_base_location = "{}/{}/{}".format(final_output_bucket, project_name, "data/04-nosubsampling_n16_joint_genotyping/sentieon/varcal")
jj_varcal_sentieon_options_dict = {
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
inputs_dict = requests.get("https://raw.githubusercontent.com/gatk-workflows/generic-wdl-script/master/generic-workflow.input.json").json()
inputs_dict["genericworkflow.docker_override"] = "gcr.io/bioskryb/sentieon-201911:latest"
inputs_dict["genericworkflow.GenericTask.machine_mem_gb"] = 6
inputs_dict["genericworkflow.GenericTask.disk_space_gb"] = 150
inputs_dict["genericworkflow.GenericTask.cpu"] = 6

indel_varcal_dcr = defaultCromwellRunner(auth_obj, indel_varcal_command_template, indel_varcal_df, jj_varcal_sentieon_options_dict, inputs_dict, wdl_text= requests.get("https://raw.githubusercontent.com/gatk-workflows/generic-wdl-script/master/generic-workflow.wdl").text)
# -

indel_varcal_dcr.submit_jobs()

indel_varcal_dcr.get_status()

# ### Sentieon - ApplyVarCal

# +
apply_varcal_df = pd.DataFrame([
    ["gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/04-nosubsampling_n16_joint_genotyping/sentieon/f7791408-5ccb-415d-94c0-72c1802d5104/call-GenericTask/glob-eb940d6b5fb9e12ecfa084f32facbc84/nosubsampling_n16_sentieon_joint_genotyping_sentieon.vcf.gz", 
     "gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/04-nosubsampling_n16_joint_genotyping/sentieon/f7791408-5ccb-415d-94c0-72c1802d5104/call-GenericTask/glob-eb940d6b5fb9e12ecfa084f32facbc84/nosubsampling_n16_sentieon_joint_genotyping_sentieon.vcf.gz.tbi",
     "gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.fasta",
     "gs://bioskryb-dev-resources-4j2g6d/Homo_sapiens_assembly38.fasta.fai",
     "gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/04-nosubsampling_n16_joint_genotyping/sentieon/varcal/7af01187-7275-45a7-ab26-9b923fc8d714/call-GenericTask/glob-eb940d6b5fb9e12ecfa084f32facbc84/indel_tranches_out.txt",
     "gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/04-nosubsampling_n16_joint_genotyping/sentieon/varcal/7af01187-7275-45a7-ab26-9b923fc8d714/call-GenericTask/glob-eb940d6b5fb9e12ecfa084f32facbc84/nosubsampling_n16_sentieon_joint_genotyping_sentieon.indel.recal.vcf.gz",
     "gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/04-nosubsampling_n16_joint_genotyping/sentieon/varcal/7af01187-7275-45a7-ab26-9b923fc8d714/call-GenericTask/glob-eb940d6b5fb9e12ecfa084f32facbc84/nosubsampling_n16_sentieon_joint_genotyping_sentieon.indel.recal.vcf.gz.tbi",
     "gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/04-nosubsampling_n16_joint_genotyping/sentieon/varcal/e8702941-94ce-430e-98b6-6cb5e1d2b25d/call-GenericTask/glob-eb940d6b5fb9e12ecfa084f32facbc84/snp_tranches_out.txt",
     "gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/04-nosubsampling_n16_joint_genotyping/sentieon/varcal/e8702941-94ce-430e-98b6-6cb5e1d2b25d/call-GenericTask/glob-eb940d6b5fb9e12ecfa084f32facbc84/nosubsampling_n16_sentieon_joint_genotyping_sentieon.snp.recal.vcf.gz",
     "gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/04-nosubsampling_n16_joint_genotyping/sentieon/varcal/e8702941-94ce-430e-98b6-6cb5e1d2b25d/call-GenericTask/glob-eb940d6b5fb9e12ecfa084f32facbc84/nosubsampling_n16_sentieon_joint_genotyping_sentieon.snp.recal.vcf.gz.tbi",
     "nosubsampling_n16_sentieon_joint_genotyping_sentieon.indel.varcal.final.vcf.gz",
     "nosubsampling_n16_sentieon_joint_genotyping_sentieon.indel.varcal.final.vcf.gz.tbi"
    ]
], columns= ['vcf=input', 'vcf_index=input', 'ref=input', 'ref_index=input', "indel_tranches=input", "indel_recal_vcf=input", "indel_recal_index_vcf=input", "snp_tranches=input", "snp_recal_vcf=input", "snp_recal_index_vcf=input", "final_varcal_vcf=output", "final_varcal_vcf_index=output"])

apply_varcal_command_template = "sentieon driver  -r {ref} -t 6  --algo ApplyVarCal --vqsr_model var_type=SNP,recal={snp_recal_vcf},tranches_file={snp_tranches},sensitivity=99.7 --vqsr_model var_type=INDEL,recal={indel_recal_vcf},tranches_file={indel_tranches},sensitivity=99.0 --vcf {vcf} {final_varcal_vcf}"
apply_varcal_command_template

# +
project_name = "ComparingGatkSentieonDRAGEN"
final_output_bucket = "gs://bioskryb-work-d8f6s9"
cromwell_runs_bucket = "gs://bioskryb-dev-cromwell-runs-8dhf5d"

    
output_base_location = "{}/{}/{}".format(final_output_bucket, project_name, "data/04-nosubsampling_n16_joint_genotyping/sentieon/applyvarcal")
jj_applyvarcal_sentieon_options_dict = {
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
inputs_dict = requests.get("https://raw.githubusercontent.com/gatk-workflows/generic-wdl-script/master/generic-workflow.input.json").json()
inputs_dict["genericworkflow.docker_override"] = "gcr.io/bioskryb/sentieon-201911:latest"
inputs_dict["genericworkflow.GenericTask.machine_mem_gb"] = 6
inputs_dict["genericworkflow.GenericTask.disk_space_gb"] = 150
inputs_dict["genericworkflow.GenericTask.cpu"] = 6

apply_varcal_dcr = defaultCromwellRunner(auth_obj, apply_varcal_command_template, apply_varcal_df, jj_applyvarcal_sentieon_options_dict, inputs_dict, wdl_text= requests.get("https://raw.githubusercontent.com/gatk-workflows/generic-wdl-script/master/generic-workflow.wdl").text)
# -

apply_varcal_dcr.submit_jobs()

apply_varcal_dcr.get_status()

# ## Copy vcfs into scratch for analysis by R

# The following was performed to acuire the vcf files for Sentieon and gatk samples:

# `rcarter@bioskryb-dt:scratch$ gsutil cp gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/04-nosubsampling_n16_joint_genotyping/sentieon/applyvarcal/95f90949-445b-4493-b80e-ddaf09353cda/call-GenericTask/**.vcf.gz* .`
#
# `rcarter@bioskryb-dt:scratch$ gsutil cp gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/04-nosubsampling_n16_joint_genotyping/JointGenotyping/e00ab3ad-4e37-4676-af07-8a274f08c741/call-FinalGatherVcf/attempt-2/**.vcf.gz* .`
#
# `rcarter@bioskryb-dt:scratch$ cp  ~/BaseSpace/Projects/.id.161953803/AppSessions/.id.231671923/AppResults.239697619.joint/Files/dragen-joint.hard-filtered.vcf.gz* .`

# GATK 4.1.3.0 (in Docker container) was used to isolate the SNP variants from each of the three platforms

# `for i in *.vcf.gz; do PREFIX=$(echo $i | sed 's/.vcf.gz//'); gatk SelectVariants --select-type SNP  -V ${PREFIX}.vcf.gz -O ${PREFIX}.snps_only.vcf.gz; done`

# The vcfs were split into indels and SNP vcfs, but were not filtered. Note that the Sentieon VCF was incorrenctly named, so it was renamed first to remove the '.indel' substring.

# `for i in nosubsampling_n16_sentieon_joint_genotyping_sentieon.indel.*; do NEW_NAME=$(echo $i | sed 's/.indel//'); mv $i $NEW_NAME; done`

# %%writefile tracked_data/run_rtg.sh
#!/usr/bin/env bash
# echo $2 | sed 's/,/\n/g' | parallel --max-args=1 --jobs 8 /home/rcarter/rtg-tools/bin/rtg-tools-3.10.1-4d58eadb/rtg RTG_MEM=6G vcfeval -b $1 -c $1 --sample=NA12878_Bulk2,{} -o scratch/NA12878_Bulk2_vs_{}-${3} -t /home/rcarter/R_bioskryb_na12878_analysis/NA12878_comparison/SDF/ --squash-ploidy

# !chmod 755 tracked_data/run_rtg.sh

# all non-bulk samples are compared against bulk  for each of the three genotyping platforms

# #### GATK

gatk_snps_vcf = vcf.Reader(open('scratch/nosubsampling_n16_gatk_joint_genotyping_gatk.vcf.gz', 'rb'))
gatk_samplenames = gatk_snps_vcf.samples
gatk_samplenames

gatk_joined_samplenames = ",".join(gatk_samplenames)
# !tracked_data/run_rtg.sh scratch/nosubsampling_n16_gatk_joint_genotyping_gatk.snps_only.vcf.gz $gatk_joined_samplenames n16_gatk

# #### Sentieon

sentieon_snps_vcf = vcf.Reader(open('scratch/nosubsampling_n16_sentieon_joint_genotyping_sentieon.varcal.final.snps_only.vcf.gz', 'rb'))
sentieon_samplenames = [_sample_name for _sample_name in sentieon_snps_vcf.samples if _sample_name != 'NA12878_Bulk2']
sentieon_samplenames

sentieon_joined_samplenames = ",".join(sentieon_samplenames)
# !tracked_data/run_rtg.sh scratch/nosubsampling_n16_sentieon_joint_genotyping_sentieon.varcal.final.snps_only.vcf.gz $sentieon_joined_samplenames n16_sentieon

# #### DRAGEN

# %%writefile tracked_data/run_dragen_rtg.sh
#!/usr/bin/env bash
# echo $2 | sed 's/,/\n/g' | parallel --max-args=1 --jobs 8 /home/rcarter/rtg-tools/bin/rtg-tools-3.10.1-4d58eadb/rtg RTG_MEM=6G vcfeval -b $1 -c $1 --sample=NA12878-Bulk2-merged-n450x10e6,{} -o scratch/NA12878-Bulk2-merged-n450x10e6_vs_{}-${3} -t /home/rcarter/R_bioskryb_na12878_analysis/NA12878_comparison/SDF/ --squash-ploidy

# !chmod 755 tracked_data/run_dragen_rtg.sh

dragen_snps_vcf = vcf.Reader(open('scratch/dragen-joint.hard-filtered.snps_only.vcf.gz', 'rb'))
dragen_samplenames = dragen_snps_vcf.samples
dragen_samplenames = [_sample_name for _sample_name in dragen_snps_vcf.samples if _sample_name != 'NA12878-Bulk2-merged-n450x10e6']
dragen_samplenames

dragen_joined_samplenames = ",".join(dragen_samplenames)
# !tracked_data/run_dragen_rtg.sh scratch/dragen-joint.hard-filtered.snps_only.vcf.gz $dragen_joined_samplenames n16_dragen

# ### Gather RTG summary files

summary_files = ls 

df_list = []
for _summary_file in summary_files:
    sample_name = re.sub("^.+vs_(.+)-idt.+", "\\1", _summary_file.split("/")[1])
    temp_df = pd.read_csv(_summary_file, skiprows=3, sep = "[\t ]+", header = None, names=["Threshold","True-pos-baseline","True-pos-call","False-pos","False-neg","Precision","Sensitivity","F-measure"])
    temp_df['sample_name'] = sample_name
    if re.search("n6e7", _summary_file):
        temp_df['subsample_level'] = 60000000
    elif re.search("ne7", _summary_file):
        temp_df['subsample_level'] = 20000000
    else:
        temp_df['subsample_level'] = 40000000
    df_list.append(temp_df)
rtg_summary_df = pd.concat(df_list)

rtg_summary_df.to_csv("data_tracked/rtg_summary_all_subsamples.tsv", sep = "\t")
