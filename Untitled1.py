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
import re, json, os, logging, io, pprint, subprocess, requests, sys
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

sample_input_files_dict = {"in2" : "gs://mojo_magoo/foo.wdl", "in1" : "gs://mojo_magoo/bar.wdl"}

internal_input_files_dict = deepcopy(sample_input_files_dict)



for _key,_value in internal_input_files_dict.items():
    internal_input_files_dict[_key] = os.path.basename(_value)
internal_input_files_dict


class defaultCromwellRunner():
    def __init__(self, auth_obj, command_template, inputs_df, options_dict, inputs_dict, wdl_text, label = "generic", move_all_outputs = True):
        self.command_template = command_template
        self.inputs_df = inputs_df
        self.move_all_outputs = move_all_outputs
        
        self.options_dict = options_dict
        self.inputs_dict = inputs_dict
        self.wdl_text = wdl_text
        self.resp_list = []
        
        (self.local_inputs_df, self.remote_files_df, self.coltypes_dict) = self.create_and_return_local_and_remote_inputs_and_coltypes_dict(inputs_df)
        self.commands_dict = self.create_run_commands_dict(self.local_inputs_df)
        self.remote_files_dict = self.remote_files_df.apply(lambda x: list(x), axis = 1).to_dict()
        self.final_args_dict = self.get_final_args_dict(self.commands_dict, self.remote_files_dict)
        
        if auth_obj:
            self.auth_obj = auth_obj
        else:
            self.auth_obj = None
            
    def authorize(self, auth_obj):
        self.auth_obj = auth_obj
        
    def get_final_args_dict(self, commands_dict, remote_files_dict):
        assert(commands_dict.keys() == remote_files_dict.keys())
        final_args_dict = {}
        for _key in commands_dict.keys():
            final_args_dict[_key] = {'command': commands_dict[_key], 'remote_input_files':  remote_files_dict[_key]}
        return final_args_dict
        
    def create_and_return_local_and_remote_inputs_and_coltypes_dict(self, inputs_df):
        coltypes_dict = self.create_coltypes_dict(inputs_df)
        remote_inputs_df = inputs_df.loc[:, coltypes_dict['remote_input']] 
        local_inputs_df = inputs_df.apply(lambda x: x.apply(os.path.basename)).rename(columns = lambda x: x.split("=")[0])
        return(local_inputs_df, remote_inputs_df, coltypes_dict)
    
    def create_run_commands_dict(self, local_inputs_df):
        commands_dict = {}
        temp_df_records_dict_dict = local_inputs_df.to_dict('index')
        for _index, _sample_args_dict in temp_df_records_dict_dict.items():
            temp_command = self.command_template.format(**_sample_args_dict)
            if self.move_all_outputs and len(self.coltypes_dict['output']) > 0:
                command_suffix = " && mv {} ./output/".format(" ".join([_sample_args_dict[_colname] for _colname in self.coltypes_dict['output']]))
                temp_command  += command_suffix
            commands_dict[_index] = temp_command
        return commands_dict
            
    def create_coltypes_dict(self, inputs_df):
        coltypes_dict = {'input': [], 'output': [], "remote_input": []}
        for _colname in inputs_df.columns:
            colname_split = _colname.split("=")
            if len(colname_split) != 2:
                sys.exit("invalid colname: {}".format(_colname))
            io_type = colname_split[1]
            local_variable_name = colname_split[0]
            try:
                coltypes_dict[io_type].append(local_variable_name)
            except:
                sys.exit("Invalid iotype: {}".format(io_type))
            if io_type == 'input':
                coltypes_dict['remote_input'].append(_colname)
            
        return(coltypes_dict)
        
        
    def submit_jobs(self):
        resp_list = []
        for _index, _args_dict in self.final_args_dict.items():
            new_inputs_dict = deepcopy(self.inputs_dict)
            new_inputs_dict["genericworkflow.GenericTask.shell_command"] = _args_dict['command']
            new_inputs_dict["genericworkflow.GenericTask.input_files"] = _args_dict['remote_input_files']
            temp_resp = cwt.submit(auth = self.auth_obj, wdl_file=io.BytesIO(self.wdl_text.encode()), inputs_files=io.BytesIO(json.dumps(new_inputs_dict).encode()), options_file=io.BytesIO(json.dumps(self.options_dict).encode()))
            resp_list.append(temp_resp.json()['id'])
            #resp_list.append(new_inputs_dict)
        self.resp_list = resp_list
        return(self.resp_list)
    
    def get_status(self):
        return pd.DataFrame([[_a, cwt.metadata(_a, self.auth_obj).json()['status']] for _a in self.resp_list])
        
    def dry_run(self):
        return(self.final_args_dict)



new_inputs_dict = deepcopy(default_inputs_dict)
#new_inputs_dict["genericworkflow.GenericTask.shell_command"] = sentieon_jj_command
#new_inputs_dict["genericworkflow.GenericTask.input_files"] = ["{}.tbi".format(i) for i in gvcf_updated_sentieon_filenames_n16[:2]] + gvcf_updated_sentieon_filenames_n16[:2]  + [ref_sequence, ref_sequence_index, dbsnp_filename, dbsnp_filename + ".tbi"]
new_inputs_dict["genericworkflow.docker_override"] = "gcr.io/bioskryb/bwa_and_samtools:argh"
new_inputs_dict["genericworkflow.GenericTask.machine_mem_gb"] = 2
new_inputs_dict["genericworkflow.GenericTask.disk_space_gb"] = 8
new_inputs_dict["genericworkflow.GenericTask.cpu"] = 1
new_inputs_dict

output_base_location = 'gs://scratch_space'
cromwell_runs_bucket = "gs://bioskryb-dev-cromwell-runs-8dhf5d"
options_dict = {
        "read_from_cache": False,
        "final_workflow_outputs_dir": output_base_location,
        "default_runtime_attributes": {
            "zones": "us-central1-a us-central1-b us-central1-c us-central1-f",
            "noAddress": False
        },
        "use_relative_output_paths": False,
        "final_call_logs_dir": "{}/call_logs".format(output_base_location),
        "jes_gcs_root": cromwell_runs_bucket,
        "google_labels": {
                "pipeline-name": "generic-workflow",
                "project-name": "generic-project"
        }
}

temp_df = pd.DataFrame([["gs://illumina-fastq-qc-output/1234_ABCD/VEGFA_11_S21_2x10e6/VEGFA_11_S21_2x10e6_md.bam", "line_count_1.txt"], ["gs://illumina-fastq-qc-output/1234_ABCD/VEGFA_11_S21_2x10e6/VEGFA_11_S21_2x10e6_md.bam", "line_count_2.txt"]], columns = ["in1=input", 'out1=output'])
command_template = "samtools view {in1} | wc -l  > {out1}"

# +
default_inputs_dict = requests.get("https://raw.githubusercontent.com/gatk-workflows/generic-wdl-script/master/generic-workflow.input.json").json()

wdl_text = requests.get("https://raw.githubusercontent.com/gatk-workflows/generic-wdl-script/master/generic-workflow.wdl").text
wdl_iobytes = io.BytesIO(wdl_text.encode())
# -

temp_cromwell_runner = defaultCromwellRunner(auth_obj, command_template, temp_df, options_dict, new_inputs_dict, wdl_text)
resp_list = temp_cromwell_runner.submit_jobs()

resp_list[0]

temp_cromwell_runner.get_status()

 #[_var.split("=")[0] for _var in temp_df.columns])


default_inputs_dict

temp_df.apply(lambda x: list(x), axis = 1).to_dict()


