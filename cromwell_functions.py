import re, json, os, logging, io, requests

from copy import deepcopy

from cromwell_tools.cromwell_api import CromwellAPI as cwt
import pandas as pd



def download_and_read_inputs_json(wdl_inputs_gcs_filename, storage_client):
    bucket_name = re.sub("^gs://", "", wdl_inputs_gcs_filename).split("/")[0]
    blob_name = "/".join(re.sub("^gs://", "", wdl_inputs_gcs_filename).split("/")[1:])
    wdl_inputs_bucket = storage_client.bucket(bucket_name)
    wdl_inputs_blob = wdl_inputs_bucket.blob(blob_name)
    local_wdl_inputs_filename = os.path.join("/tmp", os.path.basename(blob_name))
    print("Saving file to {}".format(local_wdl_inputs_filename))
    wdl_inputs_blob.download_to_filename(local_wdl_inputs_filename)
    with open(local_wdl_inputs_filename, 'r') as ifh:
        wdl_inputs_dict = json.load(ifh)
    #os.remove(local_wdl_inputs_filename)
    return(wdl_inputs_dict)

def get_wdl_iobytes(wdl_gcs_filename, storage_client):
    """
    """
    bucket_name = re.sub("^gs://", "", wdl_gcs_filename).split("/")[0]
    blob_name = "/".join(re.sub("^gs://", "", wdl_gcs_filename).split("/")[1:])
    wdl_bucket = storage_client.bucket(bucket_name)
    wdl_blob = wdl_bucket.blob(blob_name)
    #local_wdl_filename = os.path.join("/tmp", os.path.basename(blob_name))

    #logger.debug("Acquiring blob: {} ".format(blob_name))
    wdl_string = wdl_blob.download_as_string()
    #logger.debug("Acquired wdl blob as type: {} ".format(type(wdl_string)))

    return io.BytesIO(wdl_string)

def get_blob_from_gcp_location(gcp_location, storage_client):
    split_gcp_list = gcp_location.split("/")
    bucket_name = split_gcp_list[2]
    blob_name =  os.path.join(*split_gcp_list[3:])
    bucket_obj = storage_client.bucket(bucket_name)
    blob_obj = bucket_obj.get_blob(blob_name)
    return blob_obj

target_chromosomes = ["chr{}".format(str(_int+1)) for _int in range(22)] + ['chrX', 'chrY', 'chrM']

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