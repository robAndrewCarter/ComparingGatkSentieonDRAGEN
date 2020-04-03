import re, json, os, logging, io

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