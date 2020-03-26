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

# Here we gather and process known SVs in NA12878 from the svclassify manuscript

# +
import re, os, subprocess, pprint, json
import io
import requests

from google.cloud import storage
import pandas as pd

# +
import copy
from cromwell_tools.cromwell_api import CromwellAPI as cwt
from cromwell_tools import cromwell_auth
from google.cloud import storage

#import cromwell_functions
# -

storage_client = storage.Client('bioskryb')

# !wget -O scratch/Personalis_1000_Genomes_deduplicated_deletions.bed https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/technical/svclassify_Manuscript/Supplementary_Information/Personalis_1000_Genomes_deduplicated_deletions.bed

# The deletion data was liftet over to GRCh38 coordinates, although 10 of the indels failed. These are saved in tracked_data/Personalis_1000_Genomes_deduplicated_deletions.failed_liftover.txt. The successfully lifted-over deletions are found in tracked_data/Personalis_1000_Genomes_deduplicated_deletions.GRCh38.liftover.bed


