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
