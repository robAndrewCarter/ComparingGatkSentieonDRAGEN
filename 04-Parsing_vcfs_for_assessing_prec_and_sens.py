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

# Here we acquire and parse vcfs from Sentieon and Dragen for comparison to GATK

# ### Copy Sentieon's DNAscope and DRAGEN's VCFs over to scratch for analysis

dnascoped_vcfs = !gsutil ls gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/01-Sentieon_output/DNAScope/**.vcf.gz*

for _dnascope_vcf in dnascoped_vcfs:
    sample_name = _dnascope_vcf.split("/")[-3]
    output_location = "scratch/{}.{}".format(sample_name, os.path.basename(_dnascope_vcf))
#     !gsutil cp $_dnascope_vcf $output_location

dragen_vcfs = !gsutil ls gs://bioskryb-work-d8f6s9/ComparingGatkSentieonDRAGEN/data/02-DRAGEN_illumina/**.vcf.gz*

for _dragen_vcf in dragen_vcfs:
    sample_name = _dragen_vcf.split("/")[-3]
    output_location = "scratch/{}.{}".format(sample_name, os.path.basename(_dragen_vcf))
#     !gsutil cp $_dragen_vcf $output_location

# ### Split vcfs into SNPs and indels and convert DRAGEN results from GRCh37 to GRCh38

local_sentieon_dnascope_vcfs = !ls scratch/*dnascope.vcf.gz | xargs -i basename {}
local_sentieon_dnascope_vcfs

for _local_dnascope_vcf in local_sentieon_dnascope_vcfs:
    snp_output = re.sub(".vcf.gz", ".snp.vcf.gz", _local_dnascope_vcf)
    indel_output = re.sub(".vcf.gz", ".indel.vcf.gz", _local_dnascope_vcf)
    gatk_command = '"cd /home && gatk SplitVcfs --STRICT false -I {} --INDEL_OUTPUT {} --SNP_OUTPUT {}"'.format(_local_dnascope_vcf, indel_output, snp_output)
#     !docker run -v /home/rcarter/ComparingGatkSentieonDRAGEN/scratch/:/home -ti gcr.io/bioskryb/gatk:4.1.3.0 /bin/bash -c $gatk_command

# Download the chain file for crossmap and run using the py2 environment

# !wget -O scratch/hg19ToHg38.over.chain.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# The following was run on the commandline since it required running from a different environment:

local_dragen_vcfs = !ls scratch/*dragen*.vcf.gz | xargs -i basename {} | grep -v GRCh38
local_dragen_vcfs

for _local_dragen_vcf in local_dragen_vcfs:
    grch38_output = re.sub(".vcf.gz", ".GRCh38.vcf.gz", _local_dragen_vcf)
    grch38_reject_output = re.sub(".vcf.gz", ".GRCh38.reject.vcf.gz", _local_dragen_vcf)
    liftover_command = '"cd /home && gatk LiftoverVcf -C b37ToHg38.over.chain -I {} -O {} -R /home2/Homo_sapiens_assembly38.fasta --REJECT {}"'.format(_local_dragen_vcf, grch38_output, grch38_reject_output)
#     !docker run -ti -v /home/rcarter/ComparingGatkSentieonDRAGEN/scratch:/home -v /home/rcarter:/home2 gcr.io/bioskryb/gatk:4.1.3.0 /bin/bash -c $liftover_command

# Compress, index, and split the Converted DRAGEN results (GRCh38)  into snps and indels

local_dragen_grch38_vcfs = !ls scratch/*.GRCh38.vcf.gz | xargs -i basename {}
local_dragen_grch38_vcfs

for _local_dragen_grch38_vcf in local_dragen_grch38_vcfs:
    snp_output = re.sub(".vcf.gz", ".snp.vcf.gz", _local_dragen_grch38_vcf)
    indel_output = re.sub(".vcf.gz", ".indel.vcf.gz", _local_dragen_grch38_vcf)
    gatk_command = '"cd /home && gatk SplitVcfs --VALIDATION_STRINGENCY SILENT --STRICT false -I {invcf} --INDEL_OUTPUT {indel} --SNP_OUTPUT {snp}"'.format(invcf = _local_dragen_grch38_vcf, indel = indel_output, snp = snp_output)
#     !docker run -v /home/rcarter/ComparingGatkSentieonDRAGEN/scratch/:/home -ti gcr.io/bioskryb/gatk:4.1.3.0 /bin/bash -c $gatk_command

sentieon_snp_vcfs = !ls scratch/*dnascope.snp.vcf.gz | xargs -i basename {}
dragen_hf_snp_vcfs = !ls scratch/*dragen*.snp.vcf.gz | grep hard-filtered | xargs -i basename {}
dragen_snp_vcfs = !ls scratch/*dragen*.snp.vcf.gz | grep -v hard-filtered | xargs -i basename {}

for _snp_vcf in sentieon_snp_vcfs:
    sample_name = re.sub("_merged.+", "", os.path.basename(_snp_vcf))
    folder_name = "NA12878_bulk2_vs_{}_sentieon_no_ploidy".format(sample_name)
#     !cd scratch/ && ~/rtg-tools/bin/rtg-tools-3.10.1-4d58eadb/rtg RTG_MEM=4G  vcfeval --squash-ploidy  -b /home/rcarter/R_bioskryb_na12878_subsampled_with_bulk_analysis/vumc_subsampled_with_bulk.merged.snps_and_indels.vqsr.pass.snp.NA12878_Bulk2.vcf.gz -c $_snp_vcf -o $folder_name -t /home/rcarter/R_bioskryb_na12878_analysis/NA12878_comparison/SDF/
for _snp_vcf in dragen_hf_snp_vcfs:
    sample_name = re.sub("_merged.+", "", os.path.basename(_snp_vcf))
    folder_name = "NA12878_bulk2_vs_{}_dragen_hf_no_ploidy".format(sample_name)
#     !cd scratch/ && ~/rtg-tools/bin/rtg-tools-3.10.1-4d58eadb/rtg RTG_MEM=4G  vcfeval --squash-ploidy  -b /home/rcarter/R_bioskryb_na12878_subsampled_with_bulk_analysis/vumc_subsampled_with_bulk.merged.snps_and_indels.vqsr.pass.snp.NA12878_Bulk2.vcf.gz -c $_snp_vcf -o $folder_name -t /home/rcarter/R_bioskryb_na12878_analysis/NA12878_comparison/SDF/
for _snp_vcf in dragen_snp_vcfs:
    sample_name = re.sub("_merged.+", "", os.path.basename(_snp_vcf))
    folder_name = "NA12878_bulk2_vs_{}_dragen_no_ploidy".format(sample_name)
#     !cd scratch/ && ~/rtg-tools/bin/rtg-tools-3.10.1-4d58eadb/rtg RTG_MEM=4G  vcfeval --squash-ploidy  -b /home/rcarter/R_bioskryb_na12878_subsampled_with_bulk_analysis/vumc_subsampled_with_bulk.merged.snps_and_indels.vqsr.pass.snp.NA12878_Bulk2.vcf.gz -c $_snp_vcf -o $folder_name -t /home/rcarter/R_bioskryb_na12878_analysis/NA12878_comparison/SDF/


