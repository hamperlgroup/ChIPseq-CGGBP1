############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
                ############    -----   Gene annotation in OK-seq peaks   -----    ###########
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############

############    -----------------------------------------    ############
### ----------------------- About the script ------------------------ ###
############    -----------------------------------------    ############
# Author: Elizabeth Marquez-Gomez
# Date: 2024-March
# Version: 1
# Subversion: 0
# Environment: jupyter
# Redefining regions from OK-seq & DRIPc-seq overlap for final overlap with ChIP-seq data

#-----> Usage
# python define_trc-ctrl_regions.py --sampleOK [sample ID] --dataDripc [dataset DRIPc]

### -------------------------- Run log -------------------------- ###
## May 4th 2023 ##


############    -----------------------------------------    ############
### --------------------------- Libraries --------------------------- ###
############    -----------------------------------------    ############
import pandas as pd
import datatable as dt
import os
import argparse


############    -----------------------------------------    ############
### --------------------------- Functions --------------------------- ###
############    -----------------------------------------    ############

### --------------------------- RedefineRegions --------------------------- ###
## Function to redefine regions
import pandas as pandas
import datatable as dt
import numpy as np

bulk = dt.fread('data/SLAMseq_bulkRNAdata.csv').to_pandas()
bulk.rename(columns={'p value':'log10-pval'}, inplace=True)
nasc = dt.fread('data/SLAMseq_nascentRNAdata.csv').to_pandas()
nasc.rename(columns={'p value':'log10-pval'}, inplace=True)

#   rnaSet <- filter(rnaSet, (log2FC > absFC | log2FC < -absFC) & pvalue > pval)

bulk = bulk[((bulk['log2 FC'] > 1) | (bulk['log2 FC'] < -1)) & (bulk['log10-pval'] > 1)]
bulk.to_csv('data/SLAMseq_bulkRNAdata_-1>FC>1_pval<1.csv', sep='\t', index=False)

nasc = nasc[((nasc['log2 FC'] > 1) | (nasc['log2 FC'] < -1)) & (nasc['log10-pval'] > 1)]
nasc.to_csv('data/SLAMseq_nascentRNAdata_-1>FC>1_pval<1.csv', sep='\t', index=False)

############    -----------------------------------------    ############
### ----------------------- General Arguments ----------------------- ###
############    -----------------------------------------    ############

############    -----------------------------------------    ############
### ------------------------------ Main ----------------------------- ###
############    -----------------------------------------    ############
