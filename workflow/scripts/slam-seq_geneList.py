import pandas as pandas
import datatable as dt
import numpy as np

bulk = dt.fread('data/SLAMseq_bulkRNAdata.csv').to_pandas()
bulk.rename(columns={'p value':'log10-pval'}, inplace=True)
nasc = dt.fread('data/SLAMseq_nascentRNAdata.csv').to_pandas()
nasc.rename(columns={'p value':'log10-pval'}, inplace=True)

bulk = bulk[((bulk['log2 FC'] > 1) | (bulk['log2 FC'] < -1)) & (bulk['log10-pval'] < 0.05)]
bulk.to_csv('data/SLAMseq_bulkRNAdata.csv', sep='\t', index=False, )
bulk = bulk[((bulk['log2 FC'] > 1) | (bulk['log2 FC'] < -1)) & (bulk['log10-pval'] < 0.05)]

