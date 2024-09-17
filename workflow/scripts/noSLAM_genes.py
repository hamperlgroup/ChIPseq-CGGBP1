import pandas as pd
import datatable as dt

bulk = dt.fread(
    "/lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1/data/SLAMseq_bulkRNAdata.csv"
).to_pandas()
nasc = dt.fread(
    "/lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1/data/SLAMseq_nascentRNAdata.csv"
).to_pandas()

genes_slam = set(list(bulk.GeneID) + list(nasc.GeneID))

genes_chip = dt.fread(
    "/lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1/results/ChIP-peaks/narrowPeak/merged/peak_annotation/annotation_ChIP-seq_CGGBP1_ENCODE.tsv"
).to_pandas()

diff_chip = list(set(genes_chip.SYMBOL).difference(genes_slam))

final_table = genes_chip.loc[genes_chip["SYMBOL"].isin(diff_chip)]

final_table.to_csv(
    "/lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1/results/ChIP-peaks/narrowPeak/merged/peak_annotation/noSLAM_genes_ChIP-CGGBP1.tsv",
    index=False,
    sep="\t",
)
