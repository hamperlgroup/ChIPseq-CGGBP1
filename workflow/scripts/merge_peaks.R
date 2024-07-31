############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
               ############    -----   Merge replicates ChIP-seq peaks   -----    ###########
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############

############    -----------------------------------------    ############
### ----------------------- About the script ------------------------ ###
############    -----------------------------------------    ############
# Author: Elizabeth Marquez-Gomez
# Date: 2024-June
# Version: 1
# Subversion: 0
# Environment: chipSeeker
# Using ChIP-seeker peaks will be annotated to identify genes in zones of interest from OK-seq analysis

#-----> Usage
# Rscript protein_coding_genes_overlap.R --sampleName [sample name] --peaksFile [peak file]  --genesFile [protein coding genes file]


############    -----------------------------------------    ############
### --------------------------- Libraries --------------------------- ###
############    -----------------------------------------    ############

library(GenomicRanges)
library(rtracklayer)
library(Gviz)
library(optparse)


############    -----------------------------------------    ############
### ----------------------- General Arguments ----------------------- ###
############    -----------------------------------------    ############

option_list = list(
  make_option(c("--groupName"), type="character",
              help="Sample name ID", metavar="character"),
  make_option(c("--peaksFiles"), type="character",
              help="narrowPeak file", metavar="character")  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

groupName <- opt$groupName
peaksFiles <- unlist(strsplit(opt$peaksFiles, ','))

roothPath <- '/lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1'

############    -----------------------------------------    ############
### --------------------------- Functions --------------------------- ###
############    -----------------------------------------    ############


############    -----------------------------------------    ############
### ------------------------------ Main ----------------------------- ###
############    -----------------------------------------    ############

if (!dir.exists(paste0(roothPath,"/results/ChIP-peaks/narrowPeak/merged"))){
    dir.create(paste0(roothPath,"/results/ChIP-peaks/narrowPeak/merged"))}

outDir <- paste0(roothPath,"/results/ChIP-peaks/narrowPeak/merged")


# Import the narrowPeak files
peak_list <- lapply(peaksFiles, function(path) {
  import(path, format = "narrowPeak")
})

# Combine all peaks into one GRanges object
combined_peaks <- do.call(c, peak_list)

# Merge overlapping or adjacent peaks
merged_peaks <- reduce(combined_peaks)

# Export the merged peaks to a narrowPeak file
export(merged_peaks, paste0(outDir,'/merged_',groupName,'.bed'), format = "BED")
