############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    -----   Overlap coverage files with protein coding genes list   -----    ###########
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

option_list <- list(
  make_option(c("--sampleName"),
    type = "character",
    help = "Sample name ID", metavar = "character"
  ),
  make_option(c("--mergedCGGBP1"),
    type = "character",
    help = "bed mergef file", metavar = "character"
  ),
  make_option(c("--genesFile"),
    type = "character",
    help = "bed file", metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

sampleName <- opt$sampleName
mergedCGGBP1File <- opt$mergedCGGBP1
genesFile <- opt$genesFile

roothPath <- "/lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1"

############    -----------------------------------------    ############
### --------------------------- Functions --------------------------- ###
############    -----------------------------------------    ############


############    -----------------------------------------    ############
### ------------------------------ Main ----------------------------- ###
############    -----------------------------------------    ############

if (!dir.exists(paste0(roothPath, "/results/overlap_CGGBP1_all_genes"))) {
  dir.create(paste0(roothPath, "/results/overlap_CGGBP1_all_genes"))
}

outDir <- paste0(roothPath, "/results/overlap_CGGBP1_all_genes")

# Import the peaks file
mergedCGGBP1 <- import(mergedCGGBP1File, format = "BED")

# Import the BED file
genes <- import(genesFile, format = "gtf")

# Find overlaps
overlaps <- subsetByOverlaps(genes, mergedCGGBP1)
names(overlaps) <- overlaps$gene_id

export(granges(overlaps), paste0(outDir, "/overlapping_CGGBP1_all_genes_", sampleName, ".bed"), format = "BED")


## nonOverlapping
difference <- genes[!genes$gene_id %in% overlaps$gene_id]
names(difference) <- difference$gene_id

export(granges(difference), paste0(outDir, "/nonOverlapping_CGGBP1_all_genes_", sampleName, ".bed"), format = "BED")
