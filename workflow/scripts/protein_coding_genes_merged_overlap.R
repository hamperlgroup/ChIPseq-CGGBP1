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
    make_option(c("--mergedPeaks"),
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
mergedFile <- opt$mergedPeaks
genesFile <- opt$genesFile

roothPath <- "/lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1"

############    -----------------------------------------    ############
### --------------------------- Functions --------------------------- ###
############    -----------------------------------------    ############


############    -----------------------------------------    ############
### ------------------------------ Main ----------------------------- ###
############    -----------------------------------------    ############

if (!dir.exists(paste0(roothPath, "/results/overlap_protein_coding_genes"))) {
    dir.create(paste0(roothPath, "/results/overlap_protein_coding_genes"))
}

outDir <- paste0(roothPath, "/results/overlap_protein_coding_genes")

# Import the peaks file
mergedPeaks <- import(mergedFile, format = "BED")

# Import the BED file
genes <- import(genesFile, format = "BED")

# Find overlaps
overlaps <- subsetByOverlaps(genes, mergedPeaks)
names(overlaps) <- overlaps$name

export(granges(overlaps), paste0(outDir, "/overlapping_protein_coding_genes_", sampleName, ".bed"), format = "BED")


## nonOverlapping
difference <- genes[!genes$name %in% overlaps$name]
names(difference) <- difference$name

export(granges(difference), paste0(outDir, "/nonOverlapping_protein_coding_genes_", sampleName, ".bed"), format = "BED")
