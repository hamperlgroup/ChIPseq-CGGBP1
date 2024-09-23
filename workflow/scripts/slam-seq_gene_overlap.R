############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    -----   Overlap ChIP annotation - SLAM-seq gene list   -----    ###########
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############

############    -----------------------------------------    ############
### ----------------------- About the script ------------------------ ###
############    -----------------------------------------    ############
# Author: Elizabeth Marquez-Gomez
# Date: 2024-June
# Version: 1
# Subversion: 0
# Environment: overlap_OK-DRIPc
# Using ChIP-seeker peaks will be annotated to identify genes in zones of interest from OK-seq analysis

#-----> Usage
# Rscript slam-seq_gene_overlap.R --sampleName [sample name] --inputFile [annot file]  --absCutoffFC [FC cutoff abs value] --cutoffPval [pval cutoff] --peakMode [broad | narrow] --sampleMode [merged | replicates]


############    -----------------------------------------    ############
### --------------------------- Libraries --------------------------- ###
############    -----------------------------------------    ############

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(annotables)
library(optparse)
library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)
library(data.table)
library(ggplot2)
library(ggrepel)

############    -----------------------------------------    ############
### ----------------------- General Arguments ----------------------- ###
############    -----------------------------------------    ############

option_list <- list(
  make_option(c("--sampleName"),
    type = "character",
    help = "Sample name ID", metavar = "character"
  ),
  make_option(c("--inputFile"),
    type = "character",
    help = "tsv file", metavar = "character"
  ),
  make_option(c("--absCutoffFC"),
    type = "double",
    help = "Absolute cutoff value to filter FC", metavar = "FC"
  ),
  make_option(c("--cutoffPval"),
    type = "double",
    help = "cutoff value to filter P-value", metavar = "Pval"
  ),
  make_option(c("--peakMode"),
    type = "character",
    help = "Peak mode", metavar = "character"
  ),
  make_option(c("--sampleMode"),
    type = "character",
    help = "Sample mode", metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
sampleName <- opt$sampleName
inputFile <- opt$inputFile
absFC <- opt$absCutoffFC
pval <- opt$cutoffPval
peakMode <- opt$peakMode
sampleMode <- opt$sampleMode


roothPath <- "/lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1"

## SLAM-seq gene list
bulk <- as.data.frame(fread("data/SLAMseq_bulkRNAdata.csv"))
nasc <- as.data.frame(fread("data/SLAMseq_nascentRNAdata.csv"))

# Rscript workflow/scripts/slam-seq_gene_overlap.R --sampleName ChIP-seq_CGGBP1_ENCODE --inputFile results/ChIP-peaks/narrowPeak/merged/peak_annotation/annotation_ChIP-seq_CGGBP1_ENCODE.tsv --absCutoffFC 1 --cutoffPval 1 --peakMode narrowPeak --sampleMode merged

sampleName <- "ChIP-seq_CGGBP1_ENCODE"
inputFile <- "results/ChIP-peaks/narrowPeak/merged/peak_annotation/annotation_ChIP-seq_CGGBP1_ENCODE.tsv"
absFC <- 1
pval <- 1
peakMode <- "narrowPeak"
sampleMode <- "merged"

## Wes Anderson palette - Zissou1
# wes_palette("Zissou1"),
# wes_palette("AsteroidCity1")[4], wes_palette("Royal2")[5]
c(
  "#3B9AB2", "#78B7C5", "#EBCC2A",
  "#E1AF00", "#F21A00", "#6C8645", "#74A089", "#b0afa2"
)

############    -----------------------------------------    ############
### --------------------------- Functions --------------------------- ###
############    -----------------------------------------    ############

PlottingScatter <- function(sample, rnaSet, rnaType) {
  ## Read peaks gene annotation
  peaks <- as.data.frame(fread(inputFile))

  ## Classification of SLAM-seq genes
  rnaSet <- cbind(rnaSet, class = "Unchanged")
  rnaSet[rnaSet$log2FC < -absFC & rnaSet$negLog10pvalue > pval, ]$class <- "Down"
  rnaSet[rnaSet$log2FC > absFC & rnaSet$negLog10pvalue > pval, ]$class <- "Up"

  ## Merge sets
  overlap <- merge(rnaSet, peaks,
    by.x = "GeneID", by.y = "SYMBOL"
  )
  overlap <- overlap[!duplicated(overlap$GeneID), ]

  ## Check if not empty
  empty <- character(0)
  if (identical(empty, overlap)) {
    fileConn <- file(paste0(outDir, "/overlap_", sample, "-SLAM_", rnaType, ".tsv"))
    writeLines(c("Empty overlap!"), fileConn)
    close(fileConn)
    return()
  } else {
    ## Save hits
    subset <- subset(overlap, class != "Unchanged")
    unchanged <- subset(overlap, class == "Unchanged")

    ## Scatter plot genes
    xlim_value <- ceiling(max(abs(overlap$log2FC)))

    plot <- ggplot(overlap, aes(x = log2FC, y = negLog10pvalue)) +
      geom_point(aes(
        colour = class,
        show.legend = FALSE, alpha = 0.7
      )) +
      scale_color_manual(values = c(
        "Unchanged" = "#b0afa2",
        "Down" = "#3B9AB2",
        "Up" = "#F21A00"
      )) +
      geom_text_repel(
        data = subset,
        aes(log2FC, negLog10pvalue, label = GeneID),
        size = 3,
        # min.segment.length = Inf,
        min.segment.length = 0,
        seed = 42,
        box.padding = 0.2,
        max.overlaps = Inf,
        nudge_x = .15,
        nudge_y = .5,
        color = "#b0afa2"
      ) +
      geom_hline(yintercept = pval, linetype = "dashed", color = "black") +
      geom_vline(xintercept = absFC, linetype = "dashed", color = "black") +
      geom_vline(xintercept = -absFC, linetype = "dashed", color = "black") +
      xlim(-xlim_value, xlim_value) +
      ylab("-log10 q-value") +
      xlab(paste0("log2 FC ", rnaType, " RNA")) +
      ggtitle(paste0("RNA ", rnaType, " - ", sample, " - ", peakMode)) +
      theme(
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")
      )

    # pdf(paste0(outDir, "/figures/SLAM-", rnaType, "_", sample, ".pdf"))
    pdf(paste0(outDir, "/figures/SLAM-", rnaType, "_", sample, "_lines.pdf"))
    print(plot)
    dev.off()

    write.table(subset, paste0(outDir, "/", sampleMode, "/overlap_", sample, "-SLAM_", rnaType, "_filtered.tsv"), sep = "\t", quote = F, row.names = F)
    write.table(unchanged, paste0(outDir, "/", sampleMode, "/overlap_", sample, "-SLAM_", rnaType, "_unchanged.tsv"), sep = "\t", quote = F, row.names = F)

    return()
  }
}

############    -----------------------------------------    ############
### ------------------------------ Main ----------------------------- ###
############    -----------------------------------------    ############

if (!dir.exists(paste0(roothPath, "/results/ChIP-peaks/", peakMode, "/SLAM-seq_overlap/figures"))) {
  dir.create(paste0(roothPath, "/results/ChIP-peaks/", peakMode, "/SLAM-seq_overlap/figures"))
}

if (!dir.exists(paste0(roothPath, "/results/ChIP-peaks/", peakMode, "/SLAM-seq_overlap/", sampleMode))) {
  dir.create(paste0(roothPath, "/results/ChIP-peaks/", peakMode, "/SLAM-seq_overlap/", sampleMode))
}


outDir <- paste0(roothPath, "/results/ChIP-peaks/", peakMode, "/SLAM-seq_overlap")


PlottingScatter(sampleName, bulk, "bulk")

PlottingScatter(sampleName, nasc, "nascent")
