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
# Rscript slam-seq_gene_overlap.R --sampleName [sample name] --inputFile [annot file]  --absCutoffFC [FC cutoff abs value] --cutoffPval [pval cutoff] --peakMode [broad | narrow]


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

############    -----------------------------------------    ############
### --------------------------- Functions --------------------------- ###
############    -----------------------------------------    ############

PlottingScatter <- function(sample, rnaSet, rnaType) {
  ## Convert gene symbols to Entrez IDs
  gene_ids <- mapIds(org.Hs.eg.db, keys = rnaSet$GeneID, column = "ENTREZID", keytype = "SYMBOL")

  # Extract gene IDs from the peak annotations
  peaks <- as.data.frame(fread(inputFile))
  annotated_genes <- peaks$geneId

  # Find overlap with SLAM-seq gene list
  overlapping_genes <- intersect(gene_ids, annotated_genes)

  empty <- character(0)
  if (identical(empty, overlapping_genes)) {
    fileConn <- file(paste0(outDir, "/overlap_", sample, "-SLAM_", rnaType, ".tsv"))
    writeLines(c("Empty overlap!"), fileConn)
    close(fileConn)
    return()
  } else {
    overlapping_gene_symbols <- mapIds(org.Hs.eg.db, keys = overlapping_genes, column = "SYMBOL", keytype = "ENTREZID")

    # Print overlapping genes
    overlapping_gene_symbols <- as.data.frame(overlapping_gene_symbols)
    overlapping_gene_symbols$EntrezID <- rownames(overlapping_gene_symbols)
    colnames(overlapping_gene_symbols)[colnames(overlapping_gene_symbols) == "overlapping_gene_symbols"] <- "GeneID"
    overlap <- merge(rnaSet, overlapping_gene_symbols, by = "GeneID")


    subset <- subset(overlap, (log2FC > absFC | log2FC < -absFC) & pvalue > pval)
    ## Scatter plot genes
    plot <- ggplot(overlap, aes(x = log2FC, y = pvalue)) +
      geom_point(aes(colour = (log2FC > absFC | log2FC < -absFC) & pvalue > pval),
        show.legend = FALSE, alpha = 0.5
      ) +
      geom_text(
        data = subset,
        aes(log2FC, pvalue, label = GeneID),
        nudge_x = 0.25, nudge_y = 0.25,
        check_overlap = T
      ) +
      geom_hline(yintercept = pval, linetype = "dashed", color = "red") +
      geom_vline(xintercept = absFC, linetype = "dashed", color = "blue") +
      geom_vline(xintercept = -absFC, linetype = "dashed", color = "blue") +
      ylab("-log10 q-value") +
      xlab(paste0("log2 FC ", rnaType, " RNA")) +
      ggtitle(paste0("RNA ", rnaType, " - ", sample, " - ", peakMode))

    pdf(paste0(outDir, "/figures/SLAM-", rnaType, "_", sample, ".pdf"))
    print(plot)
    dev.off()

    ## Save hits
    unchanged <- overlap[!overlap$GeneID %in% subset$GeneID, ]

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
