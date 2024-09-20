############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    -----   Pie plots SLAM-seq bound/unbound  -----    ###########
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############

############    -----------------------------------------    ############
### ----------------------- About the script ------------------------ ###
############    -----------------------------------------    ############
# Author: Elizabeth Marquez-Gomez
# Date: 2024-June
# Version: 1
# Subversion: 0
# Environment: overlap_OK-DRIPc
# Using ...

#-----> Usage
# Rscript slam-seq_gene_overlap.R --sampleName [sample name] --inputFile [annot file]  --absCutoffFC [FC cutoff abs value] --cutoffPval [pval cutoff] --peakMode [broad | narrow] --sampleMode [merged | replicates]


############    -----------------------------------------    ############
### --------------------------- Libraries --------------------------- ###
############    -----------------------------------------    ############

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(annotables)
library(optparse)
library(data.table)
library(ggplot2)
library(ggrepel)
library(dplyr)

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

PlottingPieChart <- function(sample, rnaSet, rnaType) {
    peaks <- as.data.frame(fread(inputFile))

    ## Classification
    rnaSet <- cbind(rnaSet, class = "unchanged")
    rnaSet[rnaSet$log2FC < -absFC & rnaSet$negLog10pvalue > pval, ]$class <- "down"
    rnaSet[rnaSet$log2FC > absFC & rnaSet$negLog10pvalue > pval, ]$class <- "up"

    up <- cbind(rnaSet[rnaSet$class == "up", ], status = "unbound")
    up[up$GeneID %in% peaks$SYMBOL, ]$status <- "bound"

    down <- cbind(rnaSet[rnaSet$class == "down", ], status = "unbound")
    down[down$GeneID %in% peaks$SYMBOL, ]$status <- "bound"

    # Step 1: Summarize data by gene class
    up_summary <- up %>%
        group_by(status) %>%
        summarize(count = n()) %>%
        mutate(percentage = count / sum(count) * 100)

    # Step 2: Create the pie chart
    plot <- ggplot(up_summary, aes(x = "", y = percentage, fill = status)) +
        geom_bar(stat = "identity", width = 1, color = "white") +
        coord_polar(theta = "y") +
        theme_void() + # Removes background, axes, etc.
        labs(fill = "Gene Class", title = "Pie Chart of Gene Classes") +
        geom_label(aes(label = paste0(round(percentage, 1), "%")),
            color = "white",
            position = position_stack(vjust = 0.5),
            show.legend = FALSE
        ) + # Add percentage labels
        theme(legend.position = "right") +
        scale_fill_manual(values = c("#fabf7c", "#71b571"))

    pdf(paste0(outDir, "/figures/pieChart_up_SLAM-", rnaType, "_", sample, ".pdf"))
    print(plot)
    dev.off()

    # Step 1: Summarize data by gene class
    down_summary <- down %>%
        group_by(status) %>%
        summarize(count = n()) %>%
        mutate(percentage = count / sum(count) * 100)

    # Step 2: Create the pie chart
    plot <- ggplot(down_summary, aes(x = "", y = percentage, fill = status)) +
        geom_bar(stat = "identity", width = 1, color = "white") +
        coord_polar(theta = "y") +
        theme_void() + # Removes background, axes, etc.
        labs(fill = "Gene Class", title = "Pie Chart of Gene Classes") +
        geom_label(aes(label = paste0(round(percentage, 1), "%")),
            color = "white",
            position = position_stack(vjust = 0.5),
            show.legend = FALSE
        ) + # Add percentage labels
        theme(legend.position = "right") +
        scale_fill_manual(values = c("#fabf7c", "#5b85bc")) +
        annotate("text", x = 1, y = 10, label = length(down_summary????????))

    pdf(paste0(outDir, "/figures/pieChart_down_SLAM-", rnaType, "_", sample, ".pdf"))
    print(plot)
    dev.off()

    return()
}

############    -----------------------------------------    ############
### ------------------------------ Main ----------------------------- ###
############    -----------------------------------------    ############

outDir <- paste0(roothPath, "/results/ChIP-peaks/", peakMode, "/SLAM-seq_overlap")


PlottingPieChart(sampleName, bulk, "bulk")

PlottingPieChart(sampleName, nasc, "nascent")
