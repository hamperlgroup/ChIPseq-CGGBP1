############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    -----   GO term search Up & Down genes from Bulk and Nascent   -----    ###########
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############

############    -----------------------------------------    ############
### ----------------------- About the script ------------------------ ###
############    -----------------------------------------    ############
# Author: Elizabeth Marquez-Gomez
# Date: 2024-December
# Version: 1
# Subversion: 0
# Environment: downstream_RNAseq
# Using ...

#-----> Usage
# Rscript --vanilla scripts/deseq2.R {input} {output}


############    -----------------------------------------    ############
### --------------------------- Libraries --------------------------- ###
############    -----------------------------------------    ############

library(clusterProfiler)
library(org.Hs.eg.db)
library(simplifyEnrichment)
library(gridExtra)
library(ggplot2)
library(data.table)


############    -----------------------------------------    ############
### ----------------------- General Arguments ----------------------- ###
############    -----------------------------------------    ############

roothPath <- "/lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1"
prefixPlots <- "results/ChIP-peaks/narrowPeak/SLAM-seq_overlap/figures"

bulk <- "results/ChIP-peaks/narrowPeak/SLAM-seq_overlap/merged/SLAM_bulk_filtered.tsv"
nasc <- "results/ChIP-peaks/narrowPeak/SLAM-seq_overlap/merged/SLAM_nascent_filtered.tsv"


############    -----------------------------------------    ############
### -------------------------- Functions ---------------------------- ###
############    -----------------------------------------    ############

GOenrichment <- function(gene_names) {
    ## Adjust the p-values for multiple testing correction.
    ## This is essential because many enrichment tests are performed simultaneously, one for each GO term, and multiple testing increases the likelihood of false positives.
    ## pval default is 0.05, relaxed to 0.1 to get genes in nascent set
    enrichGO(gene_names, OrgDb = org.Hs.eg.db, ont = mode, keyType = "ENSEMBL", pAdjustMethod = "fdr", pvalueCutoff = 1)
}

PlotGO <- function(enrichGO, name) {
    dp <- dotplot(enrichGO)
    cn <- cnetplot(enrichGO)
    arrangeGrob(dp, cn, top = name, nrow = 2)
}

PlotClusteredGO <- function(gene_names, name, p_thr = 0.05) {
    enrichGO <- GOenrichment(gene_names)
    hits <- enrichGO@result$ID[enrichGO@result$p.adjust < p_thr]
    if (length(hits) > 1) {
        pdf(file = paste0(roothPath, "/", prefixPlots, "/GO_clustered_", name, ".pdf"), width = 10, height = 6)
        simplifyGO(GO_similarity(hits, db = "org.Hs.eg.db", ont = mode), column_title = name)
        dev.off()
    } else {
        print("Group has no or only one GO term sig. enriched")
    }
}

GOanalysis <- function(wt, mut, wt_name, mt_name, group, p_adj) {
    wildtype_plots <- PlotGO(GOenrichment(wt$ensembl), wt_name)
    mutated_plots <- PlotGO(GOenrichment(mut$ensembl), mt_name)

    go_plots <- arrangeGrob(mutated_plots, wildtype_plots, ncol = 2)
    ggsave(paste0(roothPath, "/", prefixPlots, "/GO_", group, ".pdf"), go_plots, width = 25, height = 20)

    PlotClusteredGO(mut$ensembl, mt_name)
    PlotClusteredGO(wt$ensembl, wt_name)

    genes <- list(mut$ensembl, wt$ensembl)
    names(genes) <- c(mt_name, wt_name)
    go_list <- lapply(genes, GOenrichment)
    go_list <- lapply(go_list, function(x) x$ID[x$p.adjust < p_adj])

    pdf(file = paste0(roothPath, "/", prefixPlots, "/GO_", group, "_comparative.pdf"), width = 10, height = 6)
    simplifyGOFromMultipleLists(go_list, padj_cutoff = p_adj, db = "org.Hs.eg.db", ont = mode)
    dev.off()

    return(0)
}

############    -----------------------------------------    ############
### ------------------------------ Main ----------------------------- ###
############    -----------------------------------------    ############

### ------------------------------ Inputs ----------------------------- ###
bulk <- as.data.frame(fread(bulk))
bulk <- bulk[bulk$class == "Up", ]
bulk$ensembl <- mapIds(
    x = org.Hs.eg.db,
    keys = bulk$GeneID,
    column = "ENSEMBL",
    keytype = "SYMBOL",
    multiVals = "first"
)

nasc <- as.data.frame(fread(nasc))
nasc <- nasc[nasc$class == "Up", ]
nasc$ensembl <- mapIds(
    x = org.Hs.eg.db,
    keys = nasc$GeneID,
    column = "ENSEMBL",
    keytype = "SYMBOL",
    multiVals = "first"
)

### ----------- Pathway analysis ----------- ###
mode <- "BP" # BP=biological process, CC=cellular component or MF=molecular function


GOanalysis(bulk, nasc, "Bulk", "Nascent", "bulk_nascent", 0.05)

length(intersect(nasc$GeneID, keys(org.Hs.eg.db, keytype = "SYMBOL")))


library(topGO)

gene_entrez <- na.omit(bulk$entrez)
all_genes <- keys(org.Hs.eg.db, keytype = "ENTREZID")
gene_universe <- factor(as.integer(all_genes %in% gene_entrez))
names(gene_universe) <- all_genes

GOdata <- new(
    "topGOdata",
    ontology = "BP", # Choose "BP", "MF", or "CC"
    allGenes = gene_universe,
    geneSelectionFun = function(x) x == 1,
    annot = annFUN.org,
    mapping = "org.Hs.eg.db",
    ID = "ENTREZID"
)
