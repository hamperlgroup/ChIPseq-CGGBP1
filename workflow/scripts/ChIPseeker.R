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
# Environment: chipSeeker
# Using ChIP-seeker peaks will be annotated to identify genes in zones of interest from OK-seq analysis

#-----> Usage
# Rscript ChIPseeker.R --sampleName [sample name] --inputFile [bed file]


############    -----------------------------------------    ############
### --------------------------- Libraries --------------------------- ###
############    -----------------------------------------    ############

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(annotables)
library(optparse)

############    -----------------------------------------    ############
### ----------------------- General Arguments ----------------------- ###
############    -----------------------------------------    ############

option_list = list(
  make_option(c("--sampleName"), type="character",
              help="Sample name ID", metavar="character"),
  make_option(c("--inputFile"), type="character",
              help="BED file", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

sampleName <- opt$sampleName
inputFile <- opt$inputFile

roothPath <- '/lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1'

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

############    -----------------------------------------    ############
### ------------------------------ Main ----------------------------- ###
############    -----------------------------------------    ############

if (!dir.exists(paste0(roothPath,"/results/ChIP-peaks/peak_annotation"))){
    dir.create(paste0(roothPath,"/results/ChIP-peaks/peak_annotation"))}

## Peak annotation
peakAnno <- annotatePeak(inputFile, overlap='all',
                         TxDb=txdb, annoDb="org.Hs.eg.db")

## Pie chart of genomic region annotation
pdf(paste0(roothPath,'/results/ChIP-peaks/peak_annotation/pieChart_genomicAnnot_',sampleName,'.pdf'))
  plotAnnoPie(peakAnno)
dev.off()

## Distribution of TF-binding loci relative to TSS
pdf(paste0(roothPath,'/results/ChIP-peaks/peak_annotation/distTSSbind_',sampleName,'.pdf'))
  plotDistToTSS(peakAnno, title="Distribution of transcription factor-binding loci \n relative to TSS")
dev.off()

## Get peaks annotation and adding gene symbols
annoTable <- as.data.frame(peakAnno)

# Get unique entrez gene Ids
entrezids <- unique(annoTable$geneId)

# Get hg38 entrez to gene symbol mappings
entrez2gene <- grch38 %>% 
  filter(entrez %in% entrezids) %>% 
  dplyr::select(entrez, symbol)

# Match to each annotation dataframe
m <- match(annoTable$geneId, entrez2gene$entrez)
annoTable <- cbind(annoTable[,1:14], geneSymbol=entrez2gene$symbol[m], annoTable[,15:16])


write.table(annoTable, paste0(roothPath,'/results/ChIP-peaks/peak_annotation/peaks_',sampleName,'.tsv'), sep='\t', quote=F, row.names=F)

## Write only regions features' regions
annoTable$geneChr <- sub("^", "chr", annoTable$geneChr )
write.table(unique(annoTable[,c('geneChr','geneStart','geneEnd')]),paste0(roothPath,'/results/ChIP-peaks/peak_annotation/peaks_regions_',sampleName,'.bed'), row.names = F, col.names = F, sep='\t', quote=F,)
