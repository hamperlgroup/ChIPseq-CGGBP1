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
library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)

############    -----------------------------------------    ############
### ----------------------- General Arguments ----------------------- ###
############    -----------------------------------------    ############

option_list = list(
  make_option(c("--sampleName"), type="character",
              help="Sample name ID", metavar="character"),
  make_option(c("--inputFile"), type="character",
              help="BED file", metavar="character"),
  make_option(c("--peakMode"), type="character",
              help="Peak mode", metavar="character")    
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

sampleName <- opt$sampleName
inputFile <- opt$inputFile
peakMode <- opt$peakMode

roothPath <- '/lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1'

# Load annotation data (TxDb example)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

############    -----------------------------------------    ############
### ------------------------------ Main ----------------------------- ###
############    -----------------------------------------    ############

if (!dir.exists(paste0(roothPath,"/results/ChIP-peaks/",peakMode,"/peak_annotation"))){
    dir.create(paste0(roothPath,"/results/ChIP-peaks/",peakMode,"/peak_annotation"))}

outDir <- paste0(roothPath,"/results/ChIP-peaks/",peakMode,"/peak_annotation")
## Ensembl to UCSC notation in peaks file
# Import the peak file. Function usually returns a GRanges object directly.
peaks <- import(inputFile)

# Check the class of the imported object
print(class(peaks))

# If not already a GRanges object, convert it
if (!inherits(peaks, "GRanges")) {
    peaks <- as(peaks, "GRanges")
}

# Check the current sequence levels style
print(seqlevelsStyle(peaks))

# Set the desired sequence levels style, if necessary
seqlevelsStyle(peaks) <- "UCSC"  # or "NCBI"


## Peak annotation
peakAnno <- annotatePeak(peaks, overlap='all',
                         TxDb=txdb, annoDb="org.Hs.eg.db")

## Pie chart of genomic region annotation
pdf(paste0(outDir,'/pieChart_genomicAnnot_',sampleName,'.pdf'))
  plotAnnoPie(peakAnno)
  # Add a title to the plot
  title(main=paste0("Proportion of genomic regions - ",sampleName))
dev.off()

## Distribution of TF-binding loci relative to TSS
pdf(paste0(outDir,'/distTSSbind_',sampleName,'.pdf'))
  plotDistToTSS(peakAnno, title=paste0("Distribution of transcription factor-binding loci \n relative to TSS - ",sampleName," - ",peakMode))
dev.off()

## Get peaks annotation and adding gene symbols
annoTable <- as.data.frame(peakAnno)

# # Get unique entrez gene Ids
# entrezids <- unique(annoTable$geneId)

# # Get hg38 entrez to gene symbol mappings
# entrez2gene <- grch38 %>% 
#   filter(entrez %in% entrezids) %>% 
#   dplyr::select(entrez, symbol)

# # Match to each annotation dataframe
# m <- match(annoTable$geneId, entrez2gene$entrez)
# annoTable <- cbind(annoTable[,1:14], geneSymbol=entrez2gene$symbol[m], annoTable[,15:16])


write.table(annoTable, paste0(outDir,'/annotation_',sampleName,'.tsv'), sep='\t', quote=F, row.names=F)

## Write only regions features' regions
# annoTable$geneChr <- sub("^", "chr", annoTable$geneChr )
# write.table(unique(annoTable[,c('geneChr','geneStart','geneEnd')]),paste0(roothPath,'/results/ChIP-peaks/peak_annotation/features_regions_',sampleName,'.bed'), row.names = F, col.names = F, sep='\t', quote=F,)
