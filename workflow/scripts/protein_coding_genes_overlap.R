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

option_list = list(
  make_option(c("--sampleName"), type="character",
              help="Sample name ID", metavar="character"),
  make_option(c("--peaksFile"), type="character",
              help="narrowPeak file", metavar="character"),
  make_option(c("--genesFile"), type="character",
              help="bed file", metavar="character")              
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

sampleName <- opt$sampleName
peaksFile <- opt$peaksFile
protCodGenesFile <- opt$genesFile

roothPath <- '/lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1'

############    -----------------------------------------    ############
### --------------------------- Functions --------------------------- ###
############    -----------------------------------------    ############


############    -----------------------------------------    ############
### ------------------------------ Main ----------------------------- ###
############    -----------------------------------------    ############

if (!dir.exists(paste0(roothPath,"/results/overlap_coding_genes"))){
    dir.create(paste0(roothPath,"/results/overlap_coding_genes"))}

outDir <- paste0(roothPath,"/results/overlap_coding_genes")

# Import the peaks file
peaks <- import(peaksFile, format = "narrowPeak")

# Import the BED file
genes <- import(protCodGenesFile, format = "BED")

# Find overlaps
overlaps <- findOverlaps(genes, peaks)

# Extract overlapping regions
overlappingPeaks <- peaks[subjectHits(overlaps)]
overlappingGenes <- genes[queryHits(overlaps)]

export(overlappingGenes, paste0(outDir,'/overlapping_coding_genes_',sampleName,'.bed'), format = "BED")

difference <- setdiff(genes, overlappingGenes)
export(difference, paste0(outDir,'/notOverlapping_coding_genes_',sampleName,'.bed'), format = "BED")

export(overlappingPeaks, paste0(outDir,'/overlapping_coding_genes-peaks_',sampleName,'.narrowPeak'), format = "narrowPeak")


## Additional Analysis

# Summarize the overlapping regions
summarizeOverlaps <- function(overlappingPeaks, overlappingGenes) {
  overlapsDF <- data.frame(
    gene = overlappingGenes$name,
    gene_start = start(overlappingGenes),
    gene_end = end(overlappingGenes),
    peak_start = start(overlappingPeaks),
    peak_end = end(overlappingPeaks),
    peak_score = mcols(overlappingPeaks)$score  # Assuming 'score' is a column in the narrowPeak file
  )
  return(overlapsDF)
}

overlapSummary <- summarizeOverlaps(overlappingPeaks, overlappingGenes)
write.table(as.data.frame(overlapSummary), paste0(outDir,'/summary_',sampleName,'.tsv'), sep='\t', quote=F, row.names=F)



# Visualize overlapping regions

# Create a Genome Axis Track
genomeAxisTrack <- GenomeAxisTrack()
my_chromosomes <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X")


# Create a Data Track for the narrowPeak data
peak_filtered <- keepSeqlevels(peaks, my_chromosomes, pruning.mode = "coarse")
peakTrack <- AnnotationTrack(range = peak_filtered, genome = "hg38", name = "narrowPeak")

# Create a Gene Region Track for the genes
# Filter the genes object to keep only the specified chromosomes
genes_filtered <- keepSeqlevels(genes, my_chromosomes, pruning.mode = "coarse")
genesTrack <- AnnotationTrack(range = genes_filtered, genome = "hg38", name = "Genes")

# Plot the tracks
pdf(paste0(outDir,'/overlappingRegions_',sampleName,'.pdf'))
    plotTracks(list(genomeAxisTrack, peakTrack, genesTrack))
dev.off()