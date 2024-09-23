# mamba activate givz_tracks
genes <- c("IRF2BPL", "RPS29", "C7orf50", "EGR1", "PBX1", "EIF3F", "SEC22B")

library(Gviz)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(wesanderson)

# library(biomaRt)
# library(Homo.sapiens)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

## Outdir
outDir <- "/lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1/results/genome_tracks"

## Wes Anderson palette - Zissou1
wes_palette("Zissou1")
c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")


###################### --------------> RPS29 <--------------######################
gene <- "RPS29"
# Load the BigWig file
bg_encode <- "results/covmean/ChIP-seq_CGGBP1_ENCODE.coverage.bedgraph"
bg_encode_track <- DataTrack(
    range = bg_encode,
    type = "hist", genome = "hg38",
    name = "ChIP-seq CGGBP1",
    col.histogram = "#78b7c5", fill.histogram = "#78b7c5",
    ylim = c(0, 0.6)
)

bg_ctrl <- "results/covmean/ChIP-seq_ctrl.coverage.bedgraph"
bg_ctrl_track <- DataTrack(
    range = bg_ctrl,
    type = "hist", genome = "hg38",
    name = "ChIP-seq CGGBP1 control",
    col.histogram = "#78b7c5", fill.histogram = "#78b7c5",
    ylim = c(0, 0.6)
)

# Select the gene of interest
gene_symbol <- gene
gene_info <- select(org.Hs.eg.db,
    keys = gene_symbol,
    columns = c("SYMBOL", "ENTREZID"),
    keytype = "SYMBOL"
)
entrez_id <- gene_info$ENTREZID
gene_range <- genes(txdb, filter = list(gene_id = entrez_id))

# Create the gene region track
gene_track <- GeneRegionTrack(txdb,
    chromosome = as.character(seqnames(gene_range)),
    start = start(gene_range) - 5000, end = end(gene_range) + 5000,
    name = "Gene region", strand = strand(gene_range),
    collapseTranscripts = "meta",
    transcriptAnnotation = "gene", shape = "smallArrow"
)
# Load chromosome
ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = gene_track@chromosome)
axTrack <- GenomeAxisTrack()

# Plot the tracks
pdf(paste0(outDir, "/", gene, "_track.pdf"))
plotTracks(list(ideoTrack, axTrack, gene_track, bg_encode_track, bg_ctrl_track),
    from = gene_track@start - 5000, to = gene_track@end + 5000,
    chromosome = gene_track@chromosome, showBandId = TRUE, cex.bands = 0.9
)
dev.off()

###################### --------------> C7orf50 <--------------######################
gene <- "C7orf50"
# Load the BigWig file
bg_encode <- "results/covmean/ChIP-seq_CGGBP1_ENCODE.coverage.bedgraph"
bg_encode_track <- DataTrack(
    range = bg_encode,
    type = "hist", genome = "hg38",
    name = "ChIP-seq CGGBP1",
    col.histogram = "#78b7c5", fill.histogram = "#78b7c5",
    ylim = c(0, 0.6)
)

bg_ctrl <- "results/covmean/ChIP-seq_ctrl.coverage.bedgraph"
bg_ctrl_track <- DataTrack(
    range = bg_ctrl,
    type = "hist", genome = "hg38",
    name = "ChIP-seq CGGBP1 control",
    col.histogram = "#78b7c5", fill.histogram = "#78b7c5",
    ylim = c(0, 0.6)
)

# Select the gene of interest
gene_symbol <- gene
gene_info <- select(org.Hs.eg.db,
    keys = gene_symbol,
    columns = c("SYMBOL", "ENTREZID"),
    keytype = "SYMBOL"
)
entrez_id <- gene_info$ENTREZID
gene_range <- genes(txdb, filter = list(gene_id = entrez_id))

# Create the gene region track
gene_track <- GeneRegionTrack(txdb,
    chromosome = as.character(seqnames(gene_range)),
    start = start(gene_range) - 7000, end = end(gene_range) + 7000,
    name = "Gene region", strand = strand(gene_range),
    collapseTranscripts = "meta",
    transcriptAnnotation = "gene", shape = "smallArrow"
)
# Load chromosome
ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = gene_track@chromosome)
axTrack <- GenomeAxisTrack()

# Plot the tracks
pdf(paste0(outDir, "/", gene, "_track.pdf"))
plotTracks(list(ideoTrack, axTrack, gene_track, bg_encode_track, bg_ctrl_track),
    from = gene_track@start - 7000, to = gene_track@end + 7000,
    chromosome = gene_track@chromosome, showBandId = TRUE, cex.bands = 0.9
)
dev.off()

###################### --------------> IRF2BPL <--------------######################
gene <- "IRF2BPL"
# Load the BigWig file
bg_encode <- "results/covmean/ChIP-seq_CGGBP1_ENCODE.coverage.bedgraph"
bg_encode_track <- DataTrack(
    range = bg_encode,
    type = "hist", genome = "hg38",
    name = "ChIP-seq CGGBP1",
    col.histogram = "#78b7c5", fill.histogram = "#78b7c5",
    ylim = c(0, 1)
)

bg_ctrl <- "results/covmean/ChIP-seq_ctrl.coverage.bedgraph"
bg_ctrl_track <- DataTrack(
    range = bg_ctrl,
    type = "hist", genome = "hg38",
    name = "ChIP-seq CGGBP1 control",
    col.histogram = "#78b7c5", fill.histogram = "#78b7c5",
    ylim = c(0, 1)
)

# Select the gene of interest
gene_symbol <- gene
gene_info <- select(org.Hs.eg.db,
    keys = gene_symbol,
    columns = c("SYMBOL", "ENTREZID"),
    keytype = "SYMBOL"
)
entrez_id <- gene_info$ENTREZID
gene_range <- genes(txdb, filter = list(gene_id = entrez_id))

# Create the gene region track
gene_track <- GeneRegionTrack(txdb,
    chromosome = as.character(seqnames(gene_range)),
    start = start(gene_range) - 7000, end = end(gene_range) + 7000,
    name = "Gene region", strand = strand(gene_range),
    collapseTranscripts = "meta",
    transcriptAnnotation = "gene", shape = "smallArrow"
)
# Load chromosome
ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = gene_track@chromosome)
axTrack <- GenomeAxisTrack()

# Plot the tracks
pdf(paste0(outDir, "/", gene, "_track.pdf"))
plotTracks(list(ideoTrack, axTrack, gene_track, bg_encode_track, bg_ctrl_track),
    from = gene_track@start - 7000, to = gene_track@end + 7000,
    chromosome = gene_track@chromosome, showBandId = TRUE, cex.bands = 0.9
)
dev.off()

###################### --------------> EGR1 <--------------######################
gene <- "EGR1"
# Load the BigWig file
bg_encode <- "results/covmean/ChIP-seq_CGGBP1_ENCODE.coverage.bedgraph"
bg_encode_track <- DataTrack(
    range = bg_encode,
    type = "hist", genome = "hg38",
    name = "ChIP-seq CGGBP1",
    col.histogram = "#78b7c5", fill.histogram = "#78b7c5",
    ylim = c(0, 1)
)

bg_ctrl <- "results/covmean/ChIP-seq_ctrl.coverage.bedgraph"
bg_ctrl_track <- DataTrack(
    range = bg_ctrl,
    type = "hist", genome = "hg38",
    name = "ChIP-seq CGGBP1 control",
    col.histogram = "#78b7c5", fill.histogram = "#78b7c5",
    ylim = c(0, 1)
)

# Select the gene of interest
gene_symbol <- gene
gene_info <- select(org.Hs.eg.db,
    keys = gene_symbol,
    columns = c("SYMBOL", "ENTREZID"),
    keytype = "SYMBOL"
)
entrez_id <- gene_info$ENTREZID
gene_range <- genes(txdb, filter = list(gene_id = entrez_id))

# Create the gene region track
gene_track <- GeneRegionTrack(txdb,
    chromosome = as.character(seqnames(gene_range)),
    start = start(gene_range) - 19000, end = end(gene_range) + 19000,
    name = "Gene region", strand = strand(gene_range),
    collapseTranscripts = "meta",
    transcriptAnnotation = "gene", shape = "smallArrow"
)
# Load chromosome
ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = gene_track@chromosome)
axTrack <- GenomeAxisTrack()

# Plot the tracks
pdf(paste0(outDir, "/", gene, "_track.pdf"))
plotTracks(list(ideoTrack, axTrack, gene_track, bg_encode_track, bg_ctrl_track),
    from = gene_track@start - 19000, to = gene_track@end + 19000,
    chromosome = gene_track@chromosome, showBandId = TRUE, cex.bands = 0.9
)
dev.off()


###################### --------------> PBX1 <--------------######################
gene <- "PBX1"
# Load the BigWig file
bg_encode <- "results/covmean/ChIP-seq_CGGBP1_ENCODE.coverage.bedgraph"
bg_encode_track <- DataTrack(
    range = bg_encode,
    type = "hist", genome = "hg38",
    name = "ChIP-seq CGGBP1",
    col.histogram = "#78b7c5", fill.histogram = "#78b7c5",
    ylim = c(0, 1)
)

bg_ctrl <- "results/covmean/ChIP-seq_ctrl.coverage.bedgraph"
bg_ctrl_track <- DataTrack(
    range = bg_ctrl,
    type = "hist", genome = "hg38",
    name = "ChIP-seq CGGBP1 control",
    col.histogram = "#78b7c5", fill.histogram = "#78b7c5",
    ylim = c(0, 1)
)

# Select the gene of interest
gene_symbol <- gene
gene_info <- select(org.Hs.eg.db,
    keys = gene_symbol,
    columns = c("SYMBOL", "ENTREZID"),
    keytype = "SYMBOL"
)
entrez_id <- gene_info$ENTREZID
gene_range <- genes(txdb, filter = list(gene_id = entrez_id))

# Create the gene region track
gene_track <- GeneRegionTrack(txdb,
    chromosome = as.character(seqnames(gene_range)),
    start = start(gene_range) - 90000, end = end(gene_range) + 90000,
    name = "Gene region", strand = strand(gene_range),
    collapseTranscripts = "meta",
    transcriptAnnotation = "gene", shape = "smallArrow"
)
# Load chromosome
ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = gene_track@chromosome)
axTrack <- GenomeAxisTrack()

# Plot the tracks
pdf(paste0(outDir, "/", gene, "_track.pdf"))
plotTracks(list(ideoTrack, axTrack, gene_track, bg_encode_track, bg_ctrl_track),
    from = gene_track@start - 90000, to = gene_track@end + 90000,
    chromosome = gene_track@chromosome, showBandId = TRUE, cex.bands = 0.9
)
dev.off()


###################### --------------> EIF3F <--------------######################
gene <- "EIF3F"
# Load the BigWig file
bg_encode <- "results/covmean/ChIP-seq_CGGBP1_ENCODE.coverage.bedgraph"
bg_encode_track <- DataTrack(
    range = bg_encode,
    type = "hist", genome = "hg38",
    name = "ChIP-seq CGGBP1",
    col.histogram = "#78b7c5", fill.histogram = "#78b7c5",
    ylim = c(0, 0.8)
)

bg_ctrl <- "results/covmean/ChIP-seq_ctrl.coverage.bedgraph"
bg_ctrl_track <- DataTrack(
    range = bg_ctrl,
    type = "hist", genome = "hg38",
    name = "ChIP-seq CGGBP1 control",
    col.histogram = "#78b7c5", fill.histogram = "#78b7c5",
    ylim = c(0, 0.8)
)

# Select the gene of interest
gene_symbol <- gene
gene_info <- select(org.Hs.eg.db,
    keys = gene_symbol,
    columns = c("SYMBOL", "ENTREZID"),
    keytype = "SYMBOL"
)
entrez_id <- gene_info$ENTREZID
gene_range <- genes(txdb, filter = list(gene_id = entrez_id))

# Create the gene region track
gene_track <- GeneRegionTrack(txdb,
    chromosome = as.character(seqnames(gene_range)),
    start = start(gene_range) - 10000, end = end(gene_range) + 10000,
    name = "Gene region", strand = strand(gene_range),
    collapseTranscripts = "meta",
    transcriptAnnotation = "gene", shape = "smallArrow"
)
# Load chromosome
ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = gene_track@chromosome)
axTrack <- GenomeAxisTrack()

# Plot the tracks
pdf(paste0(outDir, "/", gene, "_track.pdf"))
plotTracks(list(ideoTrack, axTrack, gene_track, bg_encode_track, bg_ctrl_track),
    from = gene_track@start - 10000, to = gene_track@end + 10000,
    chromosome = gene_track@chromosome, showBandId = TRUE, cex.bands = 0.9
)
dev.off()


###################### --------------> SEC22B <--------------######################
gene <- "SEC22B"
# Load the BigWig file
bg_encode <- "results/covmean/ChIP-seq_CGGBP1_ENCODE.coverage.bedgraph"
bg_encode_track <- DataTrack(
    range = bg_encode,
    type = "hist", genome = "hg38",
    name = "ChIP-seq CGGBP1",
    col.histogram = "#78b7c5", fill.histogram = "#78b7c5",
    ylim = c(0, 0.8)
)

bg_ctrl <- "results/covmean/ChIP-seq_ctrl.coverage.bedgraph"
bg_ctrl_track <- DataTrack(
    range = bg_ctrl,
    type = "hist", genome = "hg38",
    name = "ChIP-seq CGGBP1 control",
    col.histogram = "#78b7c5", fill.histogram = "#78b7c5",
    ylim = c(0, 0.8)
)

# Select the gene of interest
gene_symbol <- gene
gene_info <- select(org.Hs.eg.db,
    keys = gene_symbol,
    columns = c("SYMBOL", "ENTREZID"),
    keytype = "SYMBOL"
)
entrez_id <- gene_info$ENTREZID
gene_range <- genes(txdb, filter = list(gene_id = entrez_id))

# Create the gene region track
gene_track <- GeneRegionTrack(txdb,
    chromosome = as.character(seqnames(gene_range)),
    start = start(gene_range) - 70000, end = end(gene_range) + 70000,
    name = "Gene region", strand = strand(gene_range),
    collapseTranscripts = "meta",
    transcriptAnnotation = "gene", shape = "smallArrow"
)
# Load chromosome
ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = gene_track@chromosome)
axTrack <- GenomeAxisTrack()

# Plot the tracks
pdf(paste0(outDir, "/", gene, "_track.pdf"))
plotTracks(list(ideoTrack, axTrack, gene_track, bg_encode_track, bg_ctrl_track),
    from = gene_track@start - 70000, to = gene_track@end + 70000,
    chromosome = gene_track@chromosome, showBandId = TRUE, cex.bands = 0.9
)
dev.off()
