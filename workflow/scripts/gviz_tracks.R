# mamba activate givz_tracks
library(Gviz)
library(GenomicRanges)

data(cpgIslands)
class(cpgIslands)

chr <- as.character(unique(seqnames(cpgIslands)))
gen <- genome(cpgIslands)
atrack <- AnnotationTrack(cpgIslands, name = "CpG")

gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = gen, chromosome = chr)

pdf("example.pdf")
plotTracks(list(itrack, gtrack, atrack))
dev.off()


data(geneModels)
grtrack <- GeneRegionTrack(geneModels, genome = gen, chromosome = chr, name = "Gene Model")
pdf("example.pdf")
plotTracks(list(itrack, gtrack, atrack, grtrack))
dev.off()


library(BSgenome.Hsapiens.UCSC.hg38)
strack <- SequenceTrack(Hsapiens, chromosome = chr)
pdf("example.pdf")
plotTracks(list(itrack, gtrack, atrack, grtrack, strack), from = 26591822, to = 26591852, cex = 0.8)
dev.off()

grtrack <- GeneRegionTrack(geneModels,
    genome = gen,
    chromosome = chr, name = "Gene Model", transcriptAnnotation = "symbol",
    background.title = "brown"
)
head(displayPars(grtrack))
displayPars(grtrack) <- list(background.panel = "#FFFEDB", col = NULL)
head(displayPars(grtrack))
pdf("example.pdf")
plotTracks(list(itrack, gtrack, atrack, grtrack))
dev.off()



from <- 77022543
to <- 77030708
knownGenes <- UcscTrack(
    genome = "hg38", chromosome = "chr14",
    track = "knownGene", from = from, to = to, trackType = "GeneRegionTrack",
    rstarts = "exonStarts", rends = "exonEnds", gene = "name",
    symbol = "name", transcript = "name", strand = "strand",
    fill = "#8282d2", name = "UCSC Genes"
)
refGenes <- UcscTrack(
    genome = "mm9", chromosome = "chrX",
    track = "xenoRefGene", from = from, to = to,
    trackType = "GeneRegionTrack", rstarts = "exonStarts",
    rends = "exonEnds", gene = "name", symbol = "name2",
    transcript = "name", strand = "strand", fill = "#8282d2",
    stacking = "dense", name = "Other RefSeq"
)
ensGenes <- UcscTrack(
    genome = "mm9", chromosome = "chrX",
    track = "ensGene", from = from, to = to, trackType = "GeneRegionTrack",
    rstarts = "exonStarts", rends = "exonEnds", gene = "name",
    symbol = "name2", transcript = "name", strand = "strand",
    fill = "#960000", name = "Ensembl Genes"
)
cpgIslands <- UcscTrack(
    genome = "mm9", chromosome = "chrX",
    track = "cpgIslandExt", from = from, to = to,
    trackType = "AnnotationTrack", start = "chromStart",
    end = "chromEnd", id = "name", shape = "box",
    fill = "#006400", name = "CpG Islands"
)
snpLocations <- UcscTrack(
    genome = "mm9", chromosome = "chrX",
    track = "snp128", from = from, to = to, trackType = "AnnotationTrack",
    start = "chromStart", end = "chromEnd", id = "name",
    feature = "func", strand = "strand", shape = "box",
    stacking = "dense", fill = "black", name = "SNPs"
)
conservation <- UcscTrack(
    genome = "mm9", chromosome = "chrX",
    track = "Conservation", table = "phyloP30wayPlacental",
    from = from, to = to, trackType = "DataTrack",
    start = "start", end = "end", data = "score",
    type = "hist", window = "auto", col.histogram = "darkblue",
    fill.histogram = "darkblue", ylim = c(-3.7, 4),
    name = "Conservation"
)
gcContent <- UcscTrack(
    genome = "mm9", chromosome = "chrX",
    track = "GC Percent", table = "gc5Base", from = from,
    to = to, trackType = "DataTrack", start = "start",
    end = "end", data = "score", type = "hist", window = -1,
    windowSize = 1500, fill.histogram = "black",
    col.histogram = "black", ylim = c(30, 70), name = "GC Percent"
)
axTrack <- GenomeAxisTrack()
idxTrack <- IdeogramTrack(genome = "mm9", chromosome = "chrX")

pdf("example.pdf")
plotTracks(list(
    idxTrack, axTrack, knownGenes, refGenes,
    ensGenes, cpgIslands, gcContent, conservation,
    snpLocations
), from = from, to = to, showTitle = FALSE)
dev.off()


#############
bwFile <- "/lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1/results/covmean/ChIP-seq_CGGBP1_ENCODE.coverage.bigWig"
dTrack2 <- DataTrack(
    range = bwFile, genome = "hg38",
    type = "l", name = "bigWig"
)

from <- 77022543
to <- 77030708
chr <- "chr14"
pdf("example.pdf")
plotTracks(dTrack2, chromosome = chr, from = from, to = to)
dev.off()

##############
library(Gviz)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(biomaRt)

## Outdir
outDir <- "/lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1/results/genome_tracks"
# Load gene annotations
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Load the BigWig file
bg_encode <- "/lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1/results/covmean/ChIP-seq_CGGBP1_ENCODE.coverage.bedgraph"
bg_encode_track <- DataTrack(range = bg_encode, type = "hist", genome = "hg38", name = "ChIP-seq CGGBP1", col.histogram = "darkblue", fill.histogram = "darkblue", ylim = c(0, 1))

bg_ctrl <- "/lustre/groups/ies/projects/hamperl_lab/elizabeth.marquezgom/Augusto/ChIPseq-CGGBP1/results/covmean/ChIP-seq_ctrl.coverage.bedgraph"
bg_ctrl_track <- DataTrack(range = bg_ctrl, type = "hist", genome = "hg38", name = "ChIP-seq CGGBP1 control", col.histogram = "darkblue", fill.histogram = "darkblue", ylim = c(0, 1))


PlotTrack <- function(gene) {
    # Select the gene of interest
    gene_symbol <- gene
    gene_info <- select(org.Hs.eg.db, keys = gene_symbol, columns = c("SYMBOL", "ENTREZID"), keytype = "SYMBOL")
    entrez_id <- gene_info$ENTREZID
    gene_range <- genes(txdb, filter = list(gene_id = entrez_id))

    # Create the gene region track
    gene_track <- GeneRegionTrack(txdb,
        chromosome = as.character(seqnames(gene_range)),
        start = start(gene_range) - 1000000, end = end(gene_range) + 1000000,
        geneSymbol = FALSE, name = "Gene region", strand = strand(gene_range), collapseTranscripts = TRUE,
        transcriptAnnotation = FALSE
    )

    # transcript_ids <- symbol(gene_track)
    # # Use biomaRt to get gene symbols corresponding to ENST IDs
    # ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    # gene_data <- getBM(attributes = c('ensembl_transcript_id', 'hgnc_symbol'),
    #                    filters = 'ensembl_transcript_id',
    #                    values = transcript_ids,
    #                    mart = ensembl)

    # # Match ENST IDs with gene symbols
    # mapped_symbols <- gene_data$hgnc_symbol[match(transcript_ids, gene_data$ensembl_transcript_id)]

    # # Update the track to show gene symbols instead of ENST IDs
    # symbol(gene_track) <- mapped_symbols

    # gene_symbols <- mapIds(org.Hs.eg.db, keys = transcript_ids, column = "SYMBOL", keytype = "ENSEMBLTRANS", multiVals = "first")

    # Load chromosome
    ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = gene_track@chromosome, from = from, to = to)
    axTrack <- GenomeAxisTrack()

    # Plot the tracks
    pdf(paste0(outDir, "/", gene, "_track.pdf"))
    plotTracks(list(ideoTrack, axTrack, gene_track, bg_encode_track, bg_ctrl_track),
        from = start(gene_range) - 1000000, to = end(gene_range) + 1000000,
        chromosome = as.character(seqnames(gene_range)), transcriptAnnotation = FALSE, showBandId = TRUE, cex.bands = 0.9,
        collapseTranscripts = TRUE
    )
    dev.off()
}

genes <- c("IRF2BPL", "RPS29", "C7orf50", "EGR1", "PBX1", "EIF3F", "SEC22B")
genes <- c("IRF2BPL")

lapply(genes, PlotTrack)
