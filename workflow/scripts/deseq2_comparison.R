






##################################################################################################################################
##################################################################################################################################




library(DESeq2)
library(ShortRead)
library(rtracklayer)
library(sva)


library(RColorBrewer)
library(pheatmap)
library(scattermore)

library(parallel)

source("workflow/scripts/functions.R")



##################################################################################################################################
##################################################################################################################################





args = commandArgs(trailingOnly=TRUE)


input_table <- grep("Design", args, value = TRUE)
input_counts <- grep("counts_per_bin", args, value = TRUE)


output_sessionInfo <- grep("deseq2_comparison/sessionInfo", args, value = TRUE)


output_tables <- gsub("sessionInfo.txt", "tables",output_sessionInfo)
output_plots  <- gsub("sessionInfo.txt", "plots",output_sessionInfo)


dir.create(output_tables, showWarnings = FALSE, recursive = TRUE)
dir.create(output_plots,  showWarnings = FALSE, recursive = TRUE)



##################################################################################################################################
##################################################################################################################################

######################################################   Annotation    ##########################################################



my_sites <- import("results/deseq2_insertions/DMSO/final_sites.bed")



##################################################################################################################################
##################################################################################################################################













##################################################################################################################################
##################################################################################################################################

#################################################        SampleTables        #####################################################




SampleTable1 <- read.csv(input_table, header = T, stringsAsFactors = F)
SampleTable1$SampleID <- gsub("_S.*","",SampleTable1$fastq_1)
SampleTable1$Batch <- paste0("r",SampleTable1$replicate)
SampleTable1$antibody[SampleTable1$antibody == ""] <- "Input"


SampleTable2 <- read.csv("../202302_chip_seq_h3/data/Design_file.csv", header = T, stringsAsFactors = F)
SampleTable2$SampleID <- paste(SampleTable2$group, SampleTable2$replicate, sep = "_")
SampleTable2$Batch <- paste0("r",SampleTable2$replicate)
SampleTable2$antibody[SampleTable2$antibody == ""] <- "Input"



##################################################################################################################################
##################################################################################################################################




my_color_palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442")

my_color_palette2 <- brewer.pal(9,"Set1")



##################################################################################################################################
##################################################################################################################################









##################################################################################################################################
##################################################################################################################################

#################################################          tables            #####################################################


################################################################################

lnc_ser2p_DMSO  <- read.delim("results/deseq2_genomewide/DMSO/tables/lnc.ser2p.txt", row.names = 1, check.names = FALSE)
lnc_ser2p_VE821 <- read.delim("results/deseq2_genomewide/VE821/tables/lnc.ser2p.txt", row.names = 1, check.names = FALSE)

lnc_ser2p <- merge(lnc_ser2p_DMSO, lnc_ser2p_VE821, by="row.names")
rownames(lnc_ser2p) <- lnc_ser2p$Row.names
lnc_ser2p <- lnc_ser2p[,-1]

SampleTable1 <- SampleTable1[match(colnames(lnc_ser2p), SampleTable1$SampleID),]
stopifnot(identical(SampleTable1$SampleID, colnames(lnc_ser2p)))

lnc_ser2p_means <- sapply(unique(SampleTable1$group), function(x){ 
        
        rowMeans(lnc_ser2p[,SampleTable1$SampleID[SampleTable1$group == x]])  
})

colnames(lnc_ser2p_means) <- paste("ser2p", colnames(lnc_ser2p_means), sep = " ")

################################################################################




################################################################################

lnc_h3_DMSO  <- read.delim("../202302_chip_seq_h3/results/deseq2_genomewide/DMSO/tables/lnc.h3.txt", row.names = 1, check.names = FALSE)
lnc_h3_VE821 <- read.delim("../202302_chip_seq_h3/results/deseq2_genomewide/VE821/tables/lnc.h3.txt", row.names = 1, check.names = FALSE)

lnc_h3 <- merge(lnc_h3_DMSO, lnc_h3_VE821, by="row.names")
rownames(lnc_h3) <- lnc_h3$Row.names
lnc_h3 <- lnc_h3[,-1]

SampleTable2 <- SampleTable2[match(colnames(lnc_h3), SampleTable2$SampleID),]
stopifnot(identical(SampleTable2$SampleID, colnames(lnc_h3)))

lnc_h3_means <- sapply(unique(SampleTable2$group), function(x){ 
        
        rowMeans(lnc_h3[,SampleTable2$SampleID[SampleTable2$group == x]])  
})

colnames(lnc_h3_means) <- paste("h3", colnames(lnc_h3_means), sep = " ")

################################################################################




##################################################################################################################################
##################################################################################################################################


lnc_merged <- merge(lnc_h3_means, lnc_ser2p_means, by="row.names")
rownames(lnc_merged) <- lnc_merged$Row.names
lnc_merged <- lnc_merged[,-1]

lnc_merged <- scale(lnc_merged, center = TRUE, scale = TRUE)

my_bins <- data.frame(chr = gsub("_.*","",rownames(lnc_merged)),
                      start = as.numeric(gsub(".*_","",rownames(lnc_merged))),
                      end   = as.numeric(gsub(".*_","",rownames(lnc_merged)))+1e3,
                      row.names = rownames(lnc_merged))

gr_bins <- makeGRangesFromDataFrame(my_bins)


gr_bins$distance <- NA
my_distances <- distanceToNearest(gr_bins, my_sites)
gr_bins$distance[queryHits(my_distances)] <- mcols(my_distances)$distance


Less3kb <- names(subsetByOverlaps(gr_bins, my_sites, maxgap = 3e3))
Sites <- names(subsetByOverlaps(gr_bins, my_sites))
mAIRN <- grep("mAIRN", names(gr_bins), value = TRUE)


##################################################################################################################################
##################################################################################################################################



pdf(paste0(output_plots,"/comparisons.pdf"), height = 6, width = 6, useDingbats = F)

par(mfrow=c(1,1), mar = c(4,4,2,2), oma = c(4,4,4,4), mgp = c(2.5,1,0))

for(idx in colnames(lnc_merged)){
        
        for(idy in colnames(lnc_merged)){
                
                if(idx == idy){next()}
                
                scattermoreplot(x = lnc_merged[,idx],
                                y = lnc_merged[,idy],
                                xlab = paste(paste(rev(strsplit(idx, "\\.")[[1]]), collapse = " "), "(z-score)"),
                                ylab = paste(paste(rev(strsplit(idy, "\\.")[[1]]), collapse = " "), "(z-score)"),
                                xlim = c(-5,12),
                                ylim = c(-5,12),
                                col = rgb(0.7, 0.7, 0.7, 0.1),
                                pch = 19,
                                cex = 3,
                                size = c(1024, 1024))
                
                points(x = lnc_merged[rownames(lnc_merged) %in% Less3kb,idx],
                       y = lnc_merged[rownames(lnc_merged) %in% Less3kb,idy],
                       col = "black", pch = 19, cex = 0.5)
                
                points(x = lnc_merged[rownames(lnc_merged) %in% Sites,idx],
                       y = lnc_merged[rownames(lnc_merged) %in% Sites,idy],
                       col = rgb(0.9, 0.6, 0, 1), pch = 19, cex = 0.5)
                
                points(x = lnc_merged[rownames(lnc_merged) %in% mAIRN,idx],
                       y = lnc_merged[rownames(lnc_merged) %in% mAIRN,idy],
                       col = rgb(0.8, 0, 0, 1), pch = 19, cex = 0.5)
                
                # smoothScatter(x = lnc_merged[,idx], 
                #               y = lnc_merged[,idy],  
                #               xlab = paste(rev(strsplit(idx, "\\.")[[1]]), collapse = " "), 
                #               ylab = paste(rev(strsplit(idy, "\\.")[[1]]), collapse = " "), 
                #               xlim = c(0,8), 
                #               ylim = c(0,8), 
                #               colramp = colorRampPalette(c("white", blues9, "black")),
                #               cex = 0, useRaster=TRUE)
                
                abline(coef = c(0,1), col = "grey32")
                
                legend("bottomright", 
                       legend = c(paste("all n =", nrow(lnc_merged)),
                                  paste("<3kb n =", sum(rownames(lnc_merged) %in% Less3kb)),
                                  paste("sites n =", sum(rownames(lnc_merged) %in% Sites)),
                                  paste("mAIRN n =", sum(rownames(lnc_merged) %in% mAIRN))),
                       col = c(rgb(0.7, 0.7, 0.7, 0.5), "black",rgb(0.9, 0.6, 0, 1),rgb(0.8, 0, 0, 1)),
                       pch = 19, cex = 0.8)
        }
}



dev.off()









##################################################################################################################################
##################################################################################################################################



writeLines(capture.output(sessionInfo()), output_sessionInfo)




##################################################################################################################################
##################################################################################################################################


