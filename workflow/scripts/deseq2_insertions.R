






##################################################################################################################################
##################################################################################################################################




library(DESeq2)
library(ShortRead)
library(rtracklayer)
library(sva)


library(RColorBrewer)
library(pheatmap)
library(scattermore)


source("workflow/scripts/functions.R")



##################################################################################################################################
##################################################################################################################################





args = commandArgs(trailingOnly=TRUE)


input_table <- grep("Design", args, value = TRUE)
input_counts <- grep("counts_per_bin", args, value = TRUE)


output_sessionInfo <- grep("session", args, value = TRUE)


output_tables <- gsub("sessionInfo.txt", "tables",output_sessionInfo)
output_plots  <- gsub("sessionInfo.txt", "plots",output_sessionInfo)


dir.create(output_tables, showWarnings = FALSE, recursive = TRUE)
dir.create(output_plots,  showWarnings = FALSE, recursive = TRUE)



##################################################################################################################################
##################################################################################################################################

######################################################   Annotation    ##########################################################


my_site_files <- list.files(path = "data/external/mapped_sites/", pattern = "flt.bed", recursive = TRUE, full.names = TRUE)

i=1

gr <- GRanges()

for(i in seq_along(my_site_files)){
        
        gr_tmp <- import(my_site_files[i])
        gr <- c(gr, gr_tmp)
        
        rm(list = "gr_tmp")
}


my_sites <- reduce(resize(gr, width = 200, fix = "center"))

seqlevels(my_sites) <- seqlevels(my_sites)[order(as.numeric(seqlevels(my_sites)))]
my_sites <- sort(my_sites)

my_sites$n_sites <- countOverlaps(my_sites, gr)
my_sites <- my_sites[my_sites$n_sites >= 2]

seqlevelsStyle(my_sites) <- "Ensembl"

rm(list = "gr")


export.bed(my_sites, con = paste0(gsub("/tables","",output_tables),"/final_sites.bed"))


##################################################################################################################################
##################################################################################################################################













##################################################################################################################################
##################################################################################################################################

#################################################        SampleTable         #####################################################




SampleTable <- read.csv(input_table, header = T, stringsAsFactors = F)
SampleTable$SampleID <- gsub("_S.*","",SampleTable$fastq_1)

SampleTable$Batch <- paste0("r",SampleTable$replicate)
SampleTable$antibody[SampleTable$antibody == ""] <- "Input"

#SampleTable <- SampleTable[!duplicated(SampleTable$SampleID), , drop = FALSE]
#SampleTable <- SampleTable[order(SampleTable$SampleID), , drop = FALSE]


my_sampleset <- gsub(".*/","",gsub("/session.*","",output_sessionInfo))
SampleTable <- SampleTable[grep(my_sampleset, SampleTable$group),]


##################################################################################################################################
##################################################################################################################################




my_color_palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442")

my_color_palette2 <- brewer.pal(9,"Set1")



##################################################################################################################################
##################################################################################################################################









##################################################################################################################################
##################################################################################################################################

#################################################        count table         #####################################################




my_count_table <- read.delim(input_counts, check.names = FALSE)
rownames(my_count_table) <- paste(my_count_table[,1],my_count_table[,2] , sep = "_")

my_count_table <- my_count_table[!grepl("^KI|^MT|^GL|mAIRN_5000", rownames(my_count_table)),]


my_bins <- my_count_table[,1:3]
colnames(my_bins) <- c("chr","start","end")

gr_bins <- makeGRangesFromDataFrame(my_bins)


gr_bins$distance <- NA
my_distances <- distanceToNearest(gr_bins, my_sites)
gr_bins$distance[queryHits(my_distances)] <- mcols(my_distances)$distance



my_count_table <- my_count_table[,-c(1:3)]
colnames(my_count_table) <- gsub("_S.*|\\'","",colnames(my_count_table))






##################################################################################################################################
##################################################################################################################################












##################################################################################################################################
##################################################################################################################################

#################################################           Setup            #####################################################






for(my_assay in c(unique(SampleTable$antibody))){
        
        
        #########################################################
        
        my_count_tmp <- my_count_table
        
        
        SampleTable_tmp  <- SampleTable[SampleTable$antibody == my_assay,]
        my_subset <- colnames(my_count_tmp) %in% SampleTable_tmp$SampleID
        
        
        if(sum(my_subset) <= 4){next()}
        
        my_count_sub_tmp <- my_count_tmp[,my_subset]
        SampleTable_tmp  <- SampleTable_tmp[match(SampleTable_tmp$SampleID, colnames(my_count_sub_tmp)),]
        
        my_region_names <- names(gr_bins)[gr_bins$distance < 1.5e4 & !is.na(gr_bins$distance)]
        my_region_names <- c(my_region_names, grep("mAIRN", names(gr_bins), value = TRUE))
        
        my_count_sub_tmp <- my_count_sub_tmp[rownames(my_count_sub_tmp) %in% my_region_names,]
        
        #########################################################
        
        if(all(identical(colnames(my_count_sub_tmp), SampleTable_tmp$SampleID),
               identical(colnames(my_count_sub_tmp), colnames(my_count_tmp[,my_subset])))){
                
                
                dds_tmp <- setupDDS(CountTableName = "my_count_sub_tmp",
                                    SampleTableName = "SampleTable_tmp",
                                    SampleIdName = "SampleID",
                                    ConditionName = "group",
                                    BatchName = "Batch",
                                    n_samples_for_filtering = (ncol(my_count_sub_tmp)/2),
                                    min_number_of_reads = 4)
                
                sizeFactors(dds_tmp) <- estimateSizeFactorsForMatrix(my_count_tmp[,my_subset])
                
                dds_tmp <- DESeq(dds_tmp)
                
                assign(paste("dds", my_assay, sep = "."), dds_tmp)
                
        }
        
        #########################################################
        
        rm(list = ls(pattern = "_tmp"))
}


##################################################################################################################################
##################################################################################################################################












##################################################################################################################################
##################################################################################################################################

#################################################           contrasts         #####################################################



my_assays <- gsub("dds\\.","",ls(pattern = "^dds\\."))


if(length(my_assays) > 0){
        
        for(my_assay in my_assays){
                
                
                dds_tmp <- get(paste("dds", my_assay, sep = "."))
                
                contrast_list <- list(c(levels(colData(dds_tmp)$Sample)[1], levels(colData(dds_tmp)$Sample)[2]))
                
                for(contrast_name in contrast_list){
                        
                        res <- getResults(dds = dds_tmp,
                                          contrast = contrast_name,
                                          lfc_cutoff = 0,
                                          shrink = FALSE)
                        
                        res$site_id <- rownames(res)
                        res$overlap <- res$site_id %in% names(subsetByOverlaps(gr_bins, my_sites))
                        
                        res$distance <- gr_bins$distance[match(res$site_id, names(gr_bins))]
                        
                        res_name <- paste("res", my_assay,
                                          paste0(contrast_name[1],"-",contrast_name[2]), sep = ".")
                        
                        assign(res_name, res)
                        
                        write.table(res, file = paste0(output_tables,"/",res_name, ".txt"),
                                    quote = F, sep = "\t", row.names = T, col.names = NA)
                        
                        # write.table(head(res, 100), file = paste0(output_tables,"/",res_name, ".top100.txt"),
                        #             quote = F, sep = "\t", row.names = T, col.names = NA)
                        
                        rm(list = "res")
                }
                
                rm(list = ls(pattern = "_tmp"))
        }
        
}



##################################################################################################################################
##################################################################################################################################











##################################################################################################################################
##################################################################################################################################

#################################################           MA-plots         #####################################################



res_names <- ls(pattern = "^res\\.")


Sites <- names(subsetByOverlaps(gr_bins, my_sites))
Less3kb <- names(subsetByOverlaps(gr_bins, my_sites, maxgap = 3e3))
mAIRN <- grep("mAIRN", names(gr_bins), value = TRUE)

padj_cutoff <- 0.1

NoLabel <- ""

for(my_label in c("NoLabel","Sites","Less3kb","mAIRN")){
        
        if(length(res_names) > 0){
                
                for(res_name in res_names){
                        
                        pdf(paste0(output_plots,"/MA.",gsub("res.","",res_name),"_",my_label,".pdf"), height = 6, width = 6, useDingbats = F)
                        
                        par(mfrow=c(1,1), mar = c(4,4,2,2), oma = c(4,4,4,4), mgp = c(2.5,1,0))
                        
                        plottingMA(res = get(res_name),
                                   main_title = gsub(".*\\.","",res_name),
                                   main_title_size = 1.5,
                                   point_size = 3,
                                   selection_ids = get(my_label),
                                   selection_id_type = "site_id",
                                   selection_name = my_label,
                                   selection_point_size = 0.5,
                                   selection_text_label = FALSE,
                                   selection_shadow = FALSE,
                                   xlims = c(0, 2.5),
                                   ylims = c(-6,  6),
                                   x_axis_by = 0.5,
                                   padj_cutoff = padj_cutoff,
                                   show_legend = TRUE)
                        
                        dev.off()
                }
        }
}

##################################################################################################################################
##################################################################################################################################












##################################################################################################################################
##################################################################################################################################

#################################################             rld            #####################################################



my_assays <- gsub("dds\\.","",ls(pattern = "^dds\\."))


for(my_assay in my_assays){
        
        
        
        dds_tmp <- get(paste("dds", my_assay, sep = "."))
        
        modcombat <- model.matrix(~Sample, data = colData(dds_tmp))
        batchVar <- colData(dds_tmp)$Batch
        
        #########################################################
        
        
        
        lnc_tmp <- log2(counts(dds_tmp, normalized=TRUE)+1)
        
        
        
        if(!is.null(batchVar)){
                
                lnc_tmp <- ComBat(dat = lnc_tmp,
                                  batch = batchVar, mod = modcombat,
                                  par.prior = TRUE, prior.plots = FALSE) 
        }
        
        assign(paste("lnc", my_assay, sep = "."), lnc_tmp)
        
        
        
        write.table(lnc_tmp, file = paste0(output_tables,"/lnc.", my_assay,".txt"),
                    quote = F, sep = "\t", row.names = T, col.names = NA)
        
        #########################################################
        
        rm(list = ls(pattern = "_tmp"))
}



##################################################################################################################################
##################################################################################################################################













##################################################################################################################################
##################################################################################################################################

#################################################        Correlation         #####################################################



callback = function(hc, mat){
        sv = svd(t(mat))$v[,1]
        dend = (reorder(as.dendrogram(hc), wts = sv))
        as.hclust(dend)
}



for(my_assay in my_assays){
        
        
        pdf(paste0(output_plots,"/correlation.", my_assay,".pdf"), width = 10, height = 10, useDingbats = FALSE)
        par(oma=c(2,2,2,0), mar=c(4,4,4,4), mgp=c(1.5,0.75,0),cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5, pch=19)
        
        
        lnc_tmp <- get(paste("lnc", my_assay, sep = "."))
        
        lnc_scaled_tmp <- lnc_tmp #- rowMeans(lnc_tmp)
        
        my_corr <- cor(lnc_scaled_tmp, method = "spearman")
        rownames(my_corr) <- SampleTable$group[match(colnames(my_corr), SampleTable$SampleID)]
        
        
        pheatmap(my_corr,
                 main = paste0("Spearman´s R [n = ", nrow(lnc_scaled_tmp),"]"),
                 color = colorRampPalette((brewer.pal(9,name = "YlGnBu")))(100),
                 breaks = seq(0.5,1, length.out = 101),
                 clustering_callback = callback,
                 clustering_method = "ward.D2", 
                 cellheight = 12, 
                 cellwidth = 12)
        
        
        rm(list = ls(pattern = "_tmp"))
        rm(list = "my_corr")
        
        dev.off()
        
}



##################################################################################################################################
##################################################################################################################################









##################################################################################################################################
##################################################################################################################################

#################################################            PCA             #####################################################


my_assays <- gsub("dds\\.","",ls(pattern = "^dds\\."))


xycomps <- list(c(1,2),
                c(1,3),
                c(2,3))


for(my_assay in my_assays){
        
        pdf(paste0(output_plots,"/PCA.",my_assay,".pdf"), height = 6.25, width = 9, useDingbats = F)
        par(mfrow=c(2,3), mar = c(4,4,2,2), oma = c(3,3,3,3), mgp = c(2,1,0))
        
        
        #########################################################        
        
        
        lnc_tmp <- get(paste("lnc", my_assay, sep = "."))
        
        my_conditions <- SampleTable[,"group"][match(colnames(lnc_tmp), SampleTable$SampleID)]
        my_conditions <- factor(my_conditions, levels = unique(my_conditions))
        
        
        for(xycomp in xycomps){
                
                plottingPCA(lnc_tmp,
                            xcomp = xycomp[1],
                            ycomp = xycomp[2],
                            conditions = my_conditions,
                            pca_colors = my_color_palette,
                            main_title = paste("PCA -", my_assay),
                            quantiles = c(0,1),
                            show_labels = TRUE,
                            point_size = 1.1,
                            my_xlimits = c(-30,30),
                            my_ylimits = c(-30,30))
        }
        
        plot.new()
        
        legend("top",
               legend = levels(my_conditions),
               horiz = FALSE, cex = 0.8,
               fill =  my_color_palette[seq_along(levels(my_conditions))])
        
        #########################################################
        
        
        rm(list = ls(pattern = "_tmp"))
        
        
        dev.off()
        
}

##################################################################################################################################
##################################################################################################################################













##################################################################################################################################
##################################################################################################################################

#################################################          Heatmap           #####################################################


my_assays <- gsub("dds\\.","",ls(pattern = "^dds\\."))


for(my_assay in my_assays){
        
        pdf(paste0(output_plots,"/Heatmap.",my_assay,".pdf"), width = 8, height = 11, useDingbats = FALSE)
        par(oma=c(2,2,2,0), mar=c(4,4,4,4), mgp=c(1.5,0.75,0),cex.axis = 1.5, cex.main = 1.5, cex.lab=1.5, pch=19)
        
        
        ################################################################################
        
        
        lnc_tmp <- get(paste("lnc", my_assay, sep = "."))
        
        my_conditions <- SampleTable[,"group"][match(colnames(lnc_tmp), SampleTable$SampleID)]
        my_conditions <- factor(my_conditions, levels = unique(my_conditions))
        
        my_hm_data <- lnc_tmp
        my_hm_data <- my_hm_data - rowMeans(my_hm_data)
        
        pheatmap(my_hm_data, 
                 main = "relative log2 counts",
                 # annotation_col = annotation_col,
                 # annotation_colors = ann_colors,
                 color = colorRampPalette(rev(brewer.pal(9,name = "RdBu")))(100), 
                 breaks = seq(-2,2, length.out = 101),
                 clustering_callback = callback,
                 clustering_method = "ward.D2", 
                 cellwidth = 8, cellheight = 5, 
                 fontsize_col = 6, 
                 fontsize_row = 3)
        
        ################################################################################
        
        rm(list = "my_hm_data")
        
        dev.off()
        
        ################################################################################
}


##################################################################################################################################
##################################################################################################################################









##################################################################################################################################
##################################################################################################################################



writeLines(capture.output(sessionInfo()), output_sessionInfo)




##################################################################################################################################
##################################################################################################################################


