




##################################################################################################################################
##################################################################################################################################




callback = function(hc, mat){
        sv = svd(t(mat))$v[,1]
        dend = (reorder(as.dendrogram(hc), wts = sv))
        as.hclust(dend)
}





##################################################################################################################################
##################################################################################################################################


setupDDS <- function (SampleTableName, CountTableName, SampleIdName, ConditionName, 
                      BatchName = NULL, n_samples_for_filtering = 3, min_number_of_reads = 1) 
{
        SampleTable <- get(SampleTableName)
        my_counts_genes <- get(CountTableName)
        stopifnot(identical(colnames(my_counts_genes), as.character(SampleTable[, 
                                                                                SampleIdName])))
        filter <- apply(my_counts_genes, 1, function(x) length(x[x > 
                                                                         min_number_of_reads]) >= n_samples_for_filtering)
        my_counts_filtered <- my_counts_genes[filter, ]
        my_colData <- DataFrame(Sample = factor(SampleTable[, ConditionName]))
        rownames(my_colData) <- SampleTable[, SampleIdName]
        if (is.null(BatchName)) {
                my_design = formula("~Sample")
        }
        else {
                my_design = formula("~Batch+Sample")
                my_colData$Batch <- factor(SampleTable[, BatchName])
        }
        dds <- DESeqDataSetFromMatrix(countData = my_counts_filtered, 
                                      colData = my_colData, design = my_design)
        dds <- DESeq(dds)
        return(dds)
}



##################################################################################################################################
##################################################################################################################################



getResults <- function (dds, contrast, result_name = NULL, lfc_cutoff = 0, 
                        shrink = FALSE) 
{
        if (is.null(result_name)) {
                res <- results(dds, contrast = c("Sample", contrast), 
                               lfcThreshold = lfc_cutoff, independentFiltering = FALSE)
                if (shrink) {
                        res <- lfcShrink(dds, contrast = c("Sample", contrast), 
                                         type = "normal", res = res)
                }
        }
        else {
                res <- results(dds, name = result_name, lfcThreshold = lfc_cutoff, 
                               independentFiltering = FALSE)
        }
        
        res$padj[is.na(res$padj)] <- 1
        res$log2FoldChange[is.na(res$log2FoldChange)] <- 0
        res <- res[order(res$pvalue), ]
        return(res)
}







##################################################################################################################################
##################################################################################################################################


plottingMA <- function (res, main_title = "", main_title_size = 1, point_size = 0.25, 
                        point_color = rgb(0.7, 0.7, 0.7, 0.5), sign_point_color = rgb(0.9, 0.6, 0, 0.5),
                        selection_ids = NULL, selection_name = "labeled", 
                        selection_id_type = "symbol", selection_point_size = 0.5, 
                        selection_point_color = rgb(0, 0, 0, 1), selection_sign_point_color = rgb(0.8, 0, 0, 1), 
                        selection_text_label = FALSE, selection_text_size = 1, 
                        selection_text_adj = -0.5, selection_shadow = FALSE, 
                        xlims = c(0, 6), ylims = c(-5, 5), x_axis_by = 2, padj_cutoff = 0.01, 
                        show_legend = TRUE, legend_pos_label = "topright", legend_pos_all = "bottomright") 
{
        res$log10baseMean <- log10(res$baseMean + 1)
        scattermoreplot(x = res$log10baseMean, y = res$log2FoldChange, xlab = "log10 mean counts", 
             ylab = "log2 fold change", xlim = xlims, ylim = ylims, 
             col = point_color, pch = 19, cex = point_size, xaxt = "n",
             size = c(1024, 1024))
        axis(side = 1, at = seq(from = xlims[1], to = xlims[2], 
                                by = x_axis_by))
        abline(h = 0, col = "grey32")
        res.sign <- res[res$padj < padj_cutoff, ]
        points(x = res.sign$log10baseMean, y = res.sign$log2FoldChange, 
               col = sign_point_color, pch = 19, cex = selection_point_size)
        mtext(text = main_title, side = 3, line = 0.5, adj = 0.5, 
              font = 2, cex = main_title_size)
        if (!(is.null(selection_ids))) {
                selection_vector <- res[selection_id_type][, 1] %in% 
                        selection_ids
                selection_color <- ifelse(res$padj < padj_cutoff, selection_sign_point_color, 
                                          selection_point_color)
                points(x = res$log10baseMean[selection_vector], y = res$log2FoldChange[selection_vector], 
                       col = selection_color[selection_vector], pch = 16, 
                       cex = selection_point_size)
                if (selection_shadow) {
                        points(x = res$log10baseMean[selection_vector], 
                               y = res$log2FoldChange[selection_vector], col = "black", 
                               pch = 1, lwd = 0.75, cex = selection_point_size)
                }
                if (selection_text_label) {
                        if (selection_shadow) {
                                text(x = res$log10baseMean[selection_vector], 
                                     y = res$log2FoldChange[selection_vector], 
                                     labels = res[selection_id_type][, 1][selection_vector], 
                                     col = "black", adj = c(0, selection_text_adj), 
                                     font = 2, cex = selection_text_size)
                        }
                        text(x = res$log10baseMean[selection_vector], y = res$log2FoldChange[selection_vector], 
                             labels = res[selection_id_type][, 1][selection_vector], 
                             col = selection_color[selection_vector], adj = c(0, 
                                                                              selection_text_adj), cex = selection_text_size)
                }
        }
        if (show_legend) {
                if (!(is.null(selection_ids)) & sum(selection_vector) > 
                    0) {
                        legend(legend_pos_label, legend = paste(selection_name, 
                                                                c("significant", "non-significant")), col = c(selection_sign_point_color, 
                                                                                                              selection_point_color), bg = "white", border = NA, 
                               bty = "n", cex = 0.8, pch = 19)
                }
                legend(legend_pos_all, legend = c("all significant", 
                                                  "all non-significant"), col = c(sign_point_color, 
                                                                                  point_color), bg = "white", border = NA, bty = "n", 
                       cex = 0.8, pch = 19)
        }
}







##################################################################################################################################
##################################################################################################################################


plottingPCA <- function (my_data, xcomp = 1, ycomp = 2, conditions, pca_colors = c("#999999", 
                                                                    "#0072B2", "#CC79A7", "#009E73", "#E69F00", "#D55E00", "#56B4E9", 
                                                                    "#F0E442"), main_title = "PCA", quantiles = c(0, 1), show_labels = TRUE, 
          point_size = 1.1, my_xlimits = c(-100, 100), my_ylimits = c(-100, 
                                                                      100)) 
{
        my_data <- na.omit(my_data)
        rv <- rowVars(my_data)
        selection <- (rv > quantile(rv, quantiles[1]) & rv < quantile(rv, 
                                                                      quantiles[2]))
        pca <- prcomp(t(my_data[selection, ]), scale. = TRUE)
        percentVar <- round(pca$sdev^2/sum(pca$sdev^2) * 100, 1)[1:10]
        plot(pca$x[, xcomp], pca$x[, ycomp] * -1, col = pca_colors[conditions], 
             pch = 16, cex = point_size, xlab = paste("PC", xcomp, 
                                                      " (", percentVar[xcomp], "%)", sep = ""), ylab = paste("PC", 
                                                                                                             ycomp, " (", percentVar[ycomp], "%)", sep = ""), 
             xlim = my_xlimits, ylim = my_ylimits)
        points(pca$x[, xcomp], pca$x[, ycomp] * -1, cex = point_size, 
               col = "#555555", pch = 1, lwd = 0.5)
        mtext(text = main_title, side = 3, line = 0.5, adj = 0.5, 
              font = 2)
        if (show_labels) {
                text(pca$x[, xcomp], pca$x[, ycomp] * -1, labels = rownames(pca$x), 
                     adj = -0.5, col = "gray32", cex = 0.5)
        }
}



##################################################################################################################################
##################################################################################################################################













##################################################################################################################################
##################################################################################################################################

