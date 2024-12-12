> GOenrichment <- function(gene_names) {
    ## pval default is 0.05, relaxed to 0.1 to get+     ## pval default is 0.05, relaxed to 0.1 to get genes in nascent set
+     enrichGO(gene_names, OrgDb = org.Hs.eg.db, ont = mode, keyType = "ENSEMBL", pAdjustMethod = "fdr", pvalueCutoff = 0.1)
+ }

> GOenrichment(wt$ensembl)
#
# over-representation test
#
#...@organism    Homo sapiens 
#...@ontology    BP 
#...@keytype     ENSEMBL 
#...@gene        chr [1:220] "ENSG00000125378" "ENSG00000105974" "ENSG00000184009" ...
#...pvalues adjusted by 'fdr' with cutoff <0.1 
#...210 enriched terms found
'data.frame':   210 obs. of  9 variables:
 $ ID         : chr  "GO:1901990" "GO:1901991" "GO:0000082" "GO:0045930" ...
 $ Description: chr  "regulation of mitotic cell cycle phase transition" "negative regulation of mitotic cell cycle phase transition" "G1/S transition of mitotic cell cycle" "negative regulation of mitotic cell cycle" ...
 $ GeneRatio  : chr  "18/205" "13/205" "14/205" "14/205" ...
 $ BgRatio    : chr  "391/21261" "214/21261" "267/21261" "289/21261" ...
 $ pvalue     : num  5.05e-08 1.70e-07 3.44e-07 8.91e-07 1.39e-06 ...
 $ p.adjust   : num  0.000168 0.000283 0.000381 0.00074 0.000921 ...
 $ qvalue     : num  0.000142 0.000239 0.000322 0.000625 0.000778 ...
 $ geneID     : chr  "ENSG00000110092/ENSG00000089685/ENSG00000122641/ENSG00000158290/ENSG00000164109/ENSG00000168411/ENSG00000197256"| __truncated__ "ENSG00000110092/ENSG00000089685/ENSG00000122641/ENSG00000164109/ENSG00000168411/ENSG00000197256/ENSG00000175215"| __truncated__ "ENSG00000205250/ENSG00000110092/ENSG00000180198/ENSG00000122641/ENSG00000158290/ENSG00000168411/ENSG00000197256"| __truncated__ "ENSG00000125378/ENSG00000110092/ENSG00000089685/ENSG00000122641/ENSG00000164109/ENSG00000168411/ENSG00000197256"| __truncated__ ...
 $ Count      : int  18 13 14 14 14 8 13 8 15 13 ...
#...Citation
 T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu.
 clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.
 The Innovation. 2021, 2(3):100141 
 
> GOenrichment(mut$ensembl)
#
# over-representation test
#
#...@organism    Homo sapiens 
#...@ontology    BP 
#...@keytype     ENSEMBL 
#...@gene        chr [1:369] "ENSG00000138246" "ENSG00000132600" "ENSG00000131238" ...
#...pvalues adjusted by 'fdr' with cutoff <0.1 
#...4 enriched terms found
'data.frame':   4 obs. of  9 variables:
 $ ID         : chr  "GO:0034470" "GO:0006364" "GO:0016072" "GO:0006888"
 $ Description: chr  "ncRNA processing" "rRNA processing" "rRNA metabolic process" "endoplasmic reticulum to Golgi vesicle-mediated transport"
 $ GeneRatio  : chr  "21/342" "14/342" "15/342" "10/342"
 $ BgRatio    : chr  "490/21261" "250/21261" "292/21261" "144/21261"
 $ pvalue     : num  5.01e-05 5.79e-05 8.44e-05 1.16e-04
 $ p.adjust   : num  0.0938 0.0938 0.0938 0.0969
 $ qvalue     : num  0.0891 0.0891 0.0891 0.092
 $ geneID     : chr  "ENSG00000148688/ENSG00000171103/ENSG00000156697/ENSG00000105171/ENSG00000113360/ENSG00000130305/ENSG00000165526"| __truncated__ "ENSG00000148688/ENSG00000171103/ENSG00000156697/ENSG00000105171/ENSG00000113360/ENSG00000130305/ENSG00000165526"| __truncated__ "ENSG00000148688/ENSG00000171103/ENSG00000158796/ENSG00000156697/ENSG00000105171/ENSG00000113360/ENSG00000130305"| __truncated__ "ENSG00000107862/ENSG00000114988/ENSG00000265808/ENSG00000184840/ENSG00000182158/ENSG00000152700/ENSG00000151532"| __truncated__
 $ Count      : int  21 14 15 10
#...Citation
 T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu.
 clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.
 The Innovation. 2021, 2(3):100141 

> GOenrichment <- function(gene_names) {
+     ## pval default is 0.05, relaxed to 0.1 to get genes in nascent set
+     enrichGO(gene_names, OrgDb = org.Hs.eg.db, ont = mode, keyType = "ENSEMBL", pAdjustMethod = "fdr", pvalueCutoff = 0.05)
+ }
> GOenrichment(mut$ensembl)
#
# over-representation test
#
#...@organism    Homo sapiens 
#...@ontology    BP 
#...@keytype     ENSEMBL 
#...@gene        chr [1:369] "ENSG00000138246" "ENSG00000132600" "ENSG00000131238" ...
#...pvalues adjusted by 'fdr' with cutoff <0.05 
#...0 enriched terms found
#...Citation
 T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu.
 clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.
 The Innovation. 2021, 2(3):100141 









 GOenrichment <- function(gene_names) {
    ## Adjust the p-values for multiple testing correction. 
    ## This is essential because many enrichment tests are performed simultaneously, one for each GO term, and multiple testing increases the likelihood of false positives.
    ## pval default is 0.05, relaxed to 0.1 to get genes in nascent set
    enrichGO(gene_names, OrgDb = org.Hs.eg.db, ont = mode, keyType = "ENSEMBL", pAdjustMethod = "fdr", pvalueCutoff = 0.05)
}

> GOenrichment(wt$ensembl)
#
# over-representation test
#
#...@organism    Homo sapiens 
#...@ontology    BP 
#...@keytype     ENSEMBL 
#...@gene        chr [1:183] "ENSG00000131018" "ENSG00000120051" "ENSG00000211448" ...
#...pvalues adjusted by 'fdr' with cutoff <0.05 
#...0 enriched terms found
#...Citation
 T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu.
 clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.
 The Innovation. 2021, 2(3):100141 

 > GOenrichment(mut$ensembl)
#
# over-representation test
#
#...@organism    Homo sapiens 
#...@ontology    BP 
#...@keytype     ENSEMBL 
#...@gene        chr [1:190] "ENSG00000070404" "ENSG00000243819" "ENSG00000114686" ...
#...pvalues adjusted by 'fdr' with cutoff <0.05 
#...1 enriched terms found
'data.frame':   1 obs. of  9 variables:
 $ ID         : chr "GO:0010212"
 $ Description: chr "response to ionizing radiation"
 $ GeneRatio  : chr "9/171"
 $ BgRatio    : chr "145/21261"
 $ pvalue     : num 2.71e-06
 $ p.adjust   : num 0.00697
 $ qvalue     : num 0.00637
 $ geneID     : chr "ENSG00000166311/ENSG00000186104/ENSG00000175054/ENSG00000158161/ENSG00000004487/ENSG00000151876/ENSG00000067369"| __truncated__
 $ Count      : int 9
#...Citation
 T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu.
 clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.
 The Innovation. 2021, 2(3):100141 

 > GOenrichment <- function(gene_names) {
+     ## Adjust the p-values for multiple testing correction. 
+     ## This is essential because many enrichment tests are performed simultaneously, one for each GO term, and multiple testing increases the likelihood of false positives.
+     ## pval default is 0.05, relaxed to 0.1 to get genes in nascent set
+     enrichGO(gene_names, OrgDb = org.Hs.eg.db, ont = mode, keyType = "ENSEMBL", pAdjustMethod = "fdr", pvalueCutoff = 0.1)
+ }
> GOenrichment(wt$ensembl)
#
# over-representation test
#
#...@organism    Homo sapiens 
#...@ontology    BP 
#...@keytype     ENSEMBL 
#...@gene        chr [1:183] "ENSG00000131018" "ENSG00000120051" "ENSG00000211448" ...
#...pvalues adjusted by 'fdr' with cutoff <0.1 
#...0 enriched terms found
#...Citation
 T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu.
 clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.
 The Innovation. 2021, 2(3):100141 

> GOenrichment(mut$ensembl)
#
# over-representation test
#
#...@organism    Homo sapiens 
#...@ontology    BP 
#...@keytype     ENSEMBL 
#...@gene        chr [1:190] "ENSG00000070404" "ENSG00000243819" "ENSG00000114686" ...
#...pvalues adjusted by 'fdr' with cutoff <0.1 
#...1 enriched terms found
'data.frame':   1 obs. of  9 variables:
 $ ID         : chr "GO:0010212"
 $ Description: chr "response to ionizing radiation"
 $ GeneRatio  : chr "9/171"
 $ BgRatio    : chr "145/21261"
 $ pvalue     : num 2.71e-06
 $ p.adjust   : num 0.00697
 $ qvalue     : num 0.00637
 $ geneID     : chr "ENSG00000166311/ENSG00000186104/ENSG00000175054/ENSG00000158161/ENSG00000004487/ENSG00000151876/ENSG00000067369"| __truncated__
 $ Count      : int 9
#...Citation
 T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu.
 clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.
 The Innovation. 2021, 2(3):100141 