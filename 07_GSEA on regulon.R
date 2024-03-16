library(clusterProfiler)

# input data
feature_Mtb <- read.csv("./data/TP Mtb.csv", row.names = 1)
regulon_gmt <- read.gmt("./data/selected_regulon.gmt")

# run gsea
set.seed(2023)

gsea_geneList <- feature_Mtb$TP
names(gsea_geneList) <- feature_Mtb$gene
gsea_geneList <- sort(gsea_geneList, decreasing = T)

gsea_results <- GSEA(geneList = gsea_geneList, TERM2GENE = regulon_gmt, verbose = F, eps=1e-10, minGSSize = 1, pvalueCutoff = 1)

# results
gsea_results@result$ID
gsea_results@result$p.adjust
gsea_results@result$NES

# output
gsea_results_dt <- as.data.frame(gsea_results@result)
colnames(gsea_results_dt)[1] <- 'name'
write.csv(gsea_results_dt, "./data/TP of regulons_GSEA.csv")

