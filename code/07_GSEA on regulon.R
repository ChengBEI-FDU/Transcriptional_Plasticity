library(clusterProfiler)

setwd("/Users/beicheng/Desktop/MTB Transcriptome the Big Data/results9/article/Revision_Version_0127/Data and code/")

# input data
feature_Mtb <- read.csv("data/TP Mtb.csv", row.names = 1)
regulon_gmt <- read.gmt("data/original data/mSphere 2022 regulon/selected_regulon.gmt")

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
write.csv(gsea_results_dt, "./data/TP StatTest in Regulons using GSEA method.csv")
