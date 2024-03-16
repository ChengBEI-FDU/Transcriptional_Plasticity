library(openxlsx)
library(tidyverse)
library(vegan)
library(moments)
library(ggplot2)
library(rtracklayer)
library(reshape2)
library(Rmisc)
library(stringr)
library(ggsci)
library(Rtsne)
library(gtools)
library(ggpubr)
library(pheatmap)
library(edgeR)
library(limma)
library(data.table)


setwd("/Users/beicheng/Desktop/MTB Transcriptome the Big Data/results9/article/Revision_Version_0127/Data and code/")

#### Function ####
DEG_limma <- function(Treat_sample, 
                      Control_sample,
                      count_data,
                      log_FC = 0, 
                      adj_p_val = 0.05) {
  tmp_count <- count_data[,c(Treat_sample, Control_sample)]
  tmp_group <- factor(c(rep('Treat',length(Treat_sample)), rep('Control',length(Control_sample))))
  tmp_DGEList <- DGEList(counts = tmp_count, group = tmp_group)
  tmp_DGEList$counts <- tmp_DGEList$counts[rowSums(cpm(tmp_DGEList$counts)>2)>=2,]
  tmp_DGEList <- calcNormFactors(tmp_DGEList)
  tmp_design <- model.matrix(~0 + tmp_group)
  rownames(tmp_design) <- colnames(tmp_DGEList$counts)
  tmp_contrasts <- makeContrasts('Treat_vs_Control' = tmp_groupTreat-tmp_groupControl, levels = tmp_design) 
  tmp_voom <- voomWithQualityWeights(counts = tmp_DGEList, design = tmp_design, plot = F)
  tmp_fit <- lmFit(tmp_voom, design = tmp_design)
  tmp_fit <- contrasts.fit(tmp_fit, contrasts = tmp_contrasts)
  tmp_fit <- eBayes(tmp_fit)
  tmp_top <- topTable(tmp_fit, coef = 'Treat_vs_Control', adjust.method = 'fdr', number = 'Inf', sort.by = "P")
  tmp_deg <- tmp_top[(abs(tmp_top$logFC) > log_FC & tmp_top$adj.P.Val < adj_p_val),]
  tmp_deg$gene <- rownames(tmp_deg)
  return(tmp_deg)
}


#### Input data ####
phenotype_dt <- read.xlsx("data/benchmark_DE/data_CellSystems2021_BreeAldridge/mmc3-2.xlsx", sheet = 1)
exp_dt <- read.csv("data/RPKM Mtb.csv", row.names = 1)
exp_sample <- read.xlsx("data/benchmark_DE/DE_data/Sample info 894_Mtb_for_DEG.xlsx", sheet = 1) # latest version that have 129 comparisons
gene_dt <- read.csv("data/TP and feature Mtb.csv", row.names = 1)
rawcount_dt <- read.csv("data/tmp data/Raw Counts Mtb.csv", row.names = 1)


gtf_dt_Mtb <- as.data.frame(import("data/original data/genome ref/Mycobacterium_tuberculosis_H37Rv_gff_v4.gff"))


#### DEG ####
deg_sample <- exp_sample[!is.na(exp_sample$DE_index),]
DE_index <- unique(deg_sample$DE_index) 

for (i in DE_index) { 
  cat(i)
  cat('\r')
  Treat_sample <- deg_sample$Run[(deg_sample$DE_index %in% i) & (deg_sample$DE_group %in% "t")]
  Control_sample <- deg_sample$Run[(deg_sample$DE_index %in% i) & (deg_sample$DE_group %in% "c")]
  DEG <- DEG_limma(Treat_sample = Treat_sample, Control_sample = Control_sample,
                   count_data = rawcount_dt, log_FC = 0, adj_p_val = 0.05)
  if (nrow(DEG) > 0) { # remove results without DEGs
    DEG$DE_index <- i
    write.csv(DEG, paste0("data/benchmark_DE/DE_results/DE_index_",i,".csv"))
    }
  rm(Treat_sample, Control_sample, DEG)
}
rm(i)



#### Bind data ####
fileName <- list.files(path = "data/benchmark_DE/DE_results/")
DEG_bind <- data.frame()
for(i in fileName){
  cat(i)
  cat('\r')
  DEG_bind <- rbind(DEG_bind, read.csv(paste0("data/benchmark_DE/DE_results/", i, collapse = ""), header = TRUE))
}
rm(i)

paste0("Index ", setdiff(DE_index, unique(DEG_bind$DE_index)), 
       " is not included because no significant DEG was identified", collapse = "; ")


#### Dcast ####
DEG_dcast <- reshape2::dcast(DEG_bind, gene~DE_index, value.var = 'logFC')
DEG_dcast[is.na(DEG_dcast)] <- 0
rownames(DEG_dcast) <- DEG_dcast$gene
DEG_dcast <- DEG_dcast[,-1]

write.csv(DEG_dcast, "data/benchmark_DE/DE_results/DEG dcast.csv")

# heatmap
bk <- c(seq(-9,-0.01,by=0.01),seq(0,9,by=0.01))
pdf('data/benchmark_DE/DE_results/Heatmap of log2fc.pdf', width = 15, height = 12)
pheatmap(DEG_dcast, border_color = NA, show_rownames = F, angle_col = 45, 
         color = c(colorRampPalette(colors = c("navy", "navy", "white"))(length(bk)/2),
                   colorRampPalette(colors = c("white", "firebrick3", "firebrick3"))(length(bk)/2)),
         breaks = bk)
dev.off()
rm(bk)



#### q5-q95 ####
rm(list = ls(pattern = ""))

# input data
exp_dt <- read.csv("data/RPKM Mtb.csv", row.names = 1)
feature <- read.csv("data/TP and feature Mtb.csv", row.names = 1)
log2exp_dt <- log2(exp_dt+1)

#### q5-q95 ####
# function
calculating_ExprAnnotation <- function(expr, # RPKM
                                       log = T
) {
  # whether log expression data
  if (log) {expr <- log2(expr+1)}
  
  # median and IQR(log2(RPKM+1))
  expr.q5.logRPKM <- apply(expr, 1, function(x) {quantile(x,0.05)})
  expr.q95.logRPKM <- apply(expr, 1, function(x) {quantile(x,0.95)})
  
  # MERGE
  benchmark <- data.frame(gene = rownames(expr),q5 = expr.q5.logRPKM, q95 = expr.q95.logRPKM)
  return(benchmark)
}

# calculate
q95 <- calculating_ExprAnnotation(exp_dt)

benchmark_mtb <- merge(q95, feature, by = 'gene', all.y = T)
benchmark_mtb <- benchmark_mtb[,c(1,2,3,7)]

benchmark_mtb$q5_log2FC <- benchmark_mtb$q5 - benchmark_mtb$Mean_logRPKM
benchmark_mtb$q95_log2FC <- benchmark_mtb$q95 - benchmark_mtb$Mean_logRPKM

# output
write.csv(benchmark_mtb,"data/benchmark_log2FC at q5 and q95_Mtb.csv")




