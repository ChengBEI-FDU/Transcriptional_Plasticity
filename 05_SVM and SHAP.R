#### Import packages ####
library(openxlsx)
library(ggplot2)
library(ggsci)
library(reshape2)
library(ggpubr)
library(ggstream)
library(rtracklayer)
library(e1071)
library(caret)
library(ggpmisc)
library(stringr)
library(dplyr)
library(Rmisc)


#### Input data ####
sample_info <- read.xlsx("./data/sample/Sample Mtb.xlsx", sheet = 1)
TP_data <- read.csv("./data/TP and feature Mtb.csv", row.names = 1)
lgbm_feature_importance <- read.csv("./data/LGBM feature importance.csv", row.names = 1)
feature_annotation <- read.xlsx("./data/Feature Annotation.xlsx", sheet = 1)


#### SVM ####
library(e1071)

# 0-1 scaled
feature_top_Mtb <- TP_data[,c(16,18,44,73,81)]

feature_top_Mtb$width <- (feature_top_Mtb$width - min(feature_top_Mtb$width, na.rm = T)) / 
  (max(feature_top_Mtb$width, na.rm = T) - min(feature_top_Mtb$width, na.rm = T))
feature_top_Mtb$basePercent_GC <- (feature_top_Mtb$basePercent_GC - min(feature_top_Mtb$basePercent_GC, na.rm = T)) / 
  (max(feature_top_Mtb$basePercent_GC, na.rm = T) - min(feature_top_Mtb$basePercent_GC, na.rm = T))
feature_top_Mtb$operon_length <- (feature_top_Mtb$operon_length - min(feature_top_Mtb$operon_length, na.rm = T)) / 
  (max(feature_top_Mtb$operon_length, na.rm = T) - min(feature_top_Mtb$operon_length, na.rm = T))
feature_top_Mtb$Activator.number <- (feature_top_Mtb$Activator.number - min(feature_top_Mtb$Activator.number, na.rm = T)) / 
  (max(feature_top_Mtb$Activator.number, na.rm = T) - min(feature_top_Mtb$Activator.number, na.rm = T))

feature_top_rmNA_Mtb <- na.omit(feature_top_Mtb)

# SVM
fit_svm_top <- svm(TP ~ ., data = feature_top_rmNA_Mtb)
summary(fit_svm_top)
pred_svm_top <- predict(fit_svm_top, feature_top_rmNA_Mtb[,2:5])

min(pred_svm_top)
max(pred_svm_top)


#### SHAP ####
library(DALEX)
library(iBreakDown)

# DALEX explain
svm_explain <- DALEX::explain(fit_svm_top, data = feature_top_rmNA_Mtb[,2:5], y = feature_top_rmNA_Mtb$TP, type = 'regression')

# shap for all genes
svm_shap_all <- data.frame()
for (i in 1:nrow(feature_top_rmNA_Mtb)) {
  cat(paste0(i,"/",nrow(feature_top_rmNA_Mtb),collapse = ""))
  cat('\r')
  tmp <- shap(svm_explain, feature_top_rmNA_Mtb[i,], keep_distributions = T)
  tmp <- as.data.frame(tmp)
  tmp$gene_index <- i
  svm_shap_all <- rbind(svm_shap_all, tmp)
}

write.csv(svm_shap_all, "./data/SVM_shap_feature_importance_TSS added.csv")

