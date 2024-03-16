#### Import packages ####
library(openxlsx)
library(tidyverse)
library(edgeR)
library(limma)
library(vegan)
library(moments)
library(ggplot2)
library(rtracklayer)
library(mclust)
library(ggpubr)
library(sva)
library(ggstream)
library(reshape2)
library(ggsci)


#### Functions ####
ReshapeSampleForTmmRpkm <- function(sampleFile,
                                    LibraryStrategy = "RNA-Seq", 
                                    LibrarySelection = 'cDNA') {
  sample <- sampleFile[sampleFile$LibraryStrategy %in% LibraryStrategy & sampleFile$LibrarySelection %in% LibrarySelection,]
  sample <- data.frame(Run = sample$Run, BioProject = sample$BioProject)
  return(sample)
}

ReshapeGeneLengthForTmmRpkm <- function(gtfFile,
                                        Type) {
  geneLength <- gtfFile[gtfFile$type %in% Type,]
  geneLength <- data.frame(gene = geneLength$Locus, length = geneLength$width)
  geneLength <- geneLength[!duplicated(geneLength$gene),]
  return(geneLength)
} # for Msm and Mab

calculating_TMM_RPKM_allSamples <- function(counts_path, # path of htseq count files (.count)
                                            sample_study, # data frame. Column "Run" is SRA accession (refer to htseq-count file name) and column "BioProject" is bio-projects
                                            gene_length, # data frame. Column "gene" is locus tag and column "length" is gene length
                                            filter_lib_size = 1000000, # a numeric that sample with library size <= 'filter_lib_size' will be excluded
                                            filter_geneLength = 150, # genes shorter than 150 bp will be removed
                                            filter_zeroExprGenes = TRUE, # whether remove genes with zero counts in all samples
                                            filter_gene = NA # a vector contains genes that are excluded
) {
  # set path
  setwd(counts_path)
  
  # read count files
  counts_data <- readDGE(list.files(pattern =".count"), columns=c(1,2))
  counts_data$counts <- counts_data$counts[1:(nrow(counts_data$counts)-5),]
  raw.counts <- as.data.frame(counts_data$counts)
  
  # filter library size
  lib_size <- counts_data$samples
  count_sample <- rownames(lib_size)[lib_size$lib.size > filter_lib_size] 
  count_sample <- count_sample[count_sample %in% sample_study$Run]
  raw.counts <- raw.counts[,colnames(raw.counts) %in% count_sample]
  
  # filter genes
  if (filter_zeroExprGenes) {raw.counts <- raw.counts[apply(raw.counts, 1, function(x) {sum(x) > 0}),]}
  gene_geneLength <- gene_length$gene[gene_length$length >= filter_geneLength]
  count_gene <- rownames(raw.counts)
  count_gene <- count_gene[count_gene %in% gene_geneLength & (!count_gene %in% filter_gene)]
  raw.counts <- raw.counts[rownames(raw.counts) %in% count_gene,]
  
  # re-shape data
  count_gene <- rownames(raw.counts)
  
  gene_length <- gene_length[gene_length$gene %in% count_gene,]
  gene_length <- as.data.frame(gene_length[match(gene_length$gene, count_gene),])
  gene_length <- data.frame(row.names = gene_length$gene, length = gene_length$length)
  
  sample_study <- sample_study[sample_study$Run %in% count_sample,]
  sample_project <- sample_study$BioProject[!duplicated(sample_study$BioProject)]
  
  # calculate TMM normalized RPKM 
  temp.count <- raw.counts
  
  temp.counts.data <- DGEList(counts = temp.count)
  temp.counts.data <- calcNormFactors(temp.counts.data)
  temp.counts.data$genes <- gene_length
  
  RPKM <- as.data.frame(rpkm(temp.counts.data))
  return(RPKM)
}


#### Input data ####
# sample
sample_Mtb <- read.xlsx("./data/sample/Mtb_SRA.xlsx")
sample_Mab <- read.xlsx("./data/sample/Mabsc_SRA.xlsx")
sample_Msmeg <- read.xlsx("./data/sample/Msmeg_SRA.xlsx")

sample_Mtb <- sample_Mtb[sample_Mtb$`Inclusion.(RvGeneExprVar)` %in% 'y',]
sample_Msmeg <- sample_Msmeg[sample_Msmeg$Strain %in% "MC2155",]
sample_Mab <- sample_Mab[!is.na(sample_Mab$Treat),]

sample_Mtb <- sample_Mtb[,c(1,3)]
sample_Mab <- ReshapeSampleForTmmRpkm(sampleFile = sample_Mab)
sample_Msmeg <- ReshapeSampleForTmmRpkm(sampleFile = sample_Msmeg)

# gene length
gtf_Mtb <- as.data.frame(import("./data/genome_ref/Mtb_GCF_000195955.2_ASM19595v2_genomic.gff"))
gtf_Mab <- as.data.frame(import("./data/genome_ref/Mycobacterium_abscessus_ATCC_19977_gff_v4.gff"))
gtf_Msmeg <- as.data.frame(import("./data/genome_ref/Mycobacterium_smegmatis_MC2-155_gff_v4.gff"))

geneLength_Mtb <- data.frame(gene = gtf_Mtb$locus_tag[gtf_Mtb$type %in% c("gene","pseudogene")], length = gtf_Mtb$width[gtf_Mtb$type %in% c("gene","pseudogene")])
geneLength_Mab <- ReshapeGeneLengthForTmmRpkm(gtfFile = gtf_Mab, Type = 'CDS')
geneLength_Msmeg <- ReshapeGeneLengthForTmmRpkm(gtfFile = gtf_Msmeg, Type = 'CDS')


#### Calculate TMM normalized RPKM ####
# Mtb
filter_gene_mtb <- gtf_Mtb$locus_tag[gtf_Mtb$gene %in% c("rrl","rrs","rrf","rnpB")] # remove ribosome RNA
filter_gene_mtb <- unique(c(filter_gene_mtb,geneLength_Mtb$gene[grep(pattern="nc", geneLength_Mtb$gene)])) # remove noncoding RNA

rpkm_Mtb <- calculating_TMM_RPKM_allSamples(counts_path = "./data/htseqcount/Mtb/", 
                                            sample_study = sample_Mtb, gene_length = geneLength_Mtb,
                                            filter_gene = filter_gene_mtb)
setwd("../../../")
write.csv(rpkm_Mtb, "./data/unfiltered tmm rpkm Mtb.csv")


# Msm
rpkm_Msmeg <- calculating_TMM_RPKM_allSamples(counts_path = "./data/htseqcount/Msmeg/", 
                                              sample_study = sample_Msmeg, gene_length = geneLength_Msmeg)
setwd("../../../")
write.csv(rpkm_Msmeg, "./data/unfiltered tmm rpkm Msm.csv")


# Mab
pla.gene_Mab <- gtf_Mab$Locus[gtf_Mab$seqnames %in% "NC_010394.1"]
rpkm_Mab <- calculating_TMM_RPKM_allSamples(counts_path = "./data/htseqcount/Mab/", 
                                            sample_study = sample_Mab, gene_length = geneLength_Mab,
                                            filter_gene = pla.gene_Mab) # remove plasmid genes
setwd("../../../")
write.csv(rpkm_Mab, "./data/unfiltered tmm rpkm Mab.csv")


