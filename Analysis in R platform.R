# Analysis in R platform
# R v4.4.2
# Gene expression filtering
library(dplyr)
library(tidyverse)
library(DESeq2)

#-----E17 filtering---
setwd("../data")
data_E17 = read.csv("E17.csv", sep=",", header=T)
gene_mean <- rowMeans(data_E17[, -1])
gene_detection_count <- apply(data_E17[, -1] > 0, 1, sum)
data_E17_filtered <- data_E17[gene_mean > 1 & gene_detection_count >= 9, ]
data_E17_filtered_unique <- data_E17_filtered %>%
  distinct(gene_id, .keep_all = TRUE)
dim(data_E17_filtered_unique)
write.csv(data_E17_filtered_unique, "1.data_E17_filtered_unique.csv", row.names = F)

# -----------------------------E17Identifying significant genes between the groups---------------------------------------
# 
# 	Compared groups:
#		(1) CC vs. CR
#		(2) CC vs. RR
#		(3) CR vs. RR

# ---------Set functions------------

select_countData <- function(countData, group_name, comparison){
  library(dplyr)
  Gene1 <- colnames(countData)[1]
  m <- as.data.frame(t(as.matrix(countData[,-1])))
  m$group <- group_name
  m <- filter(m, group==comparison[1] | group==comparison[2])
  m$group <- NULL
  m <- as.data.frame(t(m))
  m <- cbind(countData[1], m)
  colnames(m)[1] <- Gene1
  rownames(m) <- c(1:nrow(m))
  return(m)
}
# select metadata from the groups to be compared从要比较的组中选择元数据
select_metaData <- function(metaData, group_name, comparison){
  m <- metaData
  m$group <- group_name
  m <- filter(metaData, group==comparison[1] | group==comparison[2])
  return(m)
}

# construct DESeqDataSet Object 
dds_object <- function(countData, metaData){
  dds <- DESeqDataSetFromMatrix(countData=countData, 
                                colData=metaData, 
                                design=~group, tidy = TRUE)
  dds$group <- factor(dds$group, levels = comparison)
  keep <- rowSums(counts(dds)) >= 10	# Pre-filtering
  dds <- dds[keep,]
  dds <- DESeq(dds)	
  return(dds)
}

# select significant genes between groups
select_siggene <- function(dds){
  res <- results(dds) 
  alpha <- 0.99
  siggene <- res[which(res$padj < alpha), ]
  return(siggene)
}

m <- data_E17_filtered_unique

countData_original <- m
metaData_original <- read.csv("1.group.csv", header=TRUE, sep=',', na.strings=9999, stringsAsFactor=TRUE)
metaData_original[ ,"group"] <- factor(metaData_original[ ,"group"], level=c("CC", "CR", "RR"))


# (1) CC vs. RR
comparison <- c("CC", "RR")
countData <- select_countData(countData_original, group_name=metaData_original$group, comparison=comparison)
metaData <- select_metaData(metaData_original, group_name=metaData_original$group, comparison)
dds <- dds_object(countData, metaData)  
res <- results(dds) 
res[order(res$pvalue),]
siggene <- select_siggene(dds)  


res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)

res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

res2<- res1[which(res1$padj < 1),]     
res_total <- rbind(res2)
write.csv(res_total,file="2.CC-RR_signif_gene_total.csv",quote = F)

# (2) CR vs. CC
comparison <- c("CR", "CC")
countData <- select_countData(countData_original, group_name=metaData_original$group, comparison=comparison)
metaData <- select_metaData(metaData_original, group_name=metaData_original$group, comparison)
dds <- dds_object(countData, metaData)  
res <- results(dds) 
res[order(res$pvalue),]
siggene <- select_siggene(dds)  


res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)

res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

res2<- res1[which(res1$padj <= 1.1),]     
res_total <- rbind(res2)
write.csv(res_total,file="2.CR-CC_signif_gene_total.csv",quote = F)

# (3) CR vs. RR
comparison <- c("CR", "RR")
countData <- select_countData(countData_original, group_name=metaData_original$group, comparison=comparison)
metaData <- select_metaData(metaData_original, group_name=metaData_original$group, comparison)
dds <- dds_object(countData, metaData)  
res <- results(dds) 
res[order(res$pvalue),]
siggene <- select_siggene(dds)  


res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)

res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

res2<- res1[which(res1$padj < 1),]     
res_total <- rbind(res2)
write.csv(res_total,file="2.CR-RR_signif_gene_total.csv",quote = F)








