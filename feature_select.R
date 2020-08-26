library(dplyr)        
library(tidyr)
library(tidyverse)   



# load(".RData")
setwd("D:\\E\\博士\\R_程序\\HGSOC\\GSE120196")
Data1 = read.table("GSE120196_scale.txt", header = T, check.names = FALSE)

# setwd("D:\\E\\博士\\R_程序\\HGSOC\\HGSOC的参考数据\\GSE40595")
# Data1 = read.table("GSE40595_scale_38.txt", header = T, check.names = FALSE)

setwd("D:\\E\\博士\\R_程序\\HGSOC\\Data")
gene = read.csv("gene_feature.csv", header=TRUE, sep = ',')
dim(Data1)

colnames(gene) <- c('gene')
Data2 <- cbind(rownames(Data1), Data1)
colnames(Data2) <- c('gene', colnames(Data1))

genedata <- merge(gene, Data2, by = "gene")#[,-2]
genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
genedata2 <- rbind(Data1[1,],genedata1)
rownames(genedata2) <- c('Lable', rownames(genedata1))
# write.table(genedata2,"GSE40595_38_52.txt",quote=F,sep="\t") 


# 在原数据集提取特征 ---------------------------------------------------------------

data1 = read.table("GSE69428_scale.txt", header = T, check.names = FALSE)
dim(data1)
gene <- as.matrix(genedata[,1])

colnames(gene) <- c('gene')
data2 <- cbind(rownames(data1), data1)
colnames(data2) <- c('gene', colnames(data1))

genedata <- merge(gene, data2, by = "gene")
genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
genedata2 <- rbind(data1[1,],genedata1)
rownames(genedata2) <- c('Lable', rownames(genedata1))

# write.table(genedata2,"GSE69428_40595_38_52.txt",quote=F,sep="\t") 