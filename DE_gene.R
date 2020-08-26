
library(simpleaffy)
library(affyPLM)
library(RColorBrewer)
library(graph)
library(genefilter)
library(affycoretools)
library(affy)
library(limma)

setwd("D:\\E\\��ʿ\\R_����\\HGSOC\\Data")

Data = read.table("GSE69428_scale.txt", header = T, check.names = FALSE)

# ���ͱ�ǩ ----------------------------------------------------------------------

sample <- as.matrix(colnames(Data))
label <- t(Data[1,])
class <- cbind(sample, label)
colnames(class) <- c("sample","class")
# write.csv(class, file = "phenotype.csv", row.names = F)


# ���� �� ��������----------------------------------------------------------------------

Data1 <- Data[,which(Data[1,] == 1)]

Data0 <- Data[,which(Data[1,] == 0)]

data <- cbind(Data0,Data1)[-1,]    # ע�⣡��-- Ҫ�����ݶ�Ӧ�ã��Ķ���0/1


phe = read.csv("phenotype.csv", header=TRUE, sep = ',')
phe1 <- phe[which(phe[,2] == 1),]
phe1[,2] <- c("Normal")
phe0 <- phe[which(phe[,2] == 0),]
phe0[,2] <- c("HGSOC")
phenotype <- rbind(phe0,phe1)
# write.table(phenotype, file = "phe_zh.txt", quote=F, sep="\t", row.names = F)

# ������� ----------------------------------------------------------------
disease = read.table("phe_zh.txt",header = TRUE, sep = "\t")
disease <- factor(disease[,"class"])

design <- model.matrix(~-1+disease)
contrast.matrix <- makeContrasts(contrasts = "diseaseNormal - diseaseHGSOC", levels = design)
fit <- lmFit(data,design)
fit1 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit1)
dif <- topTable(fit2,coef = "diseaseNormal - diseaseHGSOC",n = nrow(fit2),lfc = log2(2))
# dif <- topTable(fit2,coef = "diseaseNormal - diseaseHGSOC",n = nrow(fit2))
# dif <- dif[dif[,"P.Value"]<0.05,]
dif0.05 <- dif[dif[,"adj.P.Val"] < 0.05,]    #  2089    6
dif_FC <- dif0.05[abs(dif0.05[,"logFC"]) > 1.5,]     # 517   6
# dim(dif_FC)
# write.csv(dif, file = "dif.csv")  # ѡ�񲻴�lfc = log2(2)����������
# write.csv(dif_FC, file = "dif_FC.csv")

gene_adjp <- rownames(dif_FC)
# write.csv(gene_adjp, file = "gene_adjp.csv", row.names = F)


# ����gene�ı������� -------------------------------------------------------------

# merge
data <- read.table("GSE69428_scale.txt", header = T, check.names = FALSE)
data1 <- cbind(rownames(data),data)
colnames(data1) <- c("gene_symbol", colnames(data))
data2 <- data1[-1,]  

gene_adjp = read.csv("gene_adjp.csv", header=TRUE, sep = ',')
colnames(gene_adjp) <- c("gene_adjp")

data3 <- merge(gene_adjp, data2, by.x="gene_adjp",by.y = "gene_symbol",all=FALSE) 
data4 <- rbind(data1[1,-1],data3[,-1])
row.names(data4) <- c("Label",as.character(data3[,1]))
dim(data4)    # 518  20
# write.table(data4, file = "GSE69428_scale_adjp.txt",quote = F, sep = "\t")

## ɾ��label, �õ�rfe����Ҫ�ľ���
data5 <- data4[-1,]
# write.table(data5, file = "matrix_DE.txt",quote = F, sep = "\t")


# �����������о��� -----------------------------------------------------------------------

eset0.001<- Data[row.names(dif0.05),]
pearson_cor_0.001 <- cor(eset0.001)
dist.lower.0.001 <- as.dist(1 - pearson_cor_0.001)
hc0.001 <- hclust(dist.lower.0.001,"ave")
# pdf(file = "�������������ͼsample.pdf",width = 8,height = 6)
plot(hc0.001)
# dev.off()
# write.table(data,file = "matrix.txt",quote = FALSE,sep = "\t")


# ���ջ��������о���
eset0.001<- Data[row.names(dif0.05),]
rownames(eset0.001) <- c(rownames(dif0.05))
eset0.001_t <- t(eset0.001)
pearson_cor_0.001_t <- cor(eset0.001_t)
dist.lower.0.001_t <- as.dist(1 - pearson_cor_0.001_t)
hc0.001_t <- hclust(dist.lower.0.001_t,"ave")
# pdf(file = "�������������ͼgene.pdf",width = 60,height = 25)
plot(hc0.001_t)
# dev.off()

# install.packages('pheatmap')
library(pheatmap)
selected <- Data[rownames(dif0.05),]
rownames(selected) <- c(rownames(dif0.05))
# pdf(file = "������������ͼ.pdf",width = 14,height = 7)  # c("navy", "white", "firebrick3") c("green","black","red")
pheatmap(selected[1:52,],color = colorRampPalette(c("green","black","red"))(100), # ������ǰ20���������������ͼ
         fontsize_row = 4,scale = "row",cutree_cols = 2,border_color = NA)
# dev.off()
