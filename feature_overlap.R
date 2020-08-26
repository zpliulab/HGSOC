setwd("D:\\E\\博士\\R_程序\\HGSOC\\Data")
# load(".RData")    # 载入工作区

# 读入数据 --------------------------------------------------------------------

svm <- read.csv("gene_rfe_svm.csv", header=TRUE, sep = ',')
knn <- read.csv("gene_rfe_knn.csv", header=TRUE, sep = ',')
rf <- read.csv("gene_rfe_rf.csv", header=TRUE, sep = ',')
ab <- read.csv("gene_rfe_ab.csv", header=TRUE, sep = ',')
nn <- read.csv("gene_rfe_nnet.csv", header=TRUE, sep = ',')

svm <- as.matrix(svm)
knn <- as.matrix(knn)
rf <- as.matrix(rf)
ab <- as.matrix(ab)
nn <- as.matrix(nn)
  
svm_1 <- intersect(svm, knn)
svm_2 <- intersect(svm, rf)
svm_3 <- intersect(svm, ab)
svm_4 <- intersect(svm, nn)

knn_1 <- intersect(knn, rf)
knn_2 <- intersect(knn, ab)
knn_3 <- intersect(knn, nn)

rf_1 <- intersect(rf, ab)
rf_2 <- intersect(rf, nn)

ab_1 <- intersect(ab, nn)

gene <- union(svm_1, union(svm_2, union(svm_3, union(svm_4, union(knn_1, union(knn_2, union(knn_3, union(rf_1, union(rf_2, ab_1)))))))))
View(gene)

# write.csv(gene, file = "gene_feature.csv", row.names = F)


# install.packages('pheatmap')
library(pheatmap)

Data = read.table("matrix_DE.txt", header = T, check.names = FALSE)
selected <- Data[gene,]


# pdf(file = "差异表达基因热图822.pdf",width = 12,height = 12)

Label = read.table("phe_zh.txt",header = TRUE, sep = "\t")
Label <- factor(Label[,"class"])
Label <- data.frame(Label)

rownames(Label) = colnames(selected[1:52,])
pheatmap(selected[1:52,], annotation_col = Label)


# heatmap -----------------------------------------------------------------

pheatmap(selected[1:52,],annotation_col = Label, 
         color = colorRampPalette(c("blue", "white","red"))(100),
         fontsize_row = 5,scale = "row",cutree_cols = 2,border_color = NA)
         # filename = c("差异表达基因热图822.pdf"))

# dev.off()

