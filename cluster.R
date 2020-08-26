library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(carData)
library(org.Hs.eg.db)
library(GOplot)
library(stringr)


# ���õ�GO�����Ľ�� ---------------------------------------------------------------------

setwd("D:\\E\\��ʿ\\R_����\\HGSOC\\Data")
coef_Bridge <- read.csv("gene_feature.csv", header=TRUE, sep = ',')
x1 <- as.matrix(coef_Bridge)
x1 <- as.character(x1)
View(x1)

eg1 <- bitr(x1, 
            fromType="SYMBOL", 
            toType=c("ENTREZID","ENSEMBL"), 
            OrgDb="org.Hs.eg.db"); 
head(eg1)

gene <- as.matrix(eg1[,2])

ego <- enrichGO(gene          = gene,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.2,   
                qvalueCutoff  = 0.2,   
                readable      = TRUE)

ego2 = data.frame(ego)
# write.csv(ego, file = "ego.csv",row.names = F)


library(stringr)
GO=ego2[1:12,c(1,2,8,6)]  # ��������ʾ12������
# View(GO)
GO$geneID=str_replace_all(GO$geneID,"/",",")
names(GO)=c("ID","Term","Genes","adj_pval")
GO$Category="BP"


# ����ͼ -----------------------------------------------------------------------

geneList <- read.csv("dif_FC.csv", header=TRUE, sep = ',')

## merge
geneList_1 <- merge(coef_Bridge, geneList, by.x="x", by.y = "X", all=FALSE) # logreg
gene = geneList_1[,c(1,2)]
names(gene)[1]="ID"

plot_data=list(DA=GO,GE=gene)
circ2=data.frame()
circ2=circle_dat(plot_data$DA,plot_data$GE)
## chord_dat �в���ȫ0��
chord <- chord_dat(data = circ2, genes = gene, process = unique(circ2$term))

GOChord(chord)

# Display only genes which are assigned to at least three processes
GOChord(chord, limit = c(3,5))

GOHeat(chord[,-10], nlfc = 0)

# pdf(file = "go.pdf",width = 19,height = 14)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)
# dev.off()

## ���� ����logFC����
GOCluster(circ2, unique(circ2$term)[c(1:4,11:13)], clust.by = 'logFC', term.width = 2)
## ���� ����Term����

# pdf(file = "cluster.pdf",width = 16,height = 12)
GOCluster(circ2, unique(circ2$term)[c(1:12)], clust.by = 'term', lfc.col = c('darkgoldenrod1', 'black', 'cyan1'))
# dev.off()

