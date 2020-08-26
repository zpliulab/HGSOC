
# setwd("C:\\Users\\2\\Desktop\\R_study")
# install.packages("BiocGenerics")
# install.packages("parallel")
# install.packages("Biobase")
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("fdrtool")
# install.packages("tidyverse")
# install.packages("stringr")
# install.packages("data.table")

# install.packages("devtools")
# devtools::install_github("tidyverse/dplyr")

library(stringr)
library(BiocGenerics)
library(parallel)
library(Biobase)
library(dplyr)       # ��>�� �ܵ������ĵ��ã�����
library(tidyr)
library(tidyverse)   # tibble �ĵ���
library(fdrtool)     # fdrУ��
library(data.table)  # �� fread ����.soft �ļ�
library(fdrtool)     # ���� p ֵ
library(data.table)  # ʹ�� fread ��������
setwd("D:\\E\\��ʿ\\R_����\\HGSOC\\Data")

## ��������������
exprSet <- read.table("GSE69428_series_matrix.txt",header=T,sep='\t',fill=TRUE,strip.white = T)
exprSet$ID_REF <- as.character(exprSet$ID_REF)

## ����оƬ̽��ע���ļ�
anno <- read.table("GPL570.txt",header=T,sep='\t',fill=TRUE,strip.white = T,quote = "") 

library(dplyr)
library(tidyr)
anno2 <- anno %>%    
  select(ID,Gene.ID) %>%             
  filter(Gene.ID != '') 

colnames(anno2) <- c('ID_REF','EntrezID')    
anno2$ID_REF <- as.character(anno2$ID_REF)  


## ���������������оƬע���ļ���̽�������ж�Ӧ
exprset2 <- exprSet %>%                      
  inner_join(anno2,by='ID_REF') %>%        
  select(ID_REF,EntrezID, everything())   
                          
## ����оƬע���ļ���������һ��̽���Ӧ�������Ĳ�ֿ�
exprset3 <- exprset2
a <- tibble(exprset3[,1:2])

test1 <- apply(a,1, function(x){
  str_split(x[2],'///',simplify=T)      
} )

test2 <- apply(a, 1, function(x){           
  paste(x[1],str_split(x[2],'///', simplify=T), sep = "---")
})

unlist(test2)                              

x <- tibble(unlist(test2))                 
colnames(x) <- "lala"                       


x2 <- separate(x,lala,c("id","entrezID"),sep = '---')    
x3 <- merge(x2,exprset3,by.x = "id",by.y="ID_REF",all=FALSE)  
x4<-x3[,-c(1,3)]                      

zz <- as.matrix(apply(as.matrix(x4[,1]),1,function(x) as.numeric(x)))

XX <- x4[,-1]
colnames(XX)[1:3]
XX1 <- cbind(zz,XX)
colnames(XX1) <- c("entrezID",colnames(XX))


## �û���id�������õ�оƬע���ļ����л������ĸ���
homo<-read.table("homo.txt",header=T,sep='\t')

x5 <- merge(homo,XX1,by.x="GeneID",by.y = "entrezID",all=FALSE) 
# dim(x5)    # 31773    31


## ̽����ƥ���������ȡ�����̽���Ӧһ����������ݼ���IQR������IQR����̽������
expset4 <- x5 %>%
  dplyr::select(-GeneID) %>%              # ȥ��������Ϣ
  mutate(rowIQR =apply(.[,-1],1,IQR)) %>% # ����ÿ�е�IQR
  arrange(desc(rowIQR)) %>%               # �ѱ�������ƽ��ֵ���Ӵ�С����
  distinct(Symbol,.keep_all = T) %>%      # symbol���µ�һ��
  dplyr::select(-rowIQR) %>%                 # ����ѡ��ȥ��rowIQR��һ��
  tibble::column_to_rownames(colnames(.)[1]) # �ѵ�һ�б��������ɾ��
View(expset4[1:10,1:10])
dim(expset4)   # 17689    29
# write.table(expset4,"GSE69428_expr.txt",quote=F,sep="\t")  


lable = read.csv("GSE69428_all.csv", header = T, sep=',')
lable_1 <- lable[which(lable$Sample.type == "tissue"),]
expset5 <- expset4[,which(lable$Sample.type == "tissue")]
dim(expset5)          # 17689    20  


lable2 = read.csv("GSE69428_20.csv", header = T, sep=',')
dim(lable2)    # 20   2
data = rbind(as.matrix(t(lable2[,2])), as.matrix(expset5))
rownames(data) <- c('Lable', rownames(expset5))
View(data[1:10,])
dim(data)    # 17690    20
# write.table(data,"GSE69428_outcome.txt",quote=F,sep="\t")   



# �����仯��scaleһ�� ------------------------------------------------------------

setwd("D:\\E\\��ʿ\\R_����\\HGSOC\\Data")

## ��������
data = read.table("GSE69428_outcome.txt", header = T, sep='\t', fill=TRUE, strip.white = T, check.names = F)

# my_scale -------------------------------------------------------------------

my_scale <- function(x){
  x1 <- cbind(t(x[1,]), scale(t(x[-1,])))
  x2 <- t(x1)
  return(x2)
}

# ���� ----------------------------------------------------------------------

data1 <- my_scale(data)
dim(data1)    # 17690    20
# write.table(data1,"GSE69428_scale.txt",quote=F,sep="\t")  
