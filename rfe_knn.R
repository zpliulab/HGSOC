# install.packages("lattice")
library(caret)
library(lattice)
library(ggplot2)

setwd("D:\\E\\博士\\R_程序\\HGSOC\\Data")

eset <- read.table("matrix_DE.txt",header = TRUE,sep = "\t")

Control <- rfeControl(functions = caretFuncs, method = "cv",
                      verbose = FALSE , returnResamp = "final")

trControl1 <- trainControl( method = "cv",
                            classProbs=TRUE,
                            summaryFunction = twoClassSummary)


disease <- as.factor(c(rep("HGSOC",10), rep("Normal",10)))

# KNN—RFE -----------------------------------------------------------------

rf2 <- rfe(t(eset), disease, c(1:200,210,230,250,300,350,400,450,500),# sizes = c(1:517), 
           rfeControl = Control, trControl = trControl1,
           metric = "Accuracy",  
           method = "knn")
# 两种 结果 一样："ABR"   "ACKR4" "ACSL5"

feature_sele <-rf2$optVariables

# write.table(feature_sele, file = "ranklist_KNNrfe.txt", quote=F, sep="\t")

