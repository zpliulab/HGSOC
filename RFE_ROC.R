
# load(".RData")    # 载入工作区

library(e1071)
library(limma)
library(ggplot2)
library(reshape2)
library(gmodels)


# 读入数据 --------------------------------------------------------------------
setwd("D:\\E\\博士\\R_程序\\HGSOC\\Data")
eset <- read.table("matrix_DE.txt",header = TRUE,sep = "\t") 
eset <- as.matrix(eset)


# 读入排序文件 ------------------------------------------------------------------

rfe_svm <- read.table("ranklist_KNNrfe.txt",stringsAsFactors = FALSE)
colnames(rfe_svm) <- c("rank","chipID")
rfe_svm_dif<- rfe_svm[468:517,]
eset_svm <- eset[rfe_svm_dif$chipID,]

# 差异表达gene ----------------------------------------------------------------

library(limma)

data1 <- eset_svm
disease = read.table("phe_zh.txt",header = TRUE, sep = "\t")
disease <- factor(disease[,"class"])

design <- model.matrix(~-1+disease)
contrast.matrix <- makeContrasts(contrasts = "diseaseNormal - diseaseHGSOC", levels = design)
fit_feature <- lmFit(data1,design)
fit1_feature <- contrasts.fit(fit_feature,contrast.matrix)
fit2_feature <- eBayes(fit1_feature)
dif <- topTable(fit2_feature,coef = "diseaseNormal - diseaseHGSOC",n = nrow(fit2_feature),lfc = log2(2))
dif0.05 <- dif[dif[,"adj.P.Val"] < 0.05,]    #  2089    6
dif_feature <- dif0.05[abs(dif0.05[,"logFC"]) > 1.5,]
# write.csv(dif_feature, file = "dif_feature_svm.csv")

# 训练---测试 -----------------------------------------------------------------

x <- t(eset[rfe_svm_dif$chipID,])
gene_ref <- colnames(x)
# write.csv(gene_ref, file = "gene_ref_svm.csv", row.names = F)

y1 <- as.matrix(c(rep(0,10), rep(1,10)))

disease <- as.factor(c(rep("HGSOC",10), rep("Normal",10)))
y <- disease


# Acc ---------------------------------------------------------------------

svm_acc <- list()
cost <- c(rep(10,4),rep(100,4))
gamma <- rep(c(0.1,0.01,0.001,0.0001),2)
parameters <- data.frame(cost,gamma)
for(k in 1:dim(parameters)[1]){
  svm_acc_tmp <- c()
  costtmp <- parameters[k,"cost"]
  gammatmp <- parameters[k,"gamma"]
  all_test_label <- c()
  all_pred_label <- c()
  
  set.seed(666) # for reproducing results
  rowIndices <- 1 : nrow(x) # prepare row indices
  sampleSize <- 0.70 * length(rowIndices) # training sample size
  trainingRows <- sample (rowIndices, sampleSize) # random sampling
  trainingData <- x[trainingRows, ] # training data
  testData <- x[-trainingRows, ] # test data
  trainingLabel <- y[trainingRows]
  testLabel <- y[-trainingRows]
  
  # trainingLabel <- y1[trainingRows]
  # testLabel <- y1[-trainingRows]

  svmfit <- svm (trainingData,trainingLabel, kernel = "radial", cost = costtmp, gamma=gammatmp, scale = FALSE) # radial svm, scaling turned OFF
  eset_pred<- predict(svmfit, testData)
  all_test_label <- c(all_test_label,as.vector(testLabel))
  all_pred_label <- c(all_pred_label,as.vector(eset_pred))

  svm_acc <- c(svm_acc,mean(all_test_label == all_pred_label))
}

library(gmodels) # CrossTable
CrossTable(x=all_test_label,y=all_pred_label, prop.chisq=FALSE)

parametertypes <- c()
for(k in 1:dim(parameters)[1]){
  costtmp <- parameters[k,"cost"]
  costtmp <- paste("cost:",costtmp,sep = "")
  gammatmp <- parameters[k,"gamma"]
  gammatmp <- paste("gamma:",gammatmp,sep = "")
  parametertmp <- paste(costtmp,gammatmp)
  parametertypes <- c(parametertypes,parametertmp)
}
# View(parametertypes)
names(svm_acc) <- parametertypes
svm_acc<- data.frame(svm_acc)
# View(svm_acc)
library(ggplot2)
library(reshape2)
svm_melt<- melt(svm_acc)
colnames(svm_melt) <- c("Parameter","Accuracy")
svm_melt$Accuracy <- round(svm_melt$Accuracy,3)
# pdf(file = "nnet-RFE+SVM不同参数准确度图.pdf",width = 7,height = 5)
ggplot(data = svm_melt,aes(x = Parameter,y = Accuracy,fill = Parameter))+
  geom_bar(stat = 'identity', width = 0.6)+
  geom_text(aes(label = Accuracy),vjust=-0.5)+
  labs(title = "Accuracy of NNET in different parameter") +
  theme(axis.text.x = element_text(angle=30,size=10))
# dev.off()


# 画图 ----------------------------------------------------------------------

trainingLabel <- y1[trainingRows]
testLabel <- y1[-trainingRows]

svmfit <- svm (trainingData,trainingLabel, kernel = "radial", cost = 10, gamma=1e-1, scale = FALSE) # radial svm, scaling turned OFF
print(svmfit)
pred <- predict(svmfit, testData, decision.values = TRUE)
compareTable <- table(testLabel, pred)  # comparison table
compareTable

library(caret)
confusionMatrix(compareTable)

library(pROC)
# pdf(file = "ROC_SVM_knn.pdf",width = 5,height = 5)
plot.roc(testLabel, pred, print.auc=T, main="pAUC")
# Roc_svm <- smooth(roc, method="binormal")
# legend("bottomright", legend=c("Accuracy=0.942", "Precision=0.750 ", "Sensitivity=0.360", "Specificity=0.990", "F-measure=0.486"))
# dev.off()



# 性能指标 --------------------------------------------------------------------

predict = ifelse(pred > 0.5, 1, 0)
predict_value = predict 
true_value = testLabel
error = predict_value-true_value
table(true_value, predict_value) 

data <- x
accuracy = (nrow(data)-sum(abs(error)))/nrow(data) 
precision = sum(true_value & predict_value)/sum(predict_value)  #真实值预测值全为1 / 预测值全为1 --- 提取出的正确信息条数/提取出的信息条数
recall = sum(predict_value & true_value)/sum(true_value)        #真实值预测值全为1 / 真实值全为1 --- 提取出的正确信息条数 /样本中的信息条数
# P和R指标有时候会出现的矛盾的情况，这样就需要综合考虑他们，最常见的方法就是F-Measure（又称为F-Score）
F_measure= 2*precision*recall/(precision+recall)    #F-Measure是Precision和Recall加权调和平均，是一个综合评价指标
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value)) 


# 保存结果 --------------------------------------------------------------------
A_roc <- cbind(testLabel, pred)
result_svm <- c(accuracy, precision, recall, specificity, F_measure)
View(result_svm)
# write.csv(A_roc, file = "A_roc_nnet.csv", row.names = F)
# write.csv(result_svm, file = "result_ab_468-571.csv", row.names = F)
# write.csv(dif_feature, file = "dif_feature_nnet.csv")
# write.csv(gene_ref, file = "gene_rfe_nnet.csv", row.names = F)







