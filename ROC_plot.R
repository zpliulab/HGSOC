setwd("D:\\E\\博士\\R_程序\\HGSOC\\Data")
# load(".RData")    # 载入工作区


library(pROC)

A_roc_svm <- read.csv("A_test_120196.csv", header=TRUE, sep = ',')
A_roc_knn <- read.csv("A_test_26712.csv", header=TRUE, sep = ',')
A_roc_rf <- read.csv("A_test_40595_38.csv", header=TRUE, sep = ',')
A_roc_ab <- read.csv("A_test_54388.csv", header=TRUE, sep = ',')
A_roc_nnet <- read.csv("A_test_27651.csv", header=TRUE, sep = ',')

# pdf(file = "ROC.pdf",width = 5,height = 5)
plot.roc(A_roc_svm[,2], A_roc_svm[,1])
roc_knn <- lines.roc( A_roc_knn[,2], A_roc_knn[,1], col="2" )
roc_rf <- lines.roc( A_roc_rf[,2], A_roc_rf[,1], col="3" )
roc_ab <- lines.roc( A_roc_ab[,2], A_roc_ab[,1], col="6" )
roc_nnet <- lines.roc( A_roc_nnet[,2], A_roc_nnet[,1], col="4" )

legend("bottomright", legend=c("GSE120196", "GSE26712", "GSE40595", "GSE54388", "GSE27651"), col=c("1", "2", "3", "6", "4"), lwd=2)
# dev.off()
