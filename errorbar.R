#先加载包
library(ggplot2)
library(ggpubr)
library(dplyr)  


rm(list=ls())

View(ToothGrowth)
class(ToothGrowth)

p <- ggboxplot(ToothGrowth, x = "supp", y = "len",
               color = "supp", palette = "jco",
               add = "jitter")
# Add p-value
p + stat_compare_means()
# Change method
p + stat_compare_means(method = "t.test")
#####################################################################


rm(list=ls())

setwd('D:\\E\\博士\\R_程序\\HGSOC\\Data')
data <- read.table("feature__expr_data.txt", header = T, check.names = FALSE)
pred_log <- data.frame(t(data))
dim(pred_log)

# Compare two independent groups ------------------------------------------
k = 1
p <- ggboxplot(pred_log, x = "Label", y = gene[k],
               color = "Label", palette = "jco",
               add = "jitter")
## Add p-value
p + stat_compare_means()
p + stat_compare_means(method = "t.test")


# 存储显著性 -------------------------------------------------------------------

annotation=c("3.5e-07", "4.7e-05", "3.7e-11", "4.6e-15", "4.0e-09",
             "9.3e-06", "3.6e-05", "2.9e-06", "5.2e-06", "1.4e-06",
             "3.9e-07", "3.1e-07", "7.7e-07", "1.3e-05", "2.1e-07",
             "1.8e-08", "1.6e-08", "4.4e-06", "1.4e-06", "1.4e-04",
             "1.4e-04", "5.7e-05", "6.5e-05", "1.5e-05", "2.8e-07",
             "8.1e-10", "1.7e-06", "3.6e-07", "4.8e-06", "1.6e-05",
             "8.3e-05", "1.2e-06", "4.8e-05", "1.4e-05", "9.0e-05",
             "8.5e-06", "5.9e-05", "4.7e-05", "9.4e-08", "1.4e-04",
             "3.7e-05", "8.3e-05", "6.7e-08", "1.1e-04", "4.6e-05",
             "1.8e-08", "3.8e-05", "9.1e-05", "1.7e-04", "4.6e-05",
             "5.1e-05", "3.8e-06")

# annotation=c("****", "****","****","****","****","****", "****","****","****","****",
#              "****", "****","****","****","****","****", "****","****","****","***",
#              "***",  "****","****","****","****","****", "****","****","****","****",
#              "****", "****","****","****","****","****", "****","****","****","***",
#              "****", "****","****","***", "****","****", "****","****","***", "****",
#              "****","****")

## 固定显著性的位置
y <- c()
y[1] = max(max(pred_log[1:10,gene[1]]),max(pred_log[11:20,gene[1]]) ) +0.5

for (l in 2:52) {
  y[l] = max(max(pred_log[1:10,gene[l]]),max(pred_log[11:20,gene[l]]) ) +0.5
}

for (l in 2:52) {
  y[l] = max(max(pred_log[1:10,gene[l]]),max(pred_log[11:20,gene[l]]) ) +0.5
}

# 新的尝试 --------------------------------------------------------------------
# rm(list=ls())
library(data.table)
library(ggplot2)
library(ggsignif)


# 变成长格式的数据： ---------------------------------------------------------------

b <- melt(pred_log, id.vars = c("Label"))
b$Label <- as.factor(b$Label)


# 做每个箱线图的均值点及均值连线，需获得每个组每个属性的均值，并定义每个组每个属性的X坐标为固定值 ------------------

# group1的mean:
c<-copy(b)
setDF(c)
c1<-tapply(c[c$Label==1,"value"],c[c$Label==1,"variable"],mean)
c2<-tapply(c[c$Label==0,"value"],c[c$Label==0,"variable"],mean)
c3<-rbind(data.frame(variable=names(c1),value=c1,Label=1),data.frame(variable=names(c2),value=c2,Label=0))

c3$Label<-as.factor(c3$Label)
# 分别计算两组均值，用来画折线图：
c3$variable2<-NA

gene <- c3$variable 
lab1 <- c()
lab2 <- c()

## lab1,lab2与点的左右位置有关，需要调节
lab1[1] = 1.2
lab2[1] = 0.8
c3[c3$Label==1&c3$variable==gene[1],"variable2"]<-lab1[1]
c3[c3$Label==0&c3$variable==gene[1],"variable2"]<-lab2[1]

data=data.frame(x=lab2,
                xend=lab1,
                y=y,
                annotation=annotation)


for (i in 2:52) {
  # i = 2
  lab1[i] <- lab1[i-1] + 1
  lab2[i] <- lab2[i-1] + 1
  c3[c3$Label==1&c3$variable==gene[i],"variable2"]<-lab1[i]
  c3[c3$Label==0&c3$variable==gene[i],"variable2"]<-lab2[i]
}

p1<-ggplot(b)+
  geom_boxplot(aes(x=variable,y=value,fill=Label),width=0.6,
               position = position_dodge(0.8),outlier.size = 0,outlier.color = "white")+
  scale_fill_manual(values = c("red", "blue"),breaks=c("1","0"),labels=c("normal","HGSOC"))+
  geom_point(data=c3,aes(x=variable2,y=value,color=Label),shape=15,size=1)+
  geom_line(data=c3,aes(x=variable2,y=value,color=Label),size=1,linetype = "dotted")+
  # geom_smooth(data=c3,method = 'loess',formula = 'y ~ x',
  #             aes(x=variable2,y=value,color=Label),size=1,linetype = "dashed")+
  # stat_summary(fun = mean, geom = "errorbar",
  #              aes(x=variable,y=value,ymax = ..y.., ymin = ..y..,color=Label),
  #              width = .75, linetype = "dashed") +
  xlab("")+
  ylab("")+
  scale_y_continuous(limits = c(0,15),breaks=seq(0,15,5)) +
  geom_signif(stat="identity",
              data=data.frame(x=lab2,
                              xend=lab1,
                              y=y,
                              annotation=annotation),
              aes(x=x,xend=xend, y=y, yend=y, 
                  annotation=annotation, hjust=0, vjust=0.5, angle=90)) +
  xlab("Biomarker") +
  ylab("Expression") +
  theme_bw()+
  theme(
    legend.position = "top",
    legend.background=element_blank(),
    legend.key = element_blank(),
    legend.margin=margin(0,0,0,0,"mm"),
    axis.text.x=element_text(size=rel(1.1),face="bold",angle=90),
    axis.line.x = element_line(size = 0.5, colour = "black"),
    axis.line.y = element_line(size = 0.5, colour = "black"),
    legend.text=element_text(size=rel(1.1)),
    legend.title=element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(color=FALSE)
p1

# pdf(file = "p1_HGSOC1.pdf",width = 8,height = 6)
# p1
# dev.off()



# 两组整体水平的箱线图： -------------------------------------------------------------

mean(as.numeric(annotation))

p2<-ggplot(b)+
  geom_boxplot(aes(x=Label,y=value,fill=Label),width=0.8,position=position_dodge(1))+
  stat_summary(fun = mean, geom = "point", aes(x=Label,y=value,color=Label),shape=15)+
  scale_fill_manual(values = c("red", "blue"),breaks=c("1","0"),labels=c("normal","HGSOC"))+
  scale_x_discrete(breaks=c("1","0"),labels=c("normal","HGSOC"))+
  scale_y_continuous(limits = c(0,15),breaks=seq(0,15,5)) +
  geom_signif(stat="identity",
              data=data.frame(x=c(1), xend=c(2),
                              y=c(13), annotation=c("****")),
              aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))+
  theme_bw()+
  theme(
    legend.position = "top",
    legend.background=element_blank(),
    legend.key = element_blank(),
    legend.margin=margin(0,0,0,0,"mm"),
    axis.text.x=element_text(size=rel(1.1),face="bold"),
    axis.line.x = element_line(size = 0.5, colour = "black"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    legend.text=element_text(size=rel(1.1)),
    legend.title=element_blank(),
    plot.margin = margin(11.5,0,7,0,"mm"),
    panel.border = element_blank(),
    panel.grid = element_blank()
  )+
  guides(fill=F,color=F)
p2
# pdf(file = "p2_HGSOC1.pdf",width = 2,height = 4)
# p2
# dev.off()

# 合并两个图像： -----------------------------------------------------------------

# 如果没有安装则install.packages("gridExtra")
library(gridExtra)
# 合并两个图：
grid.arrange(p1, p2, nrow=1, ncol=2, widths=c(3.5,1),heights=c(4))



