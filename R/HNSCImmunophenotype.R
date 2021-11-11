#' @title HNSCImmunophenotype:A package for calculating HNSC immunophenotype.
#'
#' @docType package
#' @name HNSCImmunophenotype
NULL

#' HNSCImmunophenotype data
#'
#' @format list:
#' \describe{
#'   \item{geo.neg.gene}{the expression matrix of the GSE65858 HPV-}
#'   \item{geo.neg.type}{the immunophenotype of the GSE65858 HPV-}
#'   \item{geo.pos.gene}{the expression matrix of the GSE65858 HPV+}
#'   \item{geo.pos.type}{the immunophenotype of the GSE65858 HPV+}
#'   \item{tcga.neg.gene}{the expression matrix of the TCGA-HNSC HPV-}
#'   \item{tcga.neg.type}{the immunophenotype of the TCGA-HNSC HPV-}
#'   \item{tcga.pos.gene}{the expression matrix of the TCGA-HNSC HPV+}
#'   \item{tcga.pos.type}{the expression matrix of the TCGA-HNSC HPV+}
#' }
"dataset"


#' The HNSC immunophenotypev analysis pipeline
#'
#' \code{HNSCImmunophenotype} Returns the HNSC immunophenotype.
#'
#' @param exp the gene expression data set. A matrix with row names as symbols and columns as samples.
#' @param hpv the HPV infection information,There are "positive" and "negative" can been choiced.
#' @param ntree The number of decision trees contained in the random forest,defaults = 10000.
#' @param set.seed Set random seeds so that the results can be repeated,defaults = 1234.
#'
#' @return the immunophenotype.
#' @export
#' @importFrom caret createDataPartition
#' @importFrom randomForest randomForest
#' @importFrom pROC plot.roc
#' @importFrom FactoMineR PCA
#' @importFrom factoextra fviz_pca_ind
#' @importFrom pheatmap pheatmap
#'
#' @examples


HNSCImmunophenotype <- function(exp,hpv,ntree=10000,set.seed=123){
  ifelse(hpv %in% c('negative','positive'),gene <- exp,
  print('please enter HPV infection information'))
  #####*****************载入数据开始质控************#
  gene<- as.matrix(gene)
  gene <- na.omit(gene)
  ###data distribution
  #输入数据检测，行为基因，列为样本名，输出均值、方差和中位数绝对偏差(MAD)分布。
  gene <- scale(gene)
  gene = sweep(gene,1, apply(gene,1,median))
  rownames(gene) <- gsub('-','\\.',rownames(gene))

  ##**********读入数据*****##
  ifelse(hpv=='negative',tcga.exp <- dataset$tcga.neg.gene,
         ifelse(hpv=='positive',tcga.exp <- dataset$tcga.pos.gene,
                print('please enter HPV infection information')))

  ifelse(hpv=='negative',tcga.type <- dataset$tcga.neg.type,
         ifelse(hpv=='positive',tcga.type <- dataset$tcga.pos.type,
                print('please enter HPV infection information')))


  ifelse(hpv=='negative',geo.exp <- dataset$geo.neg.gene,
         ifelse(hpv=='positive',geo.exp <- dataset$geo.pos.gene,
                print('please enter HPV infection information')))

  ifelse(hpv=='negative',geo.type <- dataset$geo.neg.type,
         ifelse(hpv=='positive',geo.type <- dataset$geo.pos.type,
                print('please enter HPV infection information')))

  rownames(tcga.exp) <- gsub('-','\\.',rownames(tcga.exp))
  rownames(geo.exp) <- gsub('-','\\.',rownames(geo.exp))

  ######************取交集******###
  co <- intersect(rownames(gene),rownames(tcga.exp))
  model.gene <- as.data.frame(co)
  colnames(model.gene) <- 'model.gene'
  gene <- gene[co,]
  tcga.exp <- tcga.exp[co,]
  print(paste('there are ',length(co),' gene are used to build the model'))
  tcga.exp <- t(tcga.exp) %>% as.data.frame()
  tcga.exp <- tcga.exp[rownames(tcga.type),]
  tcga.exp$group <- tcga.type[,1]
  tcga.exp$group <- factor(tcga.exp$group)

  geo.exp <- t(geo.exp) %>% as.data.frame()
  geo.exp <- geo.exp[rownames(geo.type),]
  geo.exp$group <- geo.type[,1]
  geo.exp$group <- factor(geo.exp$group)

  gene <- t(gene) %>% as.data.frame()


  colnames(tcga.exp) <- gsub('-','\\.',colnames(tcga.exp))
  colnames(geo.exp) <- gsub('-','\\.',colnames(geo.exp))
  colnames(gene) <- gsub('-','\\.',colnames(gene))

  #####     随机森林分类
  #随机抽样，分组，按7:3分为训练组和验证组
  set.seed(set.seed)
  co <- createDataPartition(tcga.exp$group,p=0.7,list = F)
  train <- tcga.exp[co,]
  test1 <- tcga.exp[-co,]
  test2 <- geo.exp
  test3 <- gene
  set.seed(set.seed)
  train.forest <- randomForest(train$group ~ .,data = train,
                               importance=T,proximity=T,ntree=ntree)
  pdf('./mytree.pdf')
  plot(train.forest)
  dev.off()
  #训练集自身测试
  #训练集准确率
  my.predict.train <- predict(train.forest, train[,1:(ncol(train)-1)])
  df<- table(my.predict.train, train$group)
  predict.train <- sum(diag(df)/sum(df))

  ##使用测试集评估
  my.predict1 <- predict(train.forest, test1[,1:(ncol(test1)-1)])
  df<- table(my.predict1, test1$group)
  predict.vulue1 <- sum(diag(df)/sum(df))

  my.predict2 <- predict(train.forest, test2)
  df2<- table(my.predict2, test2$group)
  predict.vulue2 <- sum(diag(df2)/sum(df2))

  my.predict3 <- predict(train.forest, test3)
  gene$group <- my.predict3

  #构建ROC模型
  roc1 <- plot.roc(as.ordered(train$group),as.ordered(my.predict.train),
                   col='red',lwd=3,#曲线颜色、线条类型、粗细
                   print.auc=T,
                   identity.lty=2,identity.lwd=3,#对角线线条及粗细
  )
  roc2 <- plot.roc(as.ordered(test1$group),as.ordered(my.predict1),
                   col='red',lwd=3,#曲线颜色、线条类型、粗细
                   print.auc=T,
                   identity.lty=2,identity.lwd=3,#对角线线条及粗细
  )

  roc3 <- plot.roc(as.ordered(test2$group),as.ordered(my.predict2),
                   col='red',lwd=3,#曲线颜色、线条类型、粗细
                   print.auc=T,
                   identity.lty=2,identity.lwd=3,#对角线线条及粗细
  )

  n1 <- round(roc1$auc,2)
  n2 <- round(roc2$auc,2)
  n3 <- round(roc3$auc,2)

  #绘制ROC曲线
  pdf('./ROC.pdf',
      width = 4,height = 4)
  plot.roc(roc1,col="#FF6665", identity.col ="grey",lwd=4,identity.lty=3,
           identity.lwd=3)
  plot.roc(roc2,col="#356567",lwd=3,add = T)
  plot.roc(roc3,col="#659BCD",lwd=3,add = T)
  legend("bottomright",
         c(paste('TCGA_train,AUC=',n1),paste('TCGA_test,AUC=',n2),
           paste('GSE65858_test,AUC=',n3)),
         col = c("#FF6665","#356567","#659BCD"),
         lwd = 6,bty='n')
  dev.off()


  df <- gene
  df$group <- factor(df$group)
  df <- df[order(df$group),]

  # 热图
  annotation_col <- as.data.frame(df[,c('group')])
  exp <- df[,1:(ncol(df)-1)] %>% t()

  row.names(annotation_col) <- colnames(exp)
  colnames(annotation_col) <- 'subtype'

  k= round(nrow(exp)/ncol(exp),1)
  pheatmap(exp,
           method = "spearman",
           cluster_rows=T, cluster_cols=F,
           cellwidth = k, cellheight = 3, fontsize = 3,
           color = colorRampPalette(c("navy","white","red"))(64),
           scale="row",show_colnames = F,
           border_color = "NA",
           annotation_col = annotation_col,
           filename = "./pheatmap.pdf")

  write.csv(gene,'./predict.result.csv')

  #######*******PCA降维和绘图*****###
  result <- FactoMineR::PCA(df[,1:(ncol(df)-1)])
  dev.off()
  ind <- factoextra::fviz_pca_ind(
    result,
    habillage = factor(df$group),
    label = "none",
    mean.point = F
  )
  # pdf('./PCA1.pdf')
  ind
  # dev.off()
}



