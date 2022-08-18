library(Boruta)
library(corrplot)
library(dplyr)
library(data.table)
library(randomForest)
require(mlbench)
require(caret)
library(pROC)

setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202201DCRC")

boruta.imp <- function(x){
  imp <- reshape2::melt(x$ImpHistory, na.rm=T)[,-1]
  colnames(imp) <- c("Variable","Importance")
  imp <- imp[is.finite(imp$Importance),]
  
  variableGrp <- data.frame(Variable=names(x$finalDecision),
                            finalDecision=x$finalDecision)
  
  showGrp <- data.frame(Variable=c("shadowMax", "shadowMean", "shadowMin"),
                        finalDecision=c("shadowMax", "shadowMean", "shadowMin"))
  
  variableGrp <- rbind(variableGrp, showGrp)
  
  boruta.variable.imp <- merge(imp, variableGrp, all.x=T)
  
  sortedVariable <- boruta.variable.imp %>% group_by(Variable) %>%
    summarise(median=median(Importance)) %>% arrange(median)
  sortedVariable <- as.vector(sortedVariable$Variable)
  
  
  boruta.variable.imp$Variable <- factor(boruta.variable.imp$Variable, levels=sortedVariable)
  
  invisible(boruta.variable.imp)
}

# convert_name <- function(x, matname){
#   return(matname[as.numeric(gsub(x, pattern = "a", replacement=""))])
# }

# #####################################
# no convert name
# #####################################
# target
  
target <- "species" # ko
testgroup <- 1 # CTRL

meta <- read.table("data/dm_meta.txt")
meta[meta$Group == "T2DM",]$Group <- "CTRL"
meta[meta$Group == "T2DM_CRC",]$Group <- "CRC"

auto_boruta_rf <- function(target, testgroup) {
  testgroupname <- ifelse(length(testgroup)>1, "TOTAL", ifelse(testgroup == 1, "CTRL", "T2DM"))
  mat <- read.table(paste("data/dm_adj_", target, ".txt", sep=""))
  
  # differential_expression
  diff <- read.csv(paste("result/1_MMUPHin/MMUPHin_", target, "/MMUPHin_maaslin_result.csv", sep = ""))
  # colnames(diff)
  diff <- na.omit(diff)
  rownames(diff) <- diff[,2]
  diff <- diff[,c(8,18)]
  # table(diff[,1] < 0.05 | diff[,2] < 0.05)
  if(length(testgroup)>1) {
    idx <- diff[,1] < 0.05 | diff[,2] < 0.05
  } else {
    idx <- diff[,testgroup] < 0.05
  }
  s.diff <- diff[idx,]
  s.diff <- rownames(s.diff)
  s.diff
  
  # boruta
  idy <- meta$Batch %in% testgroup #选择分类组别
  feat.sig <- t(mat[s.diff, rownames(meta)[idy]])
  
  # filter
  abundance <- 0.01
  table(colSums(feat.sig)/sum(feat.sig)>=(abundance/100))
  feat.sig <- feat.sig[,colSums(feat.sig)/sum(feat.sig)>=(abundance/100)]
  data.sig <- cbind(feat.sig, data.frame(meta[rownames(feat.sig), 'Group']))
  colnames(data.sig)[ncol(data.sig)] <- 'Group'
  Group <- factor(data.sig$Group, levels = )
  
  timestart<-Sys.time()
  set.seed(1)
  boruta <- Boruta(x=feat.sig, y=Group, pValue=0.05, mcAdj=T,
                   maxRuns=1000)
  timeend <- Sys.time()
  runningtime <- timeend-timestart
  print(runningtime)
  print(table(boruta$finalDecision))
  
  #extract feature
  boruta.finalVarsWithTentative <- data.frame(Item=getSelectedAttributes(boruta, withTentative = F), Type="Boruta_with_tentative")
  feature <- boruta.finalVarsWithTentative$Item
  boruta.variable.imp.use <- boruta$ImpHistory[,feature]
  feature_importance <- apply(boruta.variable.imp.use,2,mean)
  feature_importance <- data.frame(sort(feature_importance,decreasing = TRUE))
  feature <- rownames(feature_importance)
  print(feature)
  #
  suppressWarnings(dir.create(paste("result/7_boruta", target, sep = "/")))
  write.csv(feature_importance, file = paste("result/7_boruta/", target, "/boruta_feature_", testgroupname, ".csv", sep = ''))
  
  ##importance
  boruta.variable.imp <- boruta.imp(boruta)
  head(boruta.variable.imp)
  # boruta.variable.imp$Variable<- gsub(pattern = "s__", boruta.variable.imp$Variable, replacement = "")
  # remotes::install_github("Tong-Chen/YSX")
  library(YSX)
  write.csv(boruta.variable.imp, paste("result/7_boruta/", target, "/boruta_feature_imp_", testgroupname, ".csv", sep = ''))
  
  feature_impor_plot <- sp_boxplot(boruta.variable.imp, melted=T, xvariable = "Variable", yvariable = "Importance",
                                   legend_variable = "finalDecision", # legend_variable_order = c("shadowMax", "shadowMean", "shadowMin", "Confirmed"),
                                   legend_variable_order = c("Confirmed"), manual_color_vector = "red",
                                   xtics_angle = 90, coordinate_flip =T,
                                   outlier = TRUE, notch = FALSE)
  feature_impor_plot
  pdf(paste("result/7_boruta/", target, "/boruta_feature_selected_imp_plot_", testgroupname, ".pdf", sep = ''),
      useDingbats = FALSE, width = 6, height = length(feature)/6)
  plot(feature_impor_plot)
  dev.off()
  
  # randomforest
  library(caret)
  # repforest <- lapply(1:10, function(repeats){
  set.seed(222)
  colnames(data.sig)
  
  feature[grep(feature, pattern = "X[0-9]|X\\.")] <- gsub(pattern = "X|X\\.", replacement = "", feature[grep(feature, pattern = "X[0-9]|X\\.")])
  colnames(data.sig) <- gsub(pattern = '\\[|\\]|-|\\(|\\)|\\:|\\;|\\/|\\ |\\,', replacement = ".", colnames(data.sig))
  colnames(data.sig)[substr(colnames(data.sig),1,1) == "."] <- gsub(pattern = "\\.", replacement = "", colnames(data.sig)[substr(colnames(data.sig),1,1) == "."])
  print(feature)
  print(colnames(data.sig[,feature]))
  
  idx <- colnames(data.sig) %in% feature | colnames(data.sig) == "Group"
  validateFold <- createFolds(y = Group, k=5)
  
  data <- data.sig[-validateFold[[1]], idx]
  data$Group <- factor(data$Group)
  validate.data <- data.sig[validateFold[[1]],idx]
  Group <- Group[-validateFold[[1]]]

  control <- trainControl(method="repeatedcv", number = 5, classProbs=T, summaryFunction=twoClassSummary)
  rf_default <- train(as.factor(Group)~., data=data,
                      metric="ROC",
                      # preProcess = c("center", "scale"),
                      trControl=control)
  rf_default
  best_rf <- rf_default$bestTune
  set.seed(2021)
  fold <- createFolds(y = Group, k=5)
  #return(fold)
  metrics <- matrix(NA,5,12)
  colnames(metrics) <- c("Sensitivity","Specificity","Pos Pred Value","Neg Pred Value","Precision",
                         "Recall","F1","Prevalence","Detection Rate","Detection Prevalence",
                         "Balanced Accuracy","AUC")
  library(randomForest)
  library(ROCR)
  library(pROC)
  for(j in 1:5){
    #print(j)
    fold_test <- data[fold[[j]],]
    fold_train <- data[-fold[[j]],]
    print(table(fold_train$Group))
    Group <- fold_train$Group
    set.seed(0518)
    fold_fit <- randomForest(as.factor(Group)~., data=fold_train, mtry=best_rf$mtry,
                             ntree=500, importance=TRUE)
    
    fold_pred <- predict(fold_fit,fold_test)
    #print(fold_pred)
    
    result <- confusionMatrix(factor(as.vector(fold_pred)),
                              as.factor(fold_test$Group),
                              mode = "everything", positive=levels(fold_test$Group)[2]) #######
    
    metrics[j,1:11] <- as.numeric(result$byClass)
    predicted_probs <- predict(fold_fit, fold_test, type = 'prob')
    pred <- prediction(predicted_probs[,2], fold_test$Group)
    auc <- performance(pred, 'auc')
    #print(auc@y.values[[1]])
    metrics[j,12] <- auc@y.values[[1]]
  }
  
  best_index <- which(metrics[,'AUC']==max(metrics[,'AUC']))[1]
  fold_test_best <- data[fold[[best_index]],]
  fold_train_best <- data[-fold[[best_index]],]
  
  Group <- fold_train_best$Group
  best_model <- randomForest(as.factor(Group)~., data=fold_train_best,mtry=best_rf$mtry,
                             ntree=500,importance=TRUE)
  best_model
  result_list <- list("model" = best_model,"metrics" = metrics,"best_fold"=best_index)
  print(result_list)
  
  test.pred.prob <- predict(best_model, fold_test_best, type='prob')
  
  validate.data$Group <- factor(validate.data$Group)
  validate.pred <- predict(best_model, validate.data)
  validate.pred.prob <- predict(best_model, validate.data, type='prob')
  
  return(list(fold_test_best$Group, test.pred.prob[,1], validate.data$Group, validate.pred.prob[,1]))
}

colvec <- c("#16CAB2", "#16697A", "#FF990A", "#FF0022")

ROC_Test_Validate <- function(x, testgroup, start.plot = FALSE) {
  group <- ifelse(length(testgroup)>1, 3, testgroup)
  groupname <- ifelse(length(testgroup)>1, "TOTAL", ifelse(testgroup == 1, "CTRL", "T2DM"))
  tgroup <- x[[1]]
  tprob <- x[[2]]
  vgroup <- x[[3]]
  vprob <- x[[4]]
  validate.roc <- roc(vgroup, vprob) #, levels = c("H_Ob", "C_Ob"), print.auc.y = 40)
  validate.auc <- validate.roc$auc
  if(start.plot) {
    plot(validate.roc, col = colvec[group], lwd = 2,
         lty = 1, xlab = "False Positive Rate", ylab = "True Positive Rate")
  } else {
    plot.roc(vgroup, vprob, legacy.axes = TRUE,
             col = colvec[group], lwd = 2, add = TRUE, lty = 1)
  }
  
  legend.line <- paste(groupname," (AUC=", round(validate.auc,4),")", sep="")
  return(legend.line)
}

spe.ctrl <- auto_boruta_rf("species", 1)
spe.t2dm <- auto_boruta_rf("species", 2)
spe.total <- auto_boruta_rf("species", c(1,2))

ko.ctrl <- auto_boruta_rf("ko", 1)
ko.t2dm <- auto_boruta_rf("ko", 2)
ko.total <- auto_boruta_rf("ko", c(1,2))

pdf(file = paste("result/7_boruta/species/rf.pdf", sep = ''))
auc.spe.ctrl <- ROC_Test_Validate(spe.ctrl, 1, start.plot = TRUE)
auc.spe.t2dm <- ROC_Test_Validate(spe.t2dm, 2)
auc.spe.total <- ROC_Test_Validate(spe.total, c(1,2))
legend("bottomright", c(auc.spe.ctrl, auc.spe.t2dm, auc.spe.total), 
       lwd = 2, lty=1,
       col = colvec[1:3], bty="n")
dev.off()
  
pdf(file = paste("result/7_boruta/ko/rf.pdf", sep = ''))
auc.ko.ctrl <- ROC_Test_Validate(ko.ctrl, 1, start.plot = TRUE)
auc.ko.t2dm <- ROC_Test_Validate(ko.t2dm, 2)
auc.ko.total <- ROC_Test_Validate(ko.total, c(1,2))
legend("bottomright", c(auc.ko.ctrl, auc.ko.t2dm, auc.ko.total), 
       lwd = 2, lty=1,
       col = colvec[1:3], bty="n")
dev.off()



  