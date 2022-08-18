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
meta$class <- meta$Group
meta[meta$Group == "T2DM",]$Group <- "CTRL"
meta[meta$Group == "T2DM_CRC",]$Group <- "CRC"

ootest <- function(target, testgroup) {
  testgroupname <- ifelse(length(testgroup)>1, "TOTAL", ifelse(testgroup == 1, "CTRL", "T2DM"))
  mat <- read.table(paste("data/dm_adj_", target, ".txt", sep=""))
  feature <- read.csv(paste("result/7_boruta/", target, "/boruta_feature_", testgroupname, ".csv", sep = ''))
  
  feature <- feature[,1]
  feature[grep(feature, pattern = "\\.")] <- gsub(pattern = "\\.", replacement = "", feature[grep(feature, pattern = "\\.")])
  
  feat.sig <- t(mat[, rownames(meta)])
  # feature %in% colnames(data.sig)
  colnames(feat.sig) <- gsub(pattern = '\\[|\\]|-|\\(|\\)|\\:|\\;|\\/|\\ |\\,|\\.', replacement = "", colnames(feat.sig))
  # print(feature == colnames(feat.sig[,feature]))
  
  data.sig <- cbind(feat.sig[, feature], data.frame(meta[rownames(feat.sig), 'Group']))
  colnames(data.sig)[ncol(data.sig)] <- 'Group'
  
  idy <- meta$Batch %in% testgroup # Train
  !idy # Validate
  
  library(caret)
  # repforest <- lapply(1:10, function(repeats){
  set.seed(222)
  colnames(data.sig)
  
  data <- data.sig[idy, ]
  data$Group <- factor(data$Group)
  validate.data <- data.sig[!idy,]
  # Group <- Group[-validateFold[[1]]]

  control <- trainControl(method="repeatedcv", number = 5, classProbs=T, summaryFunction=twoClassSummary)
  rf_default <- train(Group~., data=data,
                      metric="ROC",
                      # preProcess = c("center", "scale"),
                      trControl=control)
  rf_default
  best_rf <- rf_default$bestTune
  set.seed(2021)
  fold <- createFolds(y = data$Group, k=5)
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
  
  # Group <- fold_train_best$Group
  best_model <- randomForest(Group~., data=fold_train_best,mtry=best_rf$mtry,
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

spe.ctrl <- ootest("species", 1)
spe.t2dm <- ootest("species", 2)


pdf(file = paste("result/8_oo_rf/species.rf.pdf", sep = ''))
auc.spe.ctrl <- ROC_Test_Validate(spe.ctrl, 1, start.plot = TRUE)
auc.spe.t2dm <- ROC_Test_Validate(spe.t2dm, 2)
legend("bottomright", c(auc.spe.ctrl, auc.spe.t2dm), 
       lwd = 2, lty=1,
       col = colvec, bty="n")
dev.off()
  

ko.ctrl <- ootest("ko", 1)
ko.t2dm <- ootest("ko", 2)

pdf(file = paste("result/8_oo_rf/ko.rf.pdf", sep = ''))
auc.ko.ctrl <- ROC_Test_Validate(ko.ctrl, 1, start.plot = TRUE)
auc.ko.t2dm <- ROC_Test_Validate(ko.t2dm, 2)
# auc.ko.total <- ROC_Test_Validate(ko.total, c(1,2))
legend("bottomright", c(auc.ko.ctrl, auc.ko.t2dm), 
       lwd = 2, lty=1,
       col = colvec, bty="n")
dev.off()



  