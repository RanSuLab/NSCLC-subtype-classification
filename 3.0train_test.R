library("caret")
library("e1071")
library("ROCR")
library("randomForest")

set.seed(123)
for(i in 1:1){
genelistWGRFE<-read.csv('gene_WGRFE.txt',sep = '\t',header = FALSE)
genelistWGRFE<-as.vector(factor(genelistWGRFE$V1))
genelistWGRFE[length(genelistWGRFE)+1]<-"class"
genelistWGRFE2<-read.csv('gene_WGRFE2.txt',sep = '\t',header = FALSE)
genelistWGRFE2<-as.vector(factor(genelistWGRFE2$V1))
genelistWGRFE2[length(genelistWGRFE2)+1]<-"class"
genelistWGRFE3<-read.csv('gene_WGRFE3.txt',sep = '\t',header = FALSE)
genelistWGRFE3<-as.vector(factor(genelistWGRFE3$V1))
genelistWGRFE3[length(genelistWGRFE3)+1]<-"class"
}

#The genes selected by WGRFE after voting use SVM and RF classifiers respectively
for(k in 1:1){
result_WGRFEsvm<-list()
result_WGRFErf<-list()
for(i in 1:10){
  cat(i)
  train3<-train[[i]]
  train3<-train3[,genelistWGRFE]
  test3<-test[[i]]
  test3<-test3[,genelistWGRFE]
  tempsvm<-test_svm(train3,test3)
  temprf<-test_rf(train3,test3)
  result_WGRFEsvm[[i]]=tempsvm
  result_WGRFErf[[i]]=temprf
}
rm(tempsvm)
rm(temprf)
rm(train3)
rm(test3)
}


#The genes selected by WGRFE before voting use SVM and RF classifiers respectively
for(k in 1:1){
  result_WGRFEsvm2<-list() #vote berfor 
  result_WGRFErf2<-list()  #vote before
  for(i in 1:10){
    cat(i)
    train3<-train[[i]]
    train3<-train3[,genelistWGRFE3]
    test3<-test[[i]]
    test3<-test3[,genelistWGRFE3]
    tempsvm<-test_svm2(train3,test3)
    result_WGRFEsvm2[[i]]=tempsvm

  }
  for(i in 1:10){
    cat(i)
    train3<-train[[i]]
    train3<-train3[,genelistWGRFE2]
    test3<-test[[i]]
    test3<-test3[,genelistWGRFE2]
    temprf<-test_rf2(train3,test3)
    result_WGRFErf2[[i]]=temprf
  }
  rm(tempsvm)
  rm(temprf)
  rm(train3)
  rm(test3)
}

for(k in 1:1){
  ass_votedWGRFE_SVM<-vector()
  ass_votedWGRFE_RF<-vector()
  ass_WGRFE_SVM<-vector()
  ass_WGRFE_RF<-vector()
for(i in 1:10){
  ass_votedWGRFE_SVM=rbind(ass_votedWGRFE_SVM,result_WGRFEsvm[[i]][[1]])
  ass_votedWGRFE_RF=rbind(ass_votedWGRFE_RF,result_WGRFErf[[i]][[1]])
  ass_WGRFE_SVM=rbind(ass_WGRFE_SVM,result_WGRFEsvm2[[i]][[1]])
  ass_WGRFE_RF=rbind(ass_WGRFE_RF,result_WGRFErf2[[i]][[1]])

  }
}
for(k in 1:1){
mean(ass_WGRFE_SVM[,4])
mean(ass_WGRFE_RF[,4])
mean(ass_votedWGRFE_SVM[,4])
mean(ass_votedWGRFE_RF[,4])
mean(ass_RFE_SVM[,4])
mean(ass_RFE_RF[,4])
}







test_svm = function(train,test){
  result2<-list()
  result<-vector()
  traindata=train
  testdata=test
  obj<-tune.svm(class~.,data=traindata,gamma = 10^0,cost = 10^(-4:3),kernel = "linear") #gamma=-3
  best_gm <- obj$best.parameters$gamma
  best_ct <- obj$best.parameters$cost
  model <- svm(class~.,data=traindata,
               kernel = "linear", gamma = best_gm, cost = best_ct, cross = 0,probability = TRUE)
  
  pre_label <- predict(model, testdata, probability = TRUE)
  
  stat_res <- table(pre_label, testdata$class)
  accurary <- (stat_res[1, 1] + stat_res[2, 2])/length(pre_label)
  sensitivity <- stat_res[2, 2]/(stat_res[1, 2] + stat_res[2, 2])
  specificity <- stat_res[1, 1]/(stat_res[1, 1] + stat_res[2, 1])
  pre_prob <- attr(pre_label, "probabilities")
  pred <- prediction(pre_prob[, 2], testdata$class)
  perf <- performance(pred, "tpr", "fpr")
  auc.tmp <- performance(pred, "auc")
  auc.value <- auc.tmp@y.values
  AUC <- round(auc.value[[1]], digits = 4)
  result<-cbind(accurary,sensitivity,specificity,AUC)
  result2<-list(result,perf,AUC)
  return(result2)
  }

test_rf = function(train,test){
  result2<-list()
  result<-vector()
  traindata=train
  testdata=test3
  test$class<-factor(test$class)
  obj<-tune.randomForest(class~.,data=traindata,mtry=(6:11))
  best_mtry <- obj$best.parameters$mtry
  model <- randomForest(class~.,data=traindata,
                        mtry = best_mtry,cross = 0,proximity=TRUE)
  
  pre_label <- predict(model, testdata)
  
  stat_res <- table(pre_label, testdata$class)
  accurary <- (stat_res[1, 1] + stat_res[2, 2])/length(pre_label)
  sensitivity <- stat_res[2, 2]/(stat_res[1, 2] + stat_res[2, 2])
  specificity <- stat_res[1, 1]/(stat_res[1, 1] + stat_res[2, 1])
  
  pre_label2<- predict(model, testdata,type="prob")
  predictions=as.vector(pre_label2[,2])
  pred=prediction(predictions,testdata$class)
  perf <- performance(pred, "tpr", "fpr")
  #plot(perf,colorize=TRUE)
  auc.tmp <- performance(pred, "auc")
  auc.value <- auc.tmp@y.values
  AUC <- round(auc.value[[1]], digits = 4)
  result<-cbind(accurary,sensitivity,specificity,AUC)
  result2<-list(result,perf,AUC)
  return(result2) 
}



test_svm2 = function(train,test){
  result2<-list()
  result<-vector()
  traindata=train
  testdata=test
  obj<-tune.svm(class~.,data=traindata,gamma = 10^0,cost = 10^(-1:2),kernel = "linear") #gamma=-3
  best_gm <- obj$best.parameters$gamma
  best_ct <- obj$best.parameters$cost
  model <- svm(class~.,data=traindata,
               kernel = "linear", gamma = best_gm, cost = best_ct, cross = 0,probability = TRUE)
  
  pre_label <- predict(model, testdata, probability = TRUE)
  
  stat_res <- table(pre_label, testdata$class)
  accurary <- (stat_res[1, 1] + stat_res[2, 2])/length(pre_label)
  sensitivity <- stat_res[2, 2]/(stat_res[1, 2] + stat_res[2, 2])
  specificity <- stat_res[1, 1]/(stat_res[1, 1] + stat_res[2, 1])
  pre_prob <- attr(pre_label, "probabilities")
  pred <- prediction(pre_prob[, 2], testdata$class)
  perf <- performance(pred, "tpr", "fpr")
  auc.tmp <- performance(pred, "auc")
  auc.value <- auc.tmp@y.values
  AUC <- round(auc.value[[1]], digits = 4)
  result<-cbind(accurary,sensitivity,specificity,AUC)
  result2<-list(result,perf,AUC)
  return(result2)
}

test_rf2 = function(train,test){
  result2<-list()
  result<-vector()
  traindata=train
  testdata=test3
  test$class<-factor(test$class)
  obj<-tune.randomForest(class~.,data=traindata,mtry=(1:2))
  best_mtry <- obj$best.parameters$mtry
  model <- randomForest(class~.,data=traindata,
                        mtry = best_mtry,cross = 0,proximity=TRUE)
  
  pre_label <- predict(model, testdata)
  
  stat_res <- table(pre_label, testdata$class)
  accurary <- (stat_res[1, 1] + stat_res[2, 2])/length(pre_label)
  sensitivity <- stat_res[2, 2]/(stat_res[1, 2] + stat_res[2, 2])
  specificity <- stat_res[1, 1]/(stat_res[1, 1] + stat_res[2, 1])
  
  pre_label2<- predict(model, testdata,type="prob")
  predictions=as.vector(pre_label2[,2])
  pred=prediction(predictions,testdata$class)
  perf <- performance(pred, "tpr", "fpr")
  #plot(perf,colorize=TRUE)
  auc.tmp <- performance(pred, "auc")
  auc.value <- auc.tmp@y.values
  AUC <- round(auc.value[[1]], digits = 4)
  result<-cbind(accurary,sensitivity,specificity,AUC)
  result2<-list(result,perf,AUC)
  return(result2) 
}












