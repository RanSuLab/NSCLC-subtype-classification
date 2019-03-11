library("caret")
library("e1071")
library("ROCR")
library("randomForest")

set.seed(123)
for(i in 1:1){
genelistRFE<-read.csv('gene_RFE.txt',sep = '\t',header = FALSE)
genelistRFE<-as.vector(factor(genelistRFE$V1))
genelistRFE[length(genelistRFE)+1]<-"class"
genelistRFE2<-read.csv('gene_RFE2.txt',sep = '\t',header = FALSE)
genelistRFE2<-as.vector(factor(genelistRFE2$V1))
genelistRFE2[length(genelistRFE2)+1]<-"class"
}

#RFE-SVM and RFE-RF
for(k in 1:1){
  result_RFEsvm<-list() #vote berfor 
  result_RFErf<-list()  #vote before
  for(i in 1:10){
    cat(i)
    train3<-train[[i]]
    train3<-train3[,genelistRFE]
    test3<-test[[i]]
    test3<-test3[,genelistRFE]
    tempsvm<-test_svm2(train3,test3)
    result_RFEsvm[[i]]=tempsvm

  }
  for(i in 1:10){
    cat(i)
    train3<-train[[i]]
    train3<-train3[,genelistRFE2]
    test3<-test[[i]]
    test3<-test3[,genelistRFE2]
    temprf<-test_rf(train3,test3)
    result_RFErf[[i]]=temprf
  }
  rm(tempsvm)
  rm(temprf)
  rm(train3)
  rm(test3)
}

for(k in 1:1){
  ass_RFE_SVM<-vector()
  ass_RFE_RF<-vector()
for(i in 1:10){
  ass_RFE_SVM=rbind(ass_RFE_SVM,result_RFEsvm[[i]][[1]])
  ass_RFE_RF=rbind(ass_RFE_RF,result_RFErf[[i]][[1]])

  }
}





for(k in 1:1){
write.csv(ass_votedWGRFE_SVM[,1:3],"ass_votedWGRFE_SVM.csv")
write.csv(ass_votedWGRFE_RF[,1:3],"ass_votedWGRFE_RF.csv")
write.csv(ass_WGRFE_SVM[,1:3],"ass_WGRFE_SVM.csv")
write.csv(ass_WGRFE_RF[,1:3],"ass_WGRFE_RF.csv")
write.csv(ass_RFE_SVM[,1:3],"ass_RFE_SVM.csv")
write.csv(ass_RFE_RF[,1:3],"ass_RFE_RF.csv")
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
  obj<-tune.svm(class~.,data=traindata,gamma = 10^0,cost = 10^(-2:2),kernel = "linear") #gamma=-3
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










