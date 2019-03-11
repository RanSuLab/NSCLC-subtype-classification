library("caret")
library("e1071")
library("ROCR")
library("randomForest")
library("WGCNA")

set.seed(123)


for(k in 1:1){
  ##fold_change 
  rankfc=seq(length=2842, from=0, to=0)
  log2fc2<-del_list(log2fc)
  for(i in 1:10){
    rankfc=rankfc+log2fc2[[i]]
    names<-names(rankfc)
  }
  rankfc=abs(rankfc)/10
  rankfc<-sort(rankfc,decreasing = TRUE)
 
  ##WGRFE-SVM ranking
  ranksvm=seq(length=2842, from=0, to=0)
  svm_w<-list()
  for(i in 1:10){
    obj<-tune.svm(class~.,data=train2[[i]],gamma = 10^(0:0),cost = 10^(-2:2),kernel = "linear") #gamma=-3
    best_gm <- obj$best.parameters$gamma
    best_ct <- obj$best.parameters$cost
    stepmodel <- svm(class~.,data=train2[[i]],
                     kernel = "linear", gamma = best_gm, cost = best_ct, cross = 0,probability = TRUE)
    w<-t(find_SVM_weight(stepmodel))
    cat(i)
    svm_w[[i]]=w
    rm(obj)
    rm(stepmodel)
    rm(w)
  }
  for(i in 1:10){
    ranksvm=ranksvm+svm_w[[i]]*R2[[i]]*10000
  }
  ranksvm<-as.vector(abs(ranksvm))
  temp<-train2[[1]]
  names(ranksvm)<-colnames(temp[1:ncol(temp)-1])
  ranksvm<-sort(ranksvm,decreasing = TRUE)
  rm(temp)
  rm(svm_w)
 
  
  ###T-testresult
  ranktest=seq(length=2842, from=0, to=0)
  testlist<-del_test(train2)
  for(i in 1:10){
    ranktest=ranktest+testlist[[i]]
  }
  ranktest=abs(ranktest)/10
  names(ranktest) <- names
  ranktest<-sort(ranktest,decreasing = TRUE)

  
  ###pearson
  pearson1<-list()
  for(i in 1:10){
    temptrain=train[[i]]
    pearson1[[i]]<-cor(temptrain[1:ncol(temptrain)-1],
                       temptrain[ncol(temptrain)],method = "pearson")
    rm(temptrain)
  }
  rankpearson=seq(length=2842, from=0, to=0)
  pearson2<-del_list(pearson1)
  for(i in 1:10){
    rankpearson=rankpearson+pearson2[[i]]
  }
  rm(pearson1)
  rankpearson=abs(rankpearson)/10
  rankpearson<-as.vector(rankpearson)
  names(rankpearson) <- names
  rankpearson<-sort(rankpearson,decreasing = TRUE)
  
}

rankingsvm<-list()
rankingtest<-list()
rankingfc<-list()
rankingpearson<-list()

#train the svm model
for(j in 5:20){
  temp<-vector()
  for(i in 1:10){
    genesnames<-names(head(ranksvm,j))###
    genesnames<-union(genesnames,"class")
    temptrain<-train[[i]]
    temptest<-test[[i]]
    temptrain<-temptrain[,genesnames]
    temptest<-temptest[,genesnames]
    
    #get the results
    cat("\n-->",j,i)
    tempmodel<-test_svm3(temptrain,temptest)
    temp<-rbind(temp,tempmodel)
    }
  rankingsvm[[j-4]]=temp
}
 

#train the foldchange model
for(j in 5:20){
  temp<-vector()
  for(i in 1:10){
    genesnames<-names(head(rankfc,j))###
    genesnames =gsub("-", ".", genesnames, fixed = TRUE)
    genesnames<-union(genesnames,"class")
    temptrain<-train[[i]]
    temptest<-test[[i]]
    temptrain<-temptrain[,genesnames]
    temptest<-temptest[,genesnames]
    
    #get the results
    cat("\n-->",j,i)
    tempmodel<-test_svm4(temptrain,temptest)
    temp<-rbind(temp,tempmodel)
  }
  rankingfc[[j-4]]=temp
}



#train the pearson model
for(j in 5:20){
  temp<-vector()
  for(i in 1:10){
    genesnames<-names(head(rankpearson,j))###
    genesnames =gsub("-", ".", genesnames, fixed = TRUE)
    genesnames<-union(genesnames,"class")
    temptrain<-train[[i]]
    temptest<-test[[i]]
    temptrain<-temptrain[,genesnames]
    temptest<-temptest[,genesnames]
    
    #get the results
    cat("\n-->",j,i)
    tempmodel<-test_svm4(temptrain,temptest)
    temp<-rbind(temp,tempmodel)
  }
  rankingpearson[[j-4]]=temp
}

#train the t-test model
for(j in 5:20){
  temp<-vector()
  for(i in 1:10){
    genesnames<-names(head(ranktest,j))###
    genesnames =gsub("-", ".", genesnames, fixed = TRUE)
    genesnames<-union(genesnames,"class")
    temptrain<-train[[i]]
    temptest<-test[[i]]
    temptrain<-temptrain[,genesnames]
    temptest<-temptest[,genesnames]
    
    #get the results
    cat("\n-->",j,i)
    tempmodel<-test_svm4(temptrain,temptest)
    temp<-rbind(temp,tempmodel)
  }
  rankingtest[[j-4]]=temp
}







 for(k in 1:1){

  result_svmranking<-vector()
  result_testranking<-vector()
  result_fcranking<-vector()
  result_pearsonranking<-vector()
  
  for(j in 1:16){
  meanacc<-mean(rankingsvm[[j]][,1])
  temp<-union(meanacc,j+4)
  result_svmranking<-cbind(result_svmranking,temp)
  }
  for(j in 1:16){
  meanacc<-mean(rankingtest[[j]][6,])
  temp<-union(meanacc,j+4)
  result_testranking<-cbind(result_testranking,temp)
  }
  for(j in 1:16){
    meanacc<-mean(rankingfc[[j]][6,])
    temp<-union(meanacc,j+4)
    result_fcranking<-cbind(result_fcranking,temp)
  }
  for(j in 1:16){
    meanacc<-mean(rankingpearson[[j]][6,])
    temp<-union(meanacc,j+4)
    result_pearsonranking<-cbind(result_pearsonranking,temp)
  }


}

result2<-vector()
result2<-cbind(result2,result_svmranking[1,])
result2<-cbind(result2,result_fcranking[1,])
result2<-cbind(result2,result_testranking[1,])
result2<-cbind(result2,result_pearsonranking[1,])


test_svm3 = function(train,test){
  result2<-list()
  result<-vector()
  traindata=train
  testdata=test
  obj<-tune.svm(class~.,data=traindata,gamma = 10^0,cost = 10^(-3:3),kernel = "linear") #gamma=-3
  best_gm <- obj$best.parameters$gamma
  best_ct <- obj$best.parameters$cost
  model <- svm(class~.,data=traindata,
               kernel = "linear", gamma = best_gm, cost = best_ct, cross = 0,probability = TRUE)
  
  pre_label <- predict(model, testdata, probability = TRUE)
  
  stat_res <- table(pre_label, testdata$class)
  accurary <- (stat_res[1, 1] + stat_res[2, 2])/length(pre_label)
  sensitivity <- stat_res[2, 2]/(stat_res[1, 2] + stat_res[2, 2])
  specificity <- stat_res[1, 1]/(stat_res[1, 1] + stat_res[2, 1])
  result<-cbind(accurary,sensitivity,specificity)
  return(result)
}

test_svm4 = function(train,test){
  result2<-list()
  result<-vector()
  traindata=train
  testdata=test
  obj<-tune.svm(class~.,data=traindata,gamma = 10^0,cost = 10^(-1:1),kernel = "linear") #gamma=-3
  best_gm <- obj$best.parameters$gamma
  best_ct <- obj$best.parameters$cost
  model <- svm(class~.,data=traindata,
               kernel = "linear", gamma = best_gm, cost = best_ct, cross = 0,probability = TRUE)
  
  pre_label <- predict(model, testdata, probability = TRUE)
  
  stat_res <- table(pre_label, testdata$class)
  accurary <- (stat_res[1, 1] + stat_res[2, 2])/length(pre_label)
  sensitivity <- stat_res[2, 2]/(stat_res[1, 2] + stat_res[2, 2])
  specificity <- stat_res[1, 1]/(stat_res[1, 1] + stat_res[2, 1])
  result<-accurary
  return(result)
}




#Get the T statistic for a set of samples
del_test<-function(train){
  temp<-list()
  for(i in 1:10){
    temp2<-vector()
    a<-train[[i]]
    a1<-subset(a,class=="1")
    a2<-subset(a,class=="2")
    for(j in 1:2842){
      a3<-a1[,j]
      a4<-a2[,j]
      a5=t.test(a3,a4,paired = FALSE)
      temp2[j]<-a5$statistic
    }
  temp[[i]]=temp2
  }
  return(temp)
}



