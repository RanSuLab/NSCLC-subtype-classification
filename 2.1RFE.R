library("e1071")
library("randomForest")
library("dplyr")

set.seed(123)

svmresult2<-list()
rfresult2<-list()


#run the RFE-SVM
for(i in 1:10){
  cat(i)
  svmresult2[[i]]<-rfe4(train[[i]])
}


#run the RFE-RF
for(i in 1:10){
  cat(i)
  rfresult2[[i]]<-rfe5(train[[i]])
}


# get the gene selected by RFE-RF
acc_svmresult2<-t(del_result2(svmresult2))
acc_svmresult2=apply(acc_svmresult2,2,as.numeric)
genesvm2<-get_gene2(acc_svmresult2)

# get the gene selected by WGRFE-RF
acc_rfresult2<-t(del_result2(rfresult2))
acc_rfresult2=apply(acc_rfresult2,2,as.numeric)
generf2<-get_gene3(acc_rfresult2)


genesvm2<-as.data.frame(table(genesvm2))
generf2<-as.data.frame(table(generf2))



rfe4 = function(traindata)
{
  info<-vector()
  info2<-list()
  maxacc = Inf
  feat = colnames(traindata)
  j=1
  
  while(ncol(traindata) > 1){		
    
    cat("\n-->",ncol(traindata)-1, " features:")
    obj<-tune.svm(class~.,data=traindata,gamma = 10^(0:0),cost = 10^(-2:2),kernel = "linear") #gamma=-3
    best_gm <- obj$best.parameters$gamma
    best_ct <- obj$best.parameters$cost
    stepmodel <- svm(class~.,data=traindata,
                 kernel = "linear", gamma = best_gm, cost = best_ct, cross = 0,probability = TRUE)
    
    w<-t(find_SVM_weight(stepmodel))
    w<-as.vector(w)
    names(w)<-colnames(traindata[1:ncol(traindata)-1])
    if(obj$best.performance <= maxacc){
      maxacc = obj$best.performance
      feat  = colnames(traindata)			
      fit   = stepmodel
    }
    cat(obj$best.performance,"\n")
    info<-cbind(ncol(traindata)-1,obj$best.performance,obj$best.parameters$cost,colnames(traindata[1:ncol(traindata)-1]))
    
    info2[[j]]=info
    j=j+1
    
    
    stepmodel$s<- w
    
    
    ord      = order(stepmodel$s)
    
    if (ncol(traindata)>20){
    remove   = colnames(traindata)[ord[1:round(ncol(traindata)*0.10)]]}else{
    remove   = colnames(traindata)[ord[1:1]]}
    

    traindata        = traindata[, setdiff(colnames(traindata), remove), drop=FALSE]

  
  cat("===> selected ", length(feat)-1, "features: ","\n\n\n")
  svmrfe<-info2
  }
  return(svmrfe)
  
  }

find_SVM_weight<-function(model)
{
  W <- NULL
  
  if (model$nclasses==2){
    
    ### Weight for binary classification
    W <- t(model$coefs) %*% model$SV
    
  }else{
    ### Weight for multi-class classification
    start <- c(1, cumsum(model$nSV) + 1)
    start <- start[-length(start)]
    
    get_weight <- function(i,j){
      ri <- start[i] : (start[i] + model$nSV[i] - 1)
      rj <- start[j] : (start[j] + model$nSV[j] - 1)
      
      coef1 <- model$coefs[ri, j-1]
      coef2 <- model$coefs[rj, i]
      
      t(coef1) %*% model$SV[ri,] + t(coef2) %*% model$SV[rj,]
    }		
    
    for (i in 1:(model$nclasses-1)){
      for (j in (i+1):model$nclasses){
        W <- rbind(W, get_weight(i,j))
      }
    }		
  }
  return(W)
}

rfe5 = function(traindata)
{
  info<-list()
  info2<-list()
  names<-colnames(traindata)
  names =gsub("-", ".", names, fixed = TRUE)
  colnames(traindata)<-names
  maxacc = Inf
  feat = colnames(traindata)
  j=1
  while(ncol(traindata) > 1){		
    
    cat("\n-->",ncol(traindata)-1, " features:")
    stepmodel= randomForest(class~.,traindata,importance=TRUE)
    im=stepmodel$importance
    a=as.vector(importance(stepmodel,type=2))
    im=a
    names(im)<-colnames(traindata[1:ncol(traindata)-1])
    if(mean(stepmodel$err.rate) <= maxacc){
      maxacc = mean(stepmodel$err.rate)
      feat  = colnames(traindata)			
      fit   = stepmodel
    }
    info<-cbind(ncol(traindata)-1,mean(stepmodel$err.rate),colnames(traindata[1:ncol(traindata)-1]))
    
    info2[[j]]=info
    j=j+1
    cat(mean(stepmodel$err.rate),"\n")
    
    stepmodel$s<- im
    
    
    ord      = order(stepmodel$s)
    
    if (ncol(traindata)>20){
      remove   = colnames(traindata)[ord[1:round(ncol(traindata)*0.10)]]}else{
        remove   = colnames(traindata)[ord[1:1]]}
    
    
    traindata        = traindata[, setdiff(colnames(traindata), remove), drop=FALSE]
    
    
    cat("===> selected ", length(feat)-1, "features: ","\n\n\n")
    rfrfe<-info2
  }
  return(rfrfe)
}

del_result2<-function(rfresult){
  genenumber<-vector()
  temp<-rfresult[[1]]
  for(j in 1:72){
    temp2<-temp[[j]]
    temp3<-temp2[1]
    genenumber<-cbind(genenumber,temp3)
  }
  for( i in 1:10){
    temp<-rfresult[[i]]
    acc<-vector()
    for(j in 1:72){
      temp2<-temp[[j]]
      temp3<-as.numeric(temp2[1,2])
      acc<-cbind(acc,temp3)
    }
    genenumber<-rbind(genenumber,acc)
  }
  return(genenumber)
}

get_gene3<-function(acc_rfresult){
  result<-vector()
  
  for (i in 2:11){
    min1<-min(acc_rfresult[,2])
    for(j in 72:1){
      if((acc_rfresult[j,i]-min1)<0.007){
        result[i-1]=j
        #print(j)
        break
      }
    }
  }
  
  generf<-list()
  for(i in 1:10){
    temp1<-rfresult[[i]]
    k<-result[i]
    temp2<-temp1[[k]]
    generf[[i]]<-temp2[,3]
  }
  
  
  dataTemp2 <- do.call(cbind, 
                       lapply(lapply(generf, unlist), `length<-`, 
                              max(lengths(generf))))
  dataTemp2[is.na(dataTemp2)] = 0
  
  
  return(dataTemp2)
  
  
}


