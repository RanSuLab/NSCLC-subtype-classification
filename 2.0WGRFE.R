library("e1071")
library("randomForest")
library("dplyr")

set.seed(123)

svmresult<-list()
rfresult<-list()


#run the WGRFE-SVM
for(i in 1:10){
  cat(i)
  ranks=R2[[i]]
  svmresult[[i]]<-rfe2(train2[[i]],ranks)
}

#run the WGRFE-RF
for(i in 1:10){
  cat(i)
  ranks=R2[[i]]
  rfresult[[i]]<-rfe3(train2[[i]],ranks)
}
# get the gene selected by WGRFE-RF
acc_svmresult<-t(del_result(svmresult))
acc_svmresult=apply(acc_svmresult,2,as.numeric)
genesvm<-get_gene2(acc_svmresult)



# get the gene selected by WGRFE-RF
acc_rfresult<-t(del_result(rfresult))
acc_rfresult=apply(acc_rfresult,2,as.numeric)
generf<-get_gene(acc_rfresult)


geneWGRFE<-rbind(genesvm,generf)

geneWGRFE<-as.data.frame(table(genesvm))













rfe2 = function(traindata,ranks) 
{
  info<-vector()
  info2<-list()
  names<-names(ranks)
  ranks<-as.vector(ranks)
  names(ranks)<-names
  
  maxacc = Inf
  feat = colnames(traindata)
  j=1
  
  while(ncol(traindata) > 1){		
    
    cat("\n-->",ncol(traindata)-1, " features:")
    obj<-tune.svm(class~.,data=traindata,gamma = 10^(0:0),cost = 10^(-2:3),kernel = "linear") #gamma=-3
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
    stepmodel$s<- w*ranks[names(w)]*100000
    ord      = order(stepmodel$s)
    
    if (ncol(traindata)>20){
      remove   = colnames(traindata)[ord[1:round(ncol(traindata)*0.10)]]}else{
        remove   = colnames(traindata)[ord[1:1]]}
    
    
    traindata        = traindata[, setdiff(colnames(traindata), remove), drop=FALSE]
    
    
    cat("===> selected ", length(feat)-1, "features: ","\n\n\n")
    svmrfe<-info2}
  return(svmrfe)
}




rfe3 = function(traindata,ranks)
{
  info<-list()
  info2<-list()
  names<-names(ranks)
  names =gsub("-", ".", names, fixed = TRUE)
  ranks<-as.vector(ranks)
  names(ranks)<-names
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
    b=ranks[c(row.names(im))]
    a=as.vector(importance(stepmodel,type=2))
    im=b*a*100000
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

del_result<-function(rfresult){
  genenumber<-vector()
  temp<-rfresult[[1]]
  for(j in 1:66){
    temp2<-temp[[j]]
    temp3<-temp2[1]
    genenumber<-cbind(genenumber,temp3)
  }
  for( i in 1:10){
    temp<-rfresult[[i]]
    acc<-vector()
    for(j in 1:66){
      temp2<-temp[[j]]
      temp3<-as.numeric(temp2[1,2])
      acc<-cbind(acc,temp3)
     }
    genenumber<-rbind(genenumber,acc)
  }
  return(genenumber)
}


get_gene<-function(acc_rfresult){
  result<-vector()

for (i in 2:11){
  min1<-min(acc_rfresult[,2])
  for(j in 66:1){
    if((acc_rfresult[j,i]-min1)<0.005){
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

#used in SVM
get_gene2<-function(acc_rfresult){
  result<-vector()
  
  for (i in 2:11){
    min1<-min(acc_rfresult[,2])
    for(j in 66:1){
      if((acc_rfresult[j,i]-min1)<0.010){
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
