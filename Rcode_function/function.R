
#Remove abnormal samples for WGCNA
get_abnormal_samples<-function(expro.upper){
  datExpr=as.data.frame(t(expro.upper));
  gsg = goodSamplesGenes(datExpr, verbose = 3);
  gsg$allOK
  sampleTree = hclust(dist(datExpr), method = "average")
  #plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
  clust = cutreeStatic(sampleTree, cutHeight =20000, minSize = 10)
  table(clust)
  keepSamples = (clust==1)
  datExpr = datExpr[keepSamples, ]
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  return(datExpr)
}


get_tomresult<-function(tempTOM){
  probes = colnames(datExpr)
  sum=0
  for(i in 1:dim(tempTOM)[1]){
    for(j in 1:dim(tempTOM)[2])
      sum=sum+tempTOM[i,j]
  }
  number=dim(tempTOM)[1]*dim(tempTOM)[2]-dim(tempTOM)[1]
  mean=sum/number
  tempTOM<-ifelse(tempTOM >= 0.1, 1, 0)
  dimnames(tempTOM) = list(probes, probes)
  return(tempTOM)
}











#Get the training set and test set
getdata<-function(expression){
  train<-list()
  test<-list()
  for(i in 1:10) {
    train[[i]]=expression[-folds[[i]],];
    test[[i]]=expression[folds[[i]],];
  }
  data<-list(train,test)
  return(data)
}


#Get and merge foldchange resul
getlog2fc<-function(train1){
  probes = colnames(datExpr)
  tumorAD.1=train1[which(train1$class==1),]
  tumorAD.2=train1[which(train1$class==2),]
  tumor_meanAD.1=colMeans(tumorAD.1[,(1:ncol(expression)-1)])
  tumor_meanAD.2=colMeans(tumorAD.2[,(1:ncol(expression)-1)])
  fc12<-tumor_meanAD.1/tumor_meanAD.2
  fc<-fc12[probes]
  log2fc<-log2(fc)
  return( log2fc)
}

mergelog2fc<-function(train){
  log2fc<-list()
  for(i in 1:10){
    log2fc[[i]]=getlog2fc(train[[i]])
  }
  return(log2fc)
}


## the GeneRank function
geneRank <- function(W,ex,d,max.degree=Inf){
  ## Use Sparse Matrix Packagen 
  ex = abs(ex)
  norm_ex = ex/max(ex)
  degrees = pmin(max.degree, pmax(1,colSums(W), na.rm=T))
  dimW = dim(W)[1]
  A=matrix(0, nrow = dimW, ncol = dimW)
  diag(A) = 1
  D1=matrix(0, nrow = dimW, ncol = dimW)
  diag(D1) = 1.0/degrees
  A = A - d*(t(W) %*% D1)
  b  = (1-d) * norm_ex
  
  r = as.numeric(solve(A,b))
  return(r)
}

mergerank<-function(log2fc){
  R<-list()
  for(i in 1:10){
    R[[i]]<-geneRank(tomresult,log2fc[[i]],d=0.5,max.degree=Inf)
  }
  
  return(R)
} 

tempranking<-function(vector){
  probes = colnames(datExpr)
  ranks<- 1/rank(1/vector)
  names(ranks) <- probes
  return(ranks)
}

mergeranking<-function(list){
  rankingresult<-list()
  for(i in 1:length(list)){
    rankingresult[[i]]=tempranking(list[[i]])
  }
  return(rankingresult)
}

#Select genes from turquoise and grey
selecte_genes<-function(v){
  module =v
  inModule1 = (mergedColors==module[1])
  inModule2 = (mergedColors==module[2])
  inModule=vector()
  for(i in 1:length(inModule1)){
    if(inModule1[i]==TRUE || inModule2[i]==TRUE)
      inModule[i]=TRUE
    else
      inModule[i]=FALSE
    
  }
  
  return(inModule)
}

#Dimension reduction of the training set
deltrain<-function(train){
  train2<-list()
  for(i in 1:10){
    temp<-train[[i]]
    train2[[i]]<-temp[,inModule]
    
  }
  return(train2)
}


#Dimension reduction for single-column data
del_list<-function(log2fc){
  train2<-list()
  for(i in 1:10){
    temp<-log2fc[[i]]
    train2[[i]]<-temp[inModule]
    
  }
  return(train2)
  
}

#
delR<-function(R){
  R2<-list()
  for(i in 1:10){
    temp<-R[[i]]
    R2[[i]]<-temp[inModule]
    
  }
  return(R2)
}



