library("caret")
library("WGCNA")
library("stringr")


source("Rcode_function\\function.R")

#Initial setting
set.seed(123)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
corType = "pearson"
type = "unsigned"
robustY = ifelse(corType=="pearson",T,F)


#Read data
trait=read.csv('Clinical.txt',sep = '\t',row.names = 1)
expro.upper=read.csv('merge_ADSC.txt',sep = '\t',row.names = 1)
#Data preprocessing
datExpr<-get_abnormal_samples(expro.upper)
rm(expro.upper)
## WGCNA CODE
source("Rcode_WGCNA\\pickSoftThreshold.R")
source("Rcode_WGCNA\\cluster.R")
source("Rcode_WGCNA\\Eigengene cluster.R")
source("Rcode_WGCNA\\relation.R")
source("Rcode_WGCNA\\TOM.R")


#####GeneRank code  
###############
#transfor TOM into 0-1 matrix
diag(TOM) = 0
G<-get_tomresult(TOM)
rm(TOM)


#get the logfc of each gene
class<-read.csv('class.txt',sep = '\t')
expression<-cbind(datExpr,class)
expression$class<-factor(expression$class)
folds<-createFolds(y=expression$class,k=10)
rm(class)
rm(datExpr)
for(i in 1:1){
  col2<-colnames(expression)
  col2 =gsub("-", ".", col2, fixed = TRUE)
  colnames(expression)<-col2
  rm(col2)
}
data<-getdata(expression)
train<-data[[1]]
test<-data[[2]]
log2fc<-mergelog2fc(train)
R<-mergerank(log2fc)
#get the WGRFE ranking score in each fold 
rankingresult<-mergeranking(R)
rm(R)
rm(G)
#get the 2 moudule genes
module = c("turquoise","grey")
inModule<-selecte_genes(module)
train2<-deltrain(train)
R2<-delR(rankingresult)
rm(rankingresult)

##save R2,train2,train,log2fc







