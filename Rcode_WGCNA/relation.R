trait=read.csv('Clinical.txt',sep = '\t',row.names = 1)
# Read in phenotype data
traitData <- read.table(file="Clinical.txt", sep='\t', header=T, row.names=1, 
                          check.names=FALSE, comment='',quote="") 
sampleName = rownames(datExpr) 
traitData = traitData[match(sampleName, rownames(traitData)), ] 

### Module association with phenotype data
if (corType=="pearsoon") { 
  modTraitCor = cor(MEs_col, traitData, use = "p") 
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
}else { 
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY) 
  modTraitCor = modTraitCorP$bicor 
  modTraitP = modTraitCorP$p}



## Pearson correlation was used for individual columns with zero (or missing) 
## MAD. 
textMatrix = paste(signif(modTraitCor, 2)) 
dim(textMatrix) = dim(modTraitCor)
par(mar=c(2,8,2,2))
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), cex.lab = 0.6,ySymbols = colnames(MEs_col), 
               colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, 
              setStdMargins = FALSE, cex.text =0.6, zlim = c(-1,1))
 

