TOM = TOMsimilarityFromExpr(datExpr,power = sft$powerEstimate,corType=corType,networkType=type)
load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
dissTOM = 1-TOM
plotTOM = dissTOM^7
rm(dissTOM)
diag(plotTOM) = 0
TOMplot(plotTOM, net$dendrograms, mergedColors)
rm(plotTOM)
