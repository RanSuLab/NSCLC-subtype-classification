# module eigengene
MEs = net$MEs 

 
MEs_col = MEs 
colnames(MEs_col) = paste0("ME", labels2colors( 
  as.numeric(str_replace_all(colnames(MEs),"ME","")))) 
MEs_col = orderMEs(MEs_col) 

plotEigengeneNetworks(MEs_col, NULL, marDendro = c(2,4,4,4), marHeatmap = c(3,3,2,3), 
                      plotDendrograms = T, xLabelsAngle = 90)


             