net = blockwiseModules(datExpr,power = sft$powerEstimate, maxBlockSize = 7000,
                       TOMType = "unsigned", deepSplit =2 ,minModuleSize =20,
                       reassignThreshold = 0, mergeCutHeight =0.25
                       ,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "FPKM-TOM",
                       verbose = 3)
table(net$colors)
# open a graphics window
#sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    groupLabels = c("Module colors", 
                                    "GS.weight"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)



