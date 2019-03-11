library("GOplot")
library("DOSE")
library("clusterProfiler")
library("ggplot2")

gene_GO<-c("CALML3",'CERS3','DSG3','KRT16','KRT6A','KRT6B','LINC01206','SPRR1B','BUB1B','NCAPH','TMPRSS11A')
eg = bitr(gene_GO, fromType="SYMBOL", toType="ENTREZID",OrgDb="org.Hs.eg.db", drop = TRUE)

target_gene_id = as.character(eg[,2])
display_number = c(1, 13, 21)
ego_MF <- enrichGO(gene = target_gene_id,
                OrgDb = org.Hs.eg.db,
                ont = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[1], ]

ego_CC <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_id,
                   pvalueCutoff = 0.05,
                   ont = "CC",
                   readable=TRUE)
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]


ego_BP <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_id,
                   pvalueCutoff = 0.05,
                   ont = "BP",
                   readable=TRUE)
ego_result_BP <- na.omit(as.data.frame(ego_BP)[1:display_number[3], ])


category=factor(c(rep("biological process", display_number[1]), rep("cellular component", display_number[2]),
                  rep("molecular function", display_number[3])), 
                levels=c("molecular function", "cellular component", "biological process"))


#KEGG Enrichment
kk<- enrichKEGG(gene= target_gene_id,
                 organism = "hsa",keyType = "kegg", pvalueCutoff = 0.15,
                 pAdjustMethod = "BH",qvalueCutoff = 0.20)

display_number2 = c(8)
ego_result_kegg <- na.omit(as.data.frame(kk )[1:display_number2[1], ])

        