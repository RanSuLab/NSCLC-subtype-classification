library("GOplot")

ego_result=read.csv("ego_result.csv",sep = ',')
gene=read.csv("gene.csv",sep = ',')
circ<-circle_dat(ego_result,gene)
GOBar(circ,display='multiple')

chord<-chord_dat(data=circ,genes = gene)


pdf(file="plot_f.pdf",width=24,height=24)
GOChord(chord,space=0,gene.size=10)
dev.off()




