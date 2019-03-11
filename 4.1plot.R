result2<-t(result2)
for(k in 1:1){
par(bty="l",cex.axis=1.0,cex.lab=1.0)
number<-c(5:20)
accuracy=result2[1,]
plot(number,accuracy,type="l",col="red",lwd=3.5,cex=1,cex.main=1.3, cex.lab=1.2, cex.axis=1.0, font=1.0,ylim=c(0.7,0.95))
lines(number,result2[2,],type="l",col="blue",lwd=3.5,cex=1,cex.main=1.3, cex.lab=1.2, cex.axis=2.0, font=1.2)
lines(number,result2[3,],type="l",col="cyan4",lwd=3.5,cex=1,cex.main=1.3, cex.lab=1.2, cex.axis=2.0, font=1.2)
lines(number,result2[4,],type="l",col="green4",lwd=3.5,cex=1,cex.main=1.3, cex.lab=1.2, cex.axis=2.0, font=1.2)




legend("bottomright",inset = 0.01,c("Generanking","Fold change","T-statistic","Pearson rank"),
       text.col =c("red","blue","cyan4","green4"),title=,title.col="black",bty = "n",cex=1.2,horiz=FALSE)

}
