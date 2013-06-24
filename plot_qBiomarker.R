dCt<-read.table("/Users/stevensmith/bin/qpcr/qPcr_16sagreement",sep="\t",row.names=1,header=T)

#dCt<-data.frame(apply(dCt,2,function(x) ((x-min(x)))))

#dCt[dCt=="NP"]<-"NA"
abundance_16s<-read.table("/Users/stevensmith/bin/qpcr/16s_alltimes",sep="\t",row.names=1,header=T)

plot(dCt$W3D2,log(abundance_16s$W3D2),col=1,xlim=c(min(dCt,na.rm=TRUE),max(dCt,na.rm=TRUE)),xlab="dCt",ylab="log(abunance)",main="16S vs qPCR Methods",cex=2,lwd=15)
colors=c(1,2,3,4,5,6,8)
for (i in seq(from=2, to=7)){
  par(new=T)
  plot(dCt[,i],log(abundance_16s[,i]),axes=F,col=colors[i],xlab="",ylab="",cex=2, lwd=20)
}
legend(x=10,y=7,legend=colnames(dCt),fill=colors)

log(100)