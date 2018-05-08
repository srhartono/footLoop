library(ggplot2)
library(grid)
library(gridExtra)
dm = read.table("RESULTS.tsv",header=T,sep="\t")
dm$restypes = paste(dm$restype1,dm$restype2,sep="_")
groups = unique(dm$group)
restypes = unique(dm$restypes)
for (i in 1:length(groups)) {
for (j in 1:length(restypes)) {

df = dm[dm$group == groups[i] & dm$restypes == restypes[j],]
if (length(df) > 0) {
print(paste(groups[i],restypes[j],dim(df)))
origA = df[df$shuftype == "ORIG" & df$peakname == "peak1",]
origB = df[df$shuftype == "ORIG" & df$peakname == "peak2",]
shuf1A = df[df$shuftype == "SHUF1" & df$peakname == "peak1",]
shuf1B = df[df$shuftype == "SHUF1" & df$peakname == "peak2",]
shuf2A = df[df$shuftype == "SHUF2" & df$peakname == "peak1",]
shuf2B = df[df$shuftype == "SHUF2" & df$peakname == "peak2",]
origA = origA[order(as.integer(origA$beg1/100),as.integer(origA$end1/100),as.integer((origA$end1-origA$beg1)/100)),];origA$y = seq(1,dim(origA)[1])
origB = origB[order(as.integer(origB$beg2/100),as.integer(origB$end2/100),as.integer((origB$end2-origB$beg2)/100)),];origB$y = seq(1,dim(origB)[1])
shuf1A = shuf1A[order(as.integer(shuf1A$beg1/100),as.integer(shuf1A$end1/100),as.integer((shuf1A$end1-shuf1A$beg1)/100)),];shuf1A$y = seq(1,dim(shuf1A)[1])
shuf1B = shuf1B[order(as.integer(shuf1B$beg2/100),as.integer(shuf1B$end2/100),as.integer((shuf1B$end2-shuf1B$beg2)/100)),];shuf1B$y = seq(1,dim(shuf1B)[1])
shuf2A = shuf2A[order(as.integer(shuf2A$beg1/100),as.integer(shuf2A$end1/100),as.integer((shuf2A$end1-shuf2A$beg1)/100)),];shuf2A$y = seq(1,dim(shuf2A)[1])
shuf2B = shuf2B[order(as.integer(shuf2B$beg2/100),as.integer(shuf2B$end2/100),as.integer((shuf2B$end2-shuf2B$beg2)/100)),];shuf2B$y = seq(1,dim(shuf2B)[1])
origA.cor = cor(origA$ypos,origA$xpos,method="pearson")
origB.cor = cor(origB$ypos,origB$xpos,method="pearson")
shuf1A.cor = cor(shuf1A$ypos,shuf1A$xpos,method="pearson")
shuf1B.cor = cor(shuf1B$ypos,shuf1B$xpos,method="pearson")
shuf2A.cor = cor(shuf2A$ypos,shuf2A$xpos,method="pearson")
shuf2B.cor = cor(shuf2B$ypos,shuf2B$xpos,method="pearson")

pdf(paste(groups[i],"_",restypes[j],".pdf",sep=""),width=42,height=7);
p1 = ggplot(origA,aes(xmin=beg1,xmax=end1,ymin=y,ymax=y+1)) + geom_rect(fill=rgb(1,1,1,0),color="black") + theme_bw() + theme(panel.grid=element_blank()) + annotate(geom="text",x=0,y=0,label="PCB1",hjust=0)
p2 = ggplot(origA,aes(xmin=beg2,xmax=end2,ymin=y,ymax=y+1)) + geom_rect(fill=rgb(1,1,1,0),color="black") + theme_bw() + theme(panel.grid=element_blank()) + annotate(geom="text",x=0,y=0,label="PCB5",hjust=0)
p3 = ggplot(origA,aes(x=xpos,y=ypos)) + geom_point(size=0.5,alpha=0.5) + annotate(geom="segment",x=0,y=0,xend=5000,yend=5000) + coord_cartesian(xlim=c(0,6000),ylim=c(0,6000)) + stat_smooth(method="lm") + annotate(geom="text",x=2000,y=2000,label=origA.cor) + theme_bw() + theme(panel.grid=element_blank())
p4 = ggplot(origB,aes(xmin=beg2,xmax=end2,ymin=y,ymax=y+1)) + geom_rect(fill=rgb(1,1,1,0),color="black") + theme_bw() + theme(panel.grid=element_blank()) + annotate(geom="text",x=0,y=0,label="PCB5",hjust=0)
p5 = ggplot(origB,aes(xmin=beg1,xmax=end1,ymin=y,ymax=y+1)) + geom_rect(fill=rgb(1,1,1,0),color="black") + theme_bw() + theme(panel.grid=element_blank()) + annotate(geom="text",x=0,y=0,label="PCB1",hjust=0)
p6 = ggplot(origB,aes(x=xpos,y=ypos)) + geom_point(size=0.5,alpha=0.5) + annotate(geom="segment",x=0,y=0,xend=5000,yend=5000) + coord_cartesian(xlim=c(0,6000),ylim=c(0,6000)) + stat_smooth(method="lm") + annotate(geom="text",x=2000,y=2000,label=origB.cor)+ theme_bw() + theme(panel.grid=element_blank())
grid.arrange(p1,p2,p3,p4,p5,p6,nrow=1,ncol=6)
p1 = ggplot(shuf1A,aes(xmin=beg1,xmax=end1,ymin=y,ymax=y+1)) + geom_rect(fill=rgb(1,1,1,0),color="black") + theme_bw() + theme(panel.grid=element_blank()) + annotate(geom="text",x=0,y=0,label="PCB1",hjust=0)
p2 = ggplot(shuf1A,aes(xmin=beg2,xmax=end2,ymin=y,ymax=y+1)) + geom_rect(fill=rgb(1,1,1,0),color="black") + theme_bw() + theme(panel.grid=element_blank()) + annotate(geom="text",x=0,y=0,label="PCB5 shuf",hjust=0)
p3 = ggplot(shuf1A,aes(x=xpos,y=ypos)) + geom_point(size=0.5,alpha=0.5) + annotate(geom="segment",x=0,y=0,xend=5000,yend=5000) + coord_cartesian(xlim=c(0,6000),ylim=c(0,6000)) + stat_smooth(method="lm") + annotate(geom="text",x=2000,y=2000,label=shuf1A.cor) + theme_bw() + theme(panel.grid=element_blank())
p4 = ggplot(shuf1B,aes(xmin=beg2,xmax=end2,ymin=y,ymax=y+1)) + geom_rect(fill=rgb(1,1,1,0),color="black") + theme_bw() + theme(panel.grid=element_blank()) + annotate(geom="text",x=0,y=0,label="PCB5 shuf",hjust=0)
p5 = ggplot(shuf1B,aes(xmin=beg1,xmax=end1,ymin=y,ymax=y+1)) + geom_rect(fill=rgb(1,1,1,0),color="black") + theme_bw() + theme(panel.grid=element_blank()) + annotate(geom="text",x=0,y=0,label="PCB1",hjust=0)
p6 = ggplot(shuf1B,aes(x=xpos,y=ypos)) + geom_point(size=0.5,alpha=0.5) + annotate(geom="segment",x=0,y=0,xend=5000,yend=5000) + coord_cartesian(xlim=c(0,6000),ylim=c(0,6000)) + stat_smooth(method="lm") + annotate(geom="text",x=2000,y=2000,label=shuf1B.cor) +theme_bw() + theme(panel.grid=element_blank())
grid.arrange(p1,p2,p3,p4,p5,p6,nrow=1,ncol=6)
p1 = ggplot(shuf2A,aes(xmin=beg1,xmax=end1,ymin=y,ymax=y+1)) + geom_rect(fill=rgb(1,1,1,0),color="black") + theme_bw() + theme(panel.grid=element_blank()) + annotate(geom="text",x=0,y=0,label="PCB1 shuf",hjust=0)
p2 = ggplot(shuf2A,aes(xmin=beg2,xmax=end2,ymin=y,ymax=y+1)) + geom_rect(fill=rgb(1,1,1,0),color="black") + theme_bw() + theme(panel.grid=element_blank()) + annotate(geom="text",x=0,y=0,label="PCB5",hjust=0)
p3 = ggplot(shuf2A,aes(x=xpos,y=ypos)) + geom_point(size=0.5,alpha=0.5) + annotate(geom="segment",x=0,y=0,xend=5000,yend=5000) + coord_cartesian(xlim=c(0,6000),ylim=c(0,6000)) + stat_smooth(method="lm") + annotate(geom="text",x=2000,y=2000,label=shuf2A.cor) +theme_bw() + theme(panel.grid=element_blank())
p4 = ggplot(shuf2B,aes(xmin=beg2,xmax=end2,ymin=y,ymax=y+1)) + geom_rect(fill=rgb(1,1,1,0),color="black") + theme_bw() + theme(panel.grid=element_blank()) + annotate(geom="text",x=0,y=0,label="PCB5",hjust=0)
p5 = ggplot(shuf2B,aes(xmin=beg1,xmax=end1,ymin=y,ymax=y+1)) + geom_rect(fill=rgb(1,1,1,0),color="black") + theme_bw() + theme(panel.grid=element_blank()) + annotate(geom="text",x=0,y=0,label="PCB1 shuf",hjust=0)
p6 = ggplot(shuf2B,aes(x=xpos,y=ypos)) + geom_point(size=0.5,alpha=0.5) + annotate(geom="segment",x=0,y=0,xend=5000,yend=5000) + coord_cartesian(xlim=c(0,6000),ylim=c(0,6000)) + stat_smooth(method="lm") + annotate(geom="text",x=2000,y=2000,label=shuf2B.cor) +theme_bw() + theme(panel.grid=element_blank())
grid.arrange(p1,p2,p3,p4,p5,p6,nrow=1,ncol=6)
dev.off()
}
}
}
#orig=c(0.9126,0.899)
#shuf=c(0.57,0.66,0.63,0.78)
#t.test(orig,shuf)
#boxplot(orig,shuf,ylim=c(0,1),xlab="original                           shuffled",ylab="Pearson Correlation",main="CALM3 PCB1 vs PCB5")
#text(1.5,0.8,"P value=0.01049")

