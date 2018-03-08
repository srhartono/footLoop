.libPaths( c(.libPaths(), "/home/mitochi/R/x86_64-pc-linux-gnu-library/3.2/", "/home/mitochi/R/x86_64-pc-linux-gnu-library/3.4/") )
library(labeling)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(RColorBrewer)
				df = read.table("/group/stella/Work/Data/Fastq/110101_pacbio/0_Samples/PCB5/0_Fastq/test//.CALL/CALM3_Pos_20_0.65_CH.PEAK.out",skip=1,sep="\t")
				colnames(df) = c("V1",seq(1,dim(df)[2]-1))
				if (dim(df)[1] < 1000) {
					h = hclust(dist(df[,-1]))
					df = df[h$order,]
				} else {
					mysum = apply(df[,-1],1,sum)
					df = df[order(mysum),]
				}
				df$id = gsub("^.+/([0-9]+)/ccs", "\\1", df$V1, perl=T)
			
					clust = read.table("/group/stella/Work/Data/Fastq/110101_pacbio/0_Samples/PCB5/0_Fastq/test//FOOTCLUST/.TEMP/CALM3_Pos_20_0.65_CH.PEAK.local.bed.clust",header=T,sep="\t")
					clust = clust[grep("^[0-9]+\\.[0-9]+$",clust$id,perl=T,invert=T),]
					clust$y = seq(1,dim(clust)[1])
					clust = subset(clust,select=c("id","y","clust"))
					clust2 = as.data.frame(aggregate(clust$y,by=list(clust$clust),min))
					clust2$max = aggregate(clust$y,by=list(clust$clust),max)$x
					colnames(clust2) = c("clust","ymin","ymax")
					clust2$xmin = 1
					clust2$xmax = 30
					clust2$clust = clust2$clust + 10
					df3 = merge(df,clust,by="id")
					df3 = subset(df3,select=c(-y,-id,-V1))
					df3clust = df3$clust
df4 = as.data.frame(matrix(nrow=max(df3$clust),ncol=dim(df3)[2]-1))
for (x in 1:max(df3$clust)) {
	for (y in 1:(dim(df3)[2]-1)) {
		a = df3[df3$clust == x,y]
		mymax = max(a)
		peak = length(a[a == 8 | a == 9])
		conv = length(a[a == 6 | a == 7])
		none = length(a[a == 4 | a == 5])
		if (mymax == 1) {
			df4[x,y] = 99
		} else if (peak != 0 & peak/2 >= conv) {#mymax == 8 | mymax == 9) {
			df4[x,y] = as.integer(length(a[a == 8 | a == 9]) / length(a) * 9+0.5)
		} else if (conv != 0 & conv > peak/2) {#mymax == 6 | mymax == 7) {
			df4[x,y] = as.integer(length(a[a == 6 | a == 7]) / length(a) * -9+0.5)
		} else if (mymax == 4 | mymax == 5) {
			df4[x,y] = 0
		}	else if (mymax == 0) {
			df4[x,y] = -99
		} else {
			df4[x,y] = 100
		}
	}
}

					colnames(df4) = seq(1,dim(df4)[2])
					df4$clust = seq(1,max(df3clust))
					df4$y = df4$clust
					df3 = melt(df4,id.vars=c("y","clust"))
					#df3 = melt(df3,id.vars=c("y","clust")); 
					colnames(df3) = c("y","clust","x","value")
					df3$x = as.numeric(as.character(df3$x))
					df3$clust = as.numeric(as.character(df3$clust))
					df3$y = as.numeric(as.character(df3$y))
					df3$value = as.numeric(as.character(df3$value))
					df5 = df3; df5[df5$value > 10 | df5$value < -10,]$value = 0
					df4 = data.frame(clust=seq(1,max(df5$clust)),xmin=-1,xmax=-1,ymin=-1,ymax=-1)
					for (i in (min(df5$clust) : max(df5$clust))) {
						if (length(df5[df5$clust == i,]$y) > 0) {
							df4$xmin[i] = min(df5[df5$value > 0 & df5$clust == i,]$x)-0.5
							df4$xmax[i] = max(df5[df5$value > 0 & df5$clust == i,]$x)+0.5
							df4$ymin[i] = min(df5[df5$value > 0 & df5$clust == i,]$y)-0.5
							df4$ymax[i] = max(df5[df5$value > 0 & df5$clust == i,]$y)+0.5
						} else {
							df4 = df4[-i,]
						}
					}
					df5 = data.frame(clust=seq(1,max(df5$clust)),xmin=1,xmax=max(df5$x),
							ymin=seq(1,max(df5$clust))-0.5,ymax=seq(1,max(df5$clust))+0.5)
					df3$value = as.factor(df3$value)
					greens = rev(brewer.pal(9,"Greens"))
					reds = brewer.pal(9,"Reds")
					df3col=c(
"-99" = "grey",
"-9" = greens[1],
"-8" = greens[2],
"-7" = greens[3],
"-6" = greens[4],
"-5" = greens[5],
"-4" = greens[6],
"-3" = greens[7],
"-2" = greens[8],
"-1" = greens[9],
"0" = "cornsilk",
"1" = reds[1],
"2" = reds[2],
"3" = reds[3],
"4" = reds[4],
"5" = reds[5],
"6" = reds[6],
"7" = reds[7],
"8" = reds[8],
"9" = reds[9],
"99" = "white")
					clust = subset(clust,select=c("id","y"))
					df = merge(df,clust,by="id")
					df = subset(df,select=-id)
				bed = read.table("/group/stella/Work/Data/Fastq/110101_pacbio/0_Samples/PCB5/0_Fastq/test//PEAKS_LOCAL/CALM3_Pos_20_0.65_CH.PEAK.local.bed",sep="\t");
bed = merge(subset(df,select=c("V1","y")),bed,by="V1")

				dm = melt(df,id.vars=c("V1","y"))
				dm$variable = as.numeric(as.character(dm$variable))
				df3$x = as.numeric(as.character(df3$x))
				df3$y = as.numeric(as.character(df3$y))
				df3$value  = as.numeric(as.character(df3$value))
				df4$clust2 = df4$clust + 10

				p3 = ggplot(df3,aes(x,y)) + 
					geom_tile(aes(fill=as.factor(value))) +#,alpha=sqrt(abs(as.numeric(as.character(value))/9)))) +
					geom_rect(data=df5,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,x=xmin,y=ymin),
											size=0.8,fill=rgb(0,0,0,alpha=0),color=rgb(0,0,0,alpha=0.25)) +
					geom_rect(data=df4,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,x=xmin,y=ymin),
											size=2,fill=rgb(0,0,0,alpha=0),color=rgb(1,0,0,alpha=1)) +
					geom_rect(data=df4,aes(fill=as.factor(df4$clust2),x=1,y=1,xmin=1,xmax=30,ymin=ymin,ymax=ymax)) +
					geom_text(data=df4,aes(x=10,y=(ymin+ymax)/2,label=clust2-10),hjust=0,size=15) +
#					annotate(geom="segment",x=0,xend=2800,y=seq(0,5),yend=seq(0,5)) +
					theme_bw() + theme(legend.position="none") + 
					scale_fill_manual(values=c(df3col,"11"="#e41a1c","12"="#377eb8",
				"13"="#4daf4a","14"="#984ea3","15"="#ff7f00")) +
					scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
              theme(line = element_blank(),
                     axis.text = element_blank(),
                     axis.title = element_blank()
					) + coord_cartesian(ylim=c(-1,6))



				p = ggplot(dm,aes(variable,y)) +  
					geom_tile(aes(fill=as.factor(value))) + 
#					geom_tile(data=df3,aes(x=x,y=y,fill=as.factor(value),alpha=value*2)) + 
					geom_rect(data=clust2,aes(fill=as.factor(clust),x=xmin,y=ymin,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)) +
					geom_text(data=clust2,aes(group=as.factor(clust),x=10,y=(ymin+ymax)/2,label=clust-10),hjust=0,size=15) +
					theme_bw() + theme(legend.position="none") + 
#					coord_fixed() +
					scale_fill_manual(values=c("0"="grey","1"="white","4"="cornsilk","5"="cornsilk",
								  					   "6"="green4","7"="seagreen4","8"="red4","9"="maroon4",
														"11"="#e41a1c","12"="#377eb8","13"="#4daf4a","14"="#984ea3","15"="#ff7f00")) +
					scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
               theme(line = element_blank(),
                     axis.text = element_blank(),
                     axis.title = element_blank()
                    ) + ggtitle(paste("(peak=",103,"; nopk=",563,")",sep=""))
				
				if (length(bed) != 0 & dim(bed)[1] > 0) {
					bed$variable=1
					p = p + geom_rect(data=bed,aes(xmin=V2,xmax=V3,ymin=y-0.5,ymax=y+0.5),
											size=0.5,fill=rgb(0,0,0,alpha=0),color=rgb(1,0,0,alpha=0.25))
				}
				





            df2 = subset(df,select=c(-V1,-y));
				df2[df2 < 8] = 0; df2[df2 >= 8] = 1
            df2 = data.frame(x=seq(1,dim(df2)[2]), y=apply(df2,2,mean));
            if (dim(df2[df2$y > 0,])[1] > 15) {
               df2 = df2[df2$y > 0,]
               df2$x = as.numeric(as.character(df2$x));
               df2$y = as.numeric(as.character(df2$y));
               df2$x2 = df2$x
               df2$y2 = df2$y
               for (i in 1:(dim(df2)[1]-10)) {
                  a = df2[df2$x >= df2[i,]$x & df2$x <= df2[i+10-1,]$x,]
                  if (length(a) != 0 & dim(a)[1] != 0) {
                        df2[i,]$y2 = mean(a$y)
                        df2[i,]$x2 = mean(a$x)
                  }
               }
               mins = seq(1,as.integer(df2[1,]$x2)-1,10)
               maxs = seq(max(df2$x2),dim(df)[2]-2,10)
               df2 = rbind(data.frame(x=mins,y=0,x2=mins,y2=0),df2)
               df2 = rbind(df2,data.frame(x=maxs,y=0,x2=maxs,y2=0))
            } else {
               df2 = data.frame(x=seq(1,dim(df)[2]), y=0, x2=seq(1,dim(df)[2]), y2=0);
            }
            p2 = ggplot(df2,aes(x2,y2)) + geom_point(aes(x=x,y=y),size=1) + geom_line(color=rgb(1,0,0,alpha=1)) + theme_bw()+
               scale_x_continuous(expand = c(0,0)) +
               scale_y_continuous(expand = c(0,0)) +
               theme(line = element_blank(),axis.text = element_blank(),axis.title = element_blank()) +
               annotate(geom='text',x=10,y=1,label="- 100 %",size=5,hjust=0) +
               annotate(geom='text',x=10,y=0,label="- 0   %",size=5,hjust=0) +
               coord_cartesian(ylim=c(-0.05,1.05))
				ratio1 = as.integer(10*dim(df)[1]*16 / (dim(df)[1]*16 + 500 + 425)+0.5)/10
				ratio2 = as.integer(10*500           / (dim(df)[1]*16 + 500 + 425)+0.5)/10
				ratio3 = as.integer(10*425           / (dim(df)[1]*16 + 500 + 425)+0.5)/10

				png("/group/stella/Work/Data/Fastq/110101_pacbio/0_Samples/PCB5/0_Fastq/test//PNG/test_CALM3_Pos_20_0.65_CH.PEAK.out.png",width=dim(df)[2]*2,height=dim(df)[1]*16 + 500 + 425)
            grid.arrange(p,p2,p3,ncol=1,nrow=3,heights=c(ratio1,ratio2,ratio3));
				dev.off()

#				pdf("/group/stella/Work/Data/Fastq/110101_pacbio/0_Samples/PCB5/0_Fastq/test//PDF/test_CALM3_Pos_20_0.65_CH.PEAK.out.pdf",width=(dim(df)[2]*2),height=(dim(df)[1]*16 + 500 + 425 ))
#            grid.arrange(p,p2,ncol=1,nrow=2,heights=c(ratio1,ratio2));
#				dev.off()
			
