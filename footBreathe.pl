#!/usr/bin/perl

use strict; use warnings; use mitochy; use Getopt::Std; use FAlite;  use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_S $opt_G $opt_n $opt_R);
getopts("vSG:n:R");

BEGIN {

   my $softwarePath = dirname(dirname abs_path $0) . '/footLoop/softwares/';
   $ENV{PATH} = "$softwarePath/Bismark_v0.20.0/:$softwarePath/bedtools2/bin/:$softwarePath/bowtie2-2.2.6/:
$softwarePath/samtools-0.1.19/:$softwarePath/R-3.6.1/bin/:$ENV{PATH}";

}
my ($footPeakFolder) = ($opt_n);
my ($geneWant) = $opt_G;

die "
Usage: $YW $0$N -n $LCY<footPeakFolder>$N

-G: Genewant
-R: Random
-S: not implemented


" unless defined $opt_n and -d $opt_n;

my @peakFiles = <$footPeakFolder/.CALL/*.PEAK>;
die "\n\nNo .PEAK file found at $footPeakFolder/.CALL/\n\n" if @peakFiles == 0;

my $INFO = defined $opt_R ? "SHUF" : "ORIG";
my $window = 20; #fred's order 190117
my $coor = parse_footPeak_logFile($footPeakFolder);
my $outDir = defined $opt_R ? "$footPeakFolder/BREATHESHUF/" : "$footPeakFolder/BREATHE/";
my $outDirPIC = "$footPeakFolder/BREATHE/";
system("mkdir -p $outDirPIC") == 0 or die "Failed to mkdir -p $outDirPIC: $!\n";

foreach my $peakFileNTEMP (sort @peakFiles) {
	next if $peakFileNTEMP =~ /_Unk_/;

	#get gene name
	my ($gene) = $peakFileNTEMP =~ /_gene(.+)_(Pos|Neg|Unk)/; 
	my $strand = $coor->{uc($gene)}{strand};
	next if $strand eq "+" and $peakFileNTEMP !~ /Pos.+CH/;
	next if $strand eq "-" and $peakFileNTEMP !~ /Neg.+GH/;
	die "\n\nCannot prase gene from peakFile name\n$peakFileNTEMP\n\n" if not defined $gene;
	print "Processing gene=$gene\n";

	#next if gene isn't genewant
	if (defined $geneWant) {
		next if uc($geneWant) ne uc($gene);
	}

	#get nopeakfiles and revs
	my $nopkFileNTEMP = $peakFileNTEMP; $nopkFileNTEMP =~ s/\.PEAK$/\.NOPK/;
	my $peakFileTEMPL = $peakFileNTEMP;
	$peakFileTEMPL =~ s/gene(.+)_(Pos|Neg|Unk)_(.+)_(..)\.PEAK$/gene$1\_Neg\_$3\_GH.PEAK/ if $strand eq "+";
	$peakFileTEMPL =~ s/gene(.+)_(Pos|Neg|Unk)_(.+)_(..)\.PEAK$/gene$1\_Pos\_$3\_CH.PEAK/ if $strand eq "-";
	my $nopkFileTEMPL = $peakFileTEMPL; $nopkFileTEMPL =~ s/\.PEAK$/\.NOPK/;
	my ($folder1, $fileName1) = mitochy::getFilename($peakFileNTEMP, "folderfull");
	my ($folder2, $fileName2) = mitochy::getFilename($nopkFileNTEMP, "folderfull");
	
	my $data;

	my ($pos, $maxclust) = parse_clustFile($peakFileNTEMP, $coor, -1, "INFO");

	print "NEXTED AS TOTAL CLUSTER IS 0\n" and next if not defined $pos;

	# Get Beg and End of Cluster
	my ($begz, $endz, $clustz);
	foreach my $clust (sort keys %{$pos}) {
		$begz .= (not defined $begz) ?  $pos->{$clust}{beg} : ",$pos->{$clust}{beg}";
		$endz .= (not defined $endz) ? $pos->{$clust}{end} : ",$pos->{$clust}{end}";
		$clustz .= (not defined $clustz) ? $clust : ",$clust";
	}
	print "clustz=$clustz\nbegz=$begz\nendz=$endz\n";
	$begz = "c($begz)";
	$endz = "c($endz)";
	$clustz = "c($clustz)";
	my $bg;
	my $conv;
	open (my $outClust, ">", "$outDirPIC/$fileName1.$INFO.ClustConv") or die "Failed to write to $outDirPIC/$fileName1.$INFO.ClustConv: $!\n";
	($data->{FOR}) = process_PEAKNOPKFILE($peakFileNTEMP, $pos, $data->{FOR},$outClust, "INFO",0);
	($data->{REV}) = process_PEAKNOPKFILE($nopkFileTEMPL, $pos, $data->{REV},$outClust, "ORIG",0);
	if (defined $opt_R) {
		for (my $rep = 0; $rep < 25; $rep++) {
			my $datashuf;
			print "$rep\n";
			my ($pos2, $maxclust2) = parse_clustFile($peakFileNTEMP, $coor, $rep, "SHUF");
			process_PEAKNOPKFILE($nopkFileTEMPL, $pos2, $datashuf,$outClust, $INFO, $rep);
		}
	}
	close $outClust;
#	my %conv = %{$conv};
	my $undef = 0; 
	my %print;
	my $length = 0; my $lastlength = "INIT";
	foreach my $strandz (sort keys %{$data}) {
		print "STRAND = $strandz\n";
		foreach my $read (sort {$data->{$strandz}{$b}{sort}{1} <=> $data->{$strandz}{$a}{sort}{1}} keys %{$data->{$strandz}}) {
			my $clust = $data->{$strandz}{$read}{sort}{1};
			die "UNDEF CLUST $read\n" if not defined $clust;
#			$clust = $maxclust if $clust eq 0;
#			print $out1 "$read\t$clust\t$data->{$strandz}{$read}{value}\n";
			$length = scalar(split("\t", $data->{$strandz}{$read}{value}));
			die "Length = $length\n" if $lastlength ne "INIT" and $length ne $lastlength;
			my @val = split("\t", $data->{$strandz}{$read}{value});
			my $sum = 0;
			foreach my $val (@val) {
				$sum += $val >= 6 ? 1 : 0;
			}
			
			$print{$clust}{$strandz}{read}{$read}{print} = "$data->{$strandz}{$read}{value}" if $strandz eq "FOR";
			$print{$clust}{$strandz}{read}{$read}{print} = "$data->{$strandz}{$read}{value}" if $strandz eq "REV";
			$print{$clust}{$strandz}{read}{$read}{sum} = $sum;
			$lastlength = $length;
		}
	}
	foreach my $clust (sort {$a <=> $b} keys %print) {
		foreach my $strandz (sort keys %{$print{$clust}}) {
			foreach my $read (sort {$print{$clust}{$strandz}{read}{$b}{sum} <=> $print{$clust}{$strandz}{read}{$a}{sum}} keys %{$print{$clust}{$strandz}{read}}) {
				push(@{$print{$clust}{$strandz}{data}}, $print{$clust}{$strandz}{read}{$read}{print});
				push(@{$print{$clust}{$strandz}{name}}, $read);
			}
		}
	}
	my ($greycount, $blackcount) = (0,0);
	my $midlen = 3;
	open (my $out1, ">", "$outDir/$fileName1.out") or die "Cannot write to $outDir/$fileName1.out: $!\n";
	open (my $out2, ">", "$outDir/$fileName1.ALL") or die "Cannot write to $outDir/$fileName1.ALL: $!\n";
	my ($print1, $print2, $print3) = ("","","");
	foreach my $clust (sort {$a <=> $b} keys %print) {
		last if $clust eq 0 or $clust eq $maxclust;
		my @FOR   = defined $print{$clust}{FOR}{data} ? @{$print{$clust}{FOR}{data}} : ();
		my @REV   = defined $print{$clust}{REV}{data} ? @{$print{$clust}{REV}{data}} : ();
		my $totalrev = @REV;
		my $totalfor = @FOR;
		my @MID   = (-1) x $midlen;
		my $grey  = join("\t", ((-2)x$length));
		my $black = join("\t", ((-1)x$length));
		if ($totalfor < $totalrev) {
#			print "FOR < REV\n";
			my $add = 0;
			for (my $i = $totalfor; $i < $totalrev; $i++) {
				if ($totalfor == 0) {
					$FOR[$i] = $black;
					$print{$clust}{FOR}{name}[$i] = "black$greycount";
					$greycount ++;
					next;
				}
				last if $i + $add >= $totalrev;
				if ($i == $totalfor) {
					if ($totalfor < $totalrev - 3) {
						$FOR[$i+2] = $black; $add = 2;
						$print{$clust}{FOR}{name}[$i+2] = "grey$greycount";
#						print "i=" . ($i+2) . ", print=" . $print{$clust}{FOR}{name}[$i+2] . "\n";
						$greycount ++;
					}
					if ($totalfor < $totalrev - 2) {
						$FOR[$i+1] = $black; $add = 1 if $add eq 0;
						$print{$clust}{FOR}{name}[$i+1] = "grey$greycount";
						$greycount ++;
#						print "i=" . ($i+1) . ", print=" . $print{$clust}{FOR}{name}[$i+1] . "\n";
					}
					if ($totalfor < $totalrev - 1) {
						$FOR[$i+0] = $black; $add = 0 if $add eq 0;
						$print{$clust}{FOR}{name}[$i+0] = "grey$greycount";
						$greycount ++;
#						print "i=" . ($i+0) . ", print=" . $print{$clust}{FOR}{name}[$i+0] . "\n";
					}
#					print "1: i = $i + $add\n";
				}
				else {
					my $random = int(rand($totalfor));
					#$FOR[$i+$add] = $grey;
					$FOR[$i+$add] = $FOR[$random];
#					print "2: i = $i + $add\n";
					$print{$clust}{FOR}{name}[$i+$add] = $print{$clust}{FOR}{name}[$random];
				}
			}
		}
		else {
#			print "FOR > REV\n";
			my $add = 0;
			for (my $i = $totalrev; $i < $totalfor; $i++) {
				if ($totalrev == 0) {
					$REV[$i] = $black;
					$print{$clust}{REV}{name}[$i] = "black$greycount";
					$greycount ++;
					next;
				}
				last if $i + $add >= $totalfor;
				if ($i == $totalrev) {
					if ($totalrev < $totalfor - 3) {
						$REV[$i+2] = $black; $add = 2;
						$print{$clust}{REV}{name}[$i+2] = "grey$greycount";
						$greycount ++;
					}
					if ($totalrev < $totalfor - 2) {
						$REV[$i+1] = $black; $add = 1 if $add eq 0;
						$print{$clust}{REV}{name}[$i+1] = "grey$greycount";
						$greycount ++;
					}
					if ($totalrev < $totalfor - 1) {
						$REV[$i+0] = $black; $add = 0 if $add eq 0;
						$print{$clust}{REV}{name}[$i+0] = "grey$greycount";
						$greycount ++;
					}
				}
				else {
					my $random = int(rand($totalrev));
					#$REV[$i+$add] = $grey;
					$REV[$i+$add] = $REV[$random];
					$print{$clust}{REV}{name}[$i+$add] = $print{$clust}{REV}{name}[$random];
				}
			}
		}
		print "Totalfor = $totalfor, currentfor = " . scalar(@FOR) . "\ntotalrev = $totalrev, currentrev=" . scalar(@REV) . "\n\n";
		for (my $i = 0; $i < @FOR; $i++) {
			my $readName = $print{$clust}{FOR}{name}[$i];
			die "UNDEF READ NAME i=$i FOR\n" if not defined $readName;
			my @ALL = ($readName, $clust);
			my @CURRFOR = split("\t", $FOR[$i]);
			my @CURRREV = split("\t", $REV[$i]);
			for (my $j = 0; $j < @CURRFOR; $j++) {
				$CURRREV[$j] *= -1 if $CURRREV[$j] =~ /^[6789]$/;
				$ALL[$j+2] = $CURRFOR[$j] =~ /^[6789]$/ ? $CURRFOR[$j] : $CURRREV[$j] =~ /^\-?[6789]$/ ? $CURRREV[$j] : 
$CURRFOR[$j] =~ /^[45]$/ ? $CURRFOR[$j] : $CURRREV[$j] =~ /^[45]$/ ? $CURRREV[$j] : $CURRFOR[$j] ne 0 ? $CURRFOR[$j] :
$CURRREV[$j] ne 0 ? $CURRREV[$j] : $CURRFOR[$j] eq 1 ? $CURRFOR[$j] : $CURRREV[$j] eq 1 ? $CURRREV[$j] : $CURRFOR[$j];
			}
			print $out1 "$readName\t$clust\t" . $FOR[$i] . "\t" . join("\t", @MID) . "\t" . $REV[$i] . "\n";
			print $out2 join("\t", @ALL) . "\n";
		}
		print $out1 "clust$clust\t$clust\t" . $grey . "\t" . join("\t", @MID) . "\t" . $grey . "\n";
	}
	print $out1 "$print1";
	close $out1;
	
	open (my $outR, ">", "$outDir/$fileName1.R");
	print $outR "

	library(ggplot2)
	library(reshape2)
	library(GMD)
	library(Cairo)
	source(\"/home/mitochi/bin/theme_blank.R\")
	begz = $begz
	endz = $endz
	clustz = $clustz
	
	df = read.table(\"$outDir/$fileName1.out\");
	df = df[,-1]
	#,row.names=1)
	df\$clust = df\$V2
	df = df[,-1]
#	temp = subset(df,select=c(-clust))
#	temp[temp < 6 ] = 0
#	temp[temp >= 6] = 1
#	df\$mysum = apply(temp,1,sum)
	
	mycolorz=c(\"red4\",\"blue4\",\"green4\",\"purple\",\"black\",\"orange4\",\"cornflowerblue\",\"brown\",\"blue2\",\"red2\")
	mycolorz=c(mycolorz,\"red4\",\"blue4\",\"green4\",\"purple\",\"black\",\"orange4\",\"cornflowerblue\",\"brown\",\"blue2\",\"red2\")
	mycolorz=c(mycolorz,\"red4\",\"blue4\",\"green4\",\"purple\",\"black\",\"orange4\",\"cornflowerblue\",\"brown\",\"blue2\",\"red2\")
	mycolorz1 = mycolorz
	myclust = mycolorz[rev(df\$clust)]
	df = df[order(df\$clust),]#,-df\$mysum),]
	begz2 = c(0)
	endz2 = c(0)
	print(plyr::count(df\$clust))

#	ALL = read.table(\"$outDir/$fileName1.ALL\");
#	ALL = ALL[,-1]
#	ALL\$clust = ALL\$V2
#	ALL = ALL[,-1]
	
#	ALL = ALL[order(ALL\$clust),]#,-ALL\$mysum),]
#	begz2 = c(0)
#	endz2 = c(0)
#	print(plyr::count(ALL\$clust))

	dm = subset(df,select=c(-clust))#,-mysum))
#	dm = df
	#dm = dm[,apply(dm,2,median) > 1]
	mywidth = dim(dm)[2]-2
	myheight = dim(dm)[1]
	if (mywidth < 50) {mywidth = 50}
	if (myheight < 50) {myheight = 50}
	mycolorz = c(
		\"green3\",
		\"blue4\",
		\"green4\",
		\"grey\", #-2
		\"black\", #-1
		\"white\", #0, 1
		\"orange4\",#6,7
		\"red4\", #8,9
		\"white\", #10
		\"wheat2\") #12	
	mybreakz = c(
		-99.5,
		-9.5, #-7 -9
		-7.5, #-5 -6
		-5.5, #-2
		-1.5, #-1
		-0.5, #0-4
		5.5,
		7.5, #6-9
		97.5, #98.5
		98.5, #99-100
		100.5)

	#colnames(dm) = c(seq(1,dim(dm)[2]-1),\"clust\")
	#dm\$y = seq(dim(dm)[1],1,-1)
	#clustcoor = as.data.frame(aggregate(dm\$y,by=list(dm\$clust),min))
	#colnames(clustcoor) = c(\"clust\",\"ymin\")
	#clustcoor\$ymax = as.data.frame(aggregate(dm\$y,by=list(dm\$clust),max))\$x
	#clustcoor\$xmin = 1
	#clustcoor\$xmax = dim(dm)[2]
	#dm = subset(dm,select=c(-clust))
	#print(\"DM:\")
	#print(head(dm\$y))

#	dm = melt(dm,id.vars=c(\"y\"))
#	dm\$x=as.numeric(as.character(dm\$variable))
#	dm\$value=as.factor(as.character(dm\$value))
#	print(head(dm))
#	print(plyr::count(dm\$value))
#	png(\"$outDirPIC/$fileName1.reorder.$INFO.png\",width=mywidth/3,height=myheight)
#	ggplot(dm,aes(x,y)) + geom_tile(aes(fill=value)) + theme_blank +
#	theme_bw() + theme(panel.grid=element_blank()) +
#	geom_rect(data=clustcoor,aes(x=xmin,y=ymin,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,color=as.factor(clust)),fill=rgb(0,0,0,alpha=0)) +
#	scale_fill_manual(values=c(
#\"-9\" = \"blue4\",
#\"-8\" = \"blue4\",
#\"-7\" = \"blue4\",
#\"-6\" = \"blue4\",
#\"-5\" = \"green4\",
#\"-4\" = \"green4\",
#\"-2\" = \"grey4\",
#\"-1\" = \"black\",
#\"0\" = \"white\",
#\"1\" = \"white\",
#\"2\" = \"white\",
#\"4\" = \"white\",
#\"5\" = \"white\",
#\"6\" = \"orange4\",
#\"7\" = \"orange4\",
#\"8\" = \"red4\",
#\"9\" = \"red4\"))
	png(type=\"cairo\",\"$outDirPIC/$fileName1.reorder.$INFO.png\",width=4000,height=4000)
   heatmap.3(dm,dendrogram=\"none\",Rowv=F,Colv=F,cluster.by.row=F,cluster.by.col=F,labCol=F,labRow=F,
   RowIndividualColors=rev(myclust),
   labRow.by.group=TRUE,
   color.FUN=function(x){mycolorz},
   breaks = mybreakz)
   dev.off()

	mydf = subset(df,select=c(-clust))
	mydf1Beg = 1
	mydf1End = as.integer((dim(mydf)[2] - 3)/2)
	mydf2Beg = mydf1End + 3
	mydf2End = dim(mydf)[2]
	mydf1 = mydf[,seq(mydf1Beg,mydf1End)]
	mydf2 = mydf[,seq(mydf2Beg,mydf2End)]
	colnames(mydf1) = seq(1,dim(mydf1)[2])
	colnames(mydf2) = seq(1,dim(mydf2)[2])
	if (dim(mydf1)[2] < dim(mydf2)[2]) {
		mydf2 = mydf2[,seq(1,dim(mydf1)[2])]
	}
	if (dim(mydf2)[2] > dim(mydf1)[2]) {
		mydf1 = mydf1[,seq(1,dim(mydf2)[2])]
	}
   BGNONT = mydf1
   BGNONT[BGNONT < 4 | BGNONT >= 8] = NA
   BGNONT[!is.na(BGNONT) & BGNONT < 6] = 0
   BGNONT[!is.na(BGNONT) & BGNONT >= 6] = 1
   BGNONT = apply(BGNONT,2,function(x)mean(x,na.rm=T))
   BGTEMP = mydf2
   BGTEMP[BGTEMP < 4 | BGTEMP >= 8] = NA
   BGTEMP[!is.na(BGTEMP) & BGTEMP < 6] = 0
   BGTEMP[!is.na(BGTEMP) & BGTEMP >= 6] = 1
   BGTEMP = apply(BGTEMP,2,function(x)mean(x,na.rm=T))
   mybgall = data.frame(NON_TEMPLATE=BGNONT,TEMPLATE=BGTEMP,pos=seq(1,length(BGNONT)))
	if(length(mybgall[is.na(mybgall\$NON_TEMPLATE),]\$NON_TEMPLATE) > 0) {
		mybgall[is.na(mybgall\$NON_TEMPLATE),]\$NON_TEMPLATE = 0
	}
	if(length(mybgall[is.na(mybgall\$TEMPLATE),]\$TEMPLATE) > 0) {
		mybgall[is.na(mybgall\$TEMPLATE),]\$TEMPLATE = 0
	}
   mybgall=melt(mybgall,id.vars=c(\"pos\"))
   mybg = data.frame(type=c(\"NON_TEMPLATE\",\"TEMPLATE\"),bgvalue=c(mean(BGNONT,na.rm=T),mean(BGTEMP,na.rm=T)))
   colnames(mybgall) = c(\"variable\",\"type\",\"bgvalue\")
	
	mybgall\$value = mybgall\$bgvalue
	if(length(mybg[is.na(mybg\$bgvalue),]\$bgvalue) > 0) {
		mybg[is.na(mybg\$bgvalue),]\$bgvalue = 0
	}
	mydf1[mydf1 < 4] = NA
	mydf2[mydf2 < 4] = NA
   mydf1[!is.na(mydf1) & mydf1 < 6] = 0
   mydf1[!is.na(mydf1) & mydf1 >= 6 & mydf1 < 10] = 1
   mydf2[!is.na(mydf2) & mydf2 < 6] = 0
   mydf2[!is.na(mydf2) & mydf2 >= 6 & mydf2 < 10] = 1
   mydf1\$clust = df\$clust
   mydf1\$type = \"NON_TEMPLATE\"
   mydf2\$clust = df\$clust
   mydf2\$type = \"TEMPLATE\"
	print(plyr::count(mydf1\$clust))
	print(plyr::count(mydf2\$clust))
	mydf = rbind(mydf1,mydf2)
	mydf = melt(mydf,id.vars=c(\"clust\",\"type\"))
#	mydf=mydf[mydf\$value != -0.25,]
	#mydf\$variable = as.numeric(as.character(mydf\$variable))
	mydf = as.data.frame(aggregate(mydf\$value,by=list(mydf\$clust,mydf\$type,mydf\$variable),function(x)mean(x,na.rm=T)))
#,trim=0.05,na.rm=T)))
	colnames(mydf) = c(\"clust\",\"type\",\"variable\",\"value\")

 	clust = data.frame(beg=begz,end=endz,clust=clustz)
	print(\"Clust Before:\")
	print(clust)
	clusthere = unique(mydf\$clust)
	clust = clust[clust\$clust \%in\% clusthere,]
	print(\"Clust After:\")
	print(clust)

   if (dim(mydf[is.na(mydf\$value),])[1] > 0) {
      mydf[is.na(mydf\$value),]\$value = 0
   }
   print(summary(as.factor(mydf\$clust)))

	#subtract each value with average background
	mydf = merge(mydf,mybg,by=c(\"type\"),all=T)

	#subtract each value with each bacground at each positon
	#mydf = merge(mydf,mybg,by=c(\"type\",\"variable\"),all=T)
   print(summary(as.factor(mydf\$clust)))
	mydf\$value = mydf\$value - mydf\$bgvalue
	if (dim( mydf[mydf\$value < 0,])[1] > 0) {
      mydf[mydf\$value < 0,]\$value = 0
   }
	#if (length(	mydf[mydf\$value < 0,]\$value) > 0) {
	#	mydf[mydf\$value < 0,]\$value = 0
#	}
	#subtract each value with each bacground at each positon
	
	mydf[mydf\$type == \"TEMPLATE\",]\$value = 	-1 * mydf[mydf\$type == \"TEMPLATE\",]\$value
	mybg0 = mybgall
	mybg0[mybg0\$type == \"TEMPLATE\",]\$bgvalue = -1 * mybg0[mybg0\$type == \"TEMPLATE\",]\$bgvalue
   mybg0\$clust = -1
   mybg0\$value = mybg0\$bgvalue
	print(colnames(mydf))
	print(colnames(mybg0))
   mydf = rbind(mydf,mybg0)
	mydf\$variable = as.numeric(as.character(mydf\$variable))
	mydf\$value = as.numeric(as.character(mydf\$value))

	pdf(\"$outDirPIC/$fileName1.conv.$INFO.pdf\",height=$maxclust,width=14)
	p1 = ggplot(mydf,aes(variable,value)) + 
	geom_rect(data=clust,aes(x=beg,y=0,xmin=beg,xmax=end,ymin=-0.5,ymax=1),size=0.5) +
	geom_line(aes(color=type),size=1) +
	theme_bw() +
	facet_grid(clust~.) + coord_cartesian(ylim=c(-0.5,1)) +
	theme(panel.grid=element_blank()) + ylab(\"Fraction Conversion\") + xlab(\"bp from start of Amplicon\") +
	ggtitle(\"Cluster Conversion\")
	print(p1)
	dev.off()
	";

open (my $outDoR, ">", "$outDir/$fileName1.do.R");
print $outDoR "
library(ggplot2)
library(reshape2)
myformat = function(x){as.integer(x*1000+0.5)/1000;return(x)}
orig = read.table(\"$outDirPIC/$fileName1.ORIG.ClustConv\")
shuf = read.table(\"$outDirPIC/$fileName1.SHUF.ClustConv\")
colnames(orig) = c(\"clust\",\"type\",\"Value\",\"Total\",\"shuf\",\"Rep\")
colnames(shuf) = colnames(orig)
shuf = shuf[shuf\$shuf == \"SHUF\",]
shuf = shuf[shuf\$clust \%in\% orig\$clust,]
dfALL = rbind(orig,shuf)
df = dfALL
conv = df[df\$type == \"BEG\"|df\$type == \"END\",]
df = df[df\$type != \"BEG\" & df\$type != \"END\",]
df = rbind(df,data.frame(clust=99,type=\"ALL\",Value=df[df\$shuf == \"MAIN\",]\$Total[1] - sum(df[df\$shuf==\"MAIN\",]\$Value),Total=df[df\$shuf == \"MAIN\",]\$Total[1],shuf=\"MAIN\",Rep=0))
df\$perc = as.integer(df\$Value/df\$Total*10000)/100

dfMAIN = df[df\$shuf == \"MAIN\",]
dfORIG = df[df\$shuf == \"ORIG\",]
dfSHUF = df[df\$shuf == \"SHUF\",]
dfclust = unique(df\$clust)
p1 = data.frame(clust=1,p=0,mainval=1,origval=1,shufval=1,shufsd=0,shufse=0,nge=0,nle=0,ntot=0)
for (i in 1:length(dfclust)) {
	myclust = dfclust[i]
	MAIN = dfMAIN[dfMAIN\$clust == myclust,]\$perc
	ORIG = dfORIG[dfORIG\$clust == myclust,]\$perc
	SHUF = dfSHUF[dfSHUF\$clust == myclust,]
	#SHUF = dfSHUF[dfSHUF\$clust != 99,]
	Rep = unique(SHUF\$Rep)
	ntot = 0
	nge = 0
	nle = 0
	SHUFS = c(0)
	for (j in 1:length(Rep)) {
		myrep = Rep[j]
		SHUF2 = mean(SHUF[SHUF\$Rep == Rep[j],]\$perc,trim=0.05)
		SHUFS[j] = SHUF2
		if (SHUF2 >= ORIG) {nge = nge + 1}
		if (SHUF2 <= ORIG) {nle = nle + 1}
		ntot = ntot + 1
		#		print(paste(nge,nle,ntot,ORIG,SHUF2))
	}
	myp = 1
	myshufval = myformat(mean(SHUFS,trim=0.05))
	mysd = sd(SHUFS)
	myse = sd(SHUFS)/sqrt(length(SHUFS))
	if (nge > nle) {myp = (nle + 1) / (ntot + 1)}
	if (nle > nge) {myp = (nge + 1) / (ntot + 1)}
	#	print(paste(myclust,myp))
p1 = rbind(p1,data.frame(clust=myclust,p=myp,mainval=MAIN,origval=ORIG,shufval=myshufval,shufsd=myformat(mysd),shufse=myformat(myse),nge=nge,nle=nle,ntot=ntot))
}
p1 = p1[-1,]
p1\$p = as.numeric(as.character(p1\$p))
p1\$pval = \"\"
if (dim(p1[p1\$p <= 0.1,])[1] > 0) {
	p1[p1\$p <= 0.1,]\$pval = \"*\"
}
if (dim(p1[p1\$p <= 0.05,])[1] > 0) {
	p1[p1\$p <= 0.05,]\$pval = \"**\"
}
if (dim(p1[p1\$p <= 0.01,])[1] > 0) {
	p1[p1\$p <= 0.01,]\$pval = \"***\"
}


convORIG = conv[conv\$shuf == \"ORIG\",]
convSHUF = conv[conv\$shuf == \"SHUF\",]
convtype = unique(conv\$type)
convclust = unique(conv\$clust)
p2 = data.frame(clust=1,p=0,type=\"BEG\",orig=1,shuf=1,shufsd=0,shufse=0,nge=0,nle=0,ntot=0)
for (h in 1:length(convtype)) {
	mytype = convtype[h]
	for (i in 1:length(convclust)) {
		myclust = convclust[i]
		ORIG = mean(convORIG[convORIG\$clust == myclust & convORIG\$type == mytype,]\$Value,trim=0.05)
		SHUF = convSHUF[convSHUF\$clust == myclust & convSHUF\$type == mytype,]
		#SHUF = convSHUF[convSHUF\$clust != 99 & convSHUF\$type == mytype,]
		Rep = unique(SHUF\$Rep)
		ntot = 0
		nge = 0
		nle = 0
		SHUFS = c(0)
		for (j in 1:length(Rep)) {
#			print(paste(h,convtype[h],i,convclust[i],j,Rep[j]))
			myrep = Rep[j]
			SHUF2 = mean(SHUF[SHUF\$Rep == Rep[j],]\$Value,trim=0.05)
			SHUFS[j] = SHUF2
			if (SHUF2 >= ORIG) {nge = nge + 1}
			if (SHUF2 <= ORIG) {nle = nle + 1}
			ntot = ntot + 1
			#		print(paste(nge,nle,ntot,ORIG,SHUF2))
		}
		myp = 1
		mysd = sd(SHUF\$Value)
		myse = sd(SHUF\$Value)/sqrt(length(SHUF\$Value))
		if (nge > nle) {myp = (nle + 1) / (ntot + 1)}
		if (nle > nge) {myp = (nge + 1) / (ntot + 1)}
		#	print(paste(myclust,myp))
		p2 = rbind(p2,data.frame(clust=myclust,p=myp,type=mytype,orig=ORIG,shuf=myformat(mean(SHUFS,trim=0.05)),shufsd=myformat(mysd),shufse=myformat(myse),nge=nge,nle=nle,ntot=ntot))
	}
}
p2 = p2[-1,]
p2\$p = as.numeric(as.character(p2\$p))
p2\$pval = \"\"
if (dim(p2[p2\$p <= 0.1,])[1] > 0) {
	p2[p2\$p <= 0.1,]\$pval = \"*\"
}
if (dim(p2[p2\$p <= 0.05,])[1] > 0) {
	p2[p2\$p <= 0.05,]\$pval = \"**\"
}
if (dim(p2[p2\$p <= 0.01,])[1] > 0) {
	p2[p2\$p <= 0.01,]\$pval = \"***\"
}

pdf(\"$outDirPIC/$fileName1.BOXPLOT.pdf\",width=14,height=7)
ggplot(df,aes(as.factor(clust),perc)) + geom_boxplot(aes(color=shuf),fill=rgb(0,0,0,0),outlier.shape=NA) + theme_bw() +theme(panel.grid=element_blank()) +
	geom_text(data=p1,aes(x=as.factor(clust),y=100,label=pval))

ggplot(conv,aes(as.factor(clust),Value)) + geom_boxplot(aes(color=shuf),fill=rgb(0,0,0,0),outlier.shape=NA) + theme_bw() +theme(panel.grid=element_blank()) + facet_grid(type~.) +
	geom_text(data=p2,aes(x=as.factor(clust),y=100,label=pval))
dev.off()

write.table(p1,file=\"$outDirPIC/$fileName1.BOXPLOT.p1\",quote=F,row.names=F,col.names=T,sep=\"\\t\")
write.table(p2,file=\"$outDirPIC/$fileName1.BOXPLOT.p2\",quote=F,row.names=F,col.names=T,sep=\"\\t\")
# 
# df1 = plyr::count(df,vars=c(\"clust\",\"type\",\"shuf\",\"Rep\"))
# df2 = plyr::count(df,vars=c(\"type\",\"shuf\"))
# df = merge(df1,df2,by=c(\"type\",\"shuf\"),all=T)
# df\$perc = as.integer(df\$freq.x / df\$freq.y*10000)/100
	";
	close $outDoR;
#	system("Rscript $outDir/$fileName1.R");
	my ($pdffullPath) = getFullpath("$outDir/$fileName1.conv.pdf");
	my ($pngfullPath) = getFullpath("$outDir/$fileName1.reorder.png");
	print("\nscp mitochi\@crick.cse.ucdavis.edu:$pdffullPath ./\n");
	print("scp mitochi\@crick.cse.ucdavis.edu:$pngfullPath ./\n\n");
}

sub parse_footPeak_logFile {
   my ($footPeakFolder) = @_;
#	die "Index File not found in $LCY$footPeakFolder/footPeak_logFile.txt$N\n" if not defined $indexFile or not -e $indexFile;
#	print "DINEX = $indexFile\n";
	my $coor;
	open (my $in , "<", "$footPeakFolder/footPeak_logFile.txt") or die;
	while (my $line = <$in>) {
		chomp($line);
		if ($line =~ /^Index File =/) {
			$line = <$in>;
			while ($line =~ /^def=/) {
				chomp($line);
				my ($def, $chr, $beg, $end, $gene, $value, $strand) = $line =~ /^def=(.+), coor=(.+), (\d+), (\d+), (.+), (\d+), (.)$/;
				$coor->{uc($gene)}{chr} = $chr;
				$coor->{uc($gene)}{beg} = $beg;
				$coor->{uc($gene)}{end} = $end;
				$coor->{uc($gene)}{strand} = $strand;
				$coor->{uc($gene)}{coor} = "$chr\t$beg\t$end\t$gene\t$value\t$strand";
				$line = <$in>;
				print "Gnee=$LCY$gene$N, coor=$chr,$beg,$end,$strand\n" if not defined $geneWant;
				print "Gnee=$LCY$gene$N, coor=$chr,$beg,$end,$strand\n" if defined $geneWant and uc($geneWant) eq uc($gene);
			}
			last;
		}
	}
	return $coor;
#   my $paramsFile = "$footFolder/.PARAMS";
#   my $geneIndexFile;
#   my @parline = `cat $paramsFile`;
#   foreach my $parline (@parline) {
#      if ($parline =~ /footLoop.pl,geneIndexFile,/) {
#         ($geneIndexFile) = $parline =~ /geneIndexFile,(.+)$/;
#      }
#   }
#   return($geneIndexFile);
}

sub parse_clustFile {
	my ($peakFileNTEMP, $coor, $srand, $info) = @_;
	my ($footFolder, $fileName) = getFilename($peakFileNTEMP,"folderfull");
	my ($gene, $strand, $conv) = $fileName =~ /gene(.+)_(Pos|Neg|Unk)_.+_(CG|GC|GH|CH)/; die if not defined $gene or not defined $strand or not defined $conv;
	if (defined $opt_S) {
		$conv = $conv eq "CG" ? "GC" : $conv eq "GC" ? "CG" : $conv eq "GH" ? "CH" : $conv eq "CH" ? "GH" : die "Can't find conv for $fileName!\n";
		$strand = $strand eq "Pos" ? "Neg" : $strand eq "Neg" ? "Pos" : die "Can't find strand for $fileName!\n";
	}
	my ($clustFile) = <$footFolder/../FOOTCLUST/CLUST_LOCAL/*gene$gene*$strand*$conv*local.bed.indiv.clust>;
	print "clustFile $clustFile doesn't exist!\n
ls $footFolder/../FOOTCLUST/CLUST_LOCAL/*gene$gene*$strand*$conv*local.bed.indiv.clust
\n" and return if not defined $clustFile;
	print "Getting clust from peakFile=$peakFileNTEMP\nclustFile=$YW$clustFile$N\n\n";
	open (my $in, "<", $clustFile);
	my $pos;
	my $clustCount = 1;
	my %mydf;
	my $lc = 0;
	my $lastclust = "INIT";
	my %temp;
	while (my $line = <$in>) {
		chomp($line);
		$lc ++;
		my ($name, $beg, $end, $a, $b, $clust) = split("\t", $line);
		my ($date, $num, $num2) = $name =~ /^(\d{6})(\d{6})(\d+)$/;
		$num2 =~ s/0$//;
		$temp{$clust}{beg} += $beg;
		$temp{$clust}{end} += $end;
		$temp{$clust}{len} += $end - $beg;
		$temp{$clust}{total} ++;
		$temp{$clust}{name}{$name} = 1;
		$lastclust = $clust;
	}
	my $coorBeg = $coor->{uc($gene)}{beg};
	my $coorEnd = $coor->{uc($gene)}{end};
	die if not defined $coorBeg;
	my $coorLen = $coorEnd - $coorBeg - 200;
#	$coorBeg += 100;
	$coorBeg = 100;
	$srand = (defined $srand and $srand eq -1) ? 42 : defined $srand ? $srand : 42;
	srand($srand);
	foreach my $clust (sort {$temp{$b}{len} <=> $temp{$a}{len}} keys %temp) {
		if (defined $info and $info ne "SHUF") {
			print "$YW CLUST=$clust$N, total = $LGN$temp{$clust}{total}$N\n";
			print "$LRD(Skipped)$N\n" if $temp{$clust}{total} < 20;
			print "$LGN(USED)$N\n" if $temp{$clust}{total} >= 20;
		}
		
next if $temp{$clust}{total} < 20;
		if (not defined $opt_R) {
			$pos->{$clust}{beg} = int($temp{$clust}{beg}/$temp{$clust}{total}+0.5);
			$pos->{$clust}{end} = int($temp{$clust}{end}/$temp{$clust}{total}+0.5);
		}
		else {
			my $rand1 = find_random($coorBeg, $coorLen, \%temp);
			my $rand2 = find_random($coorBeg, $coorLen, \%temp, $rand1);
			$pos->{$clust}{beg} = $rand1 < $rand2 ? $rand1 : $rand2;
			$pos->{$clust}{end} = $rand1 < $rand2 ? $rand2 : $rand1;
			print "CLUSTER $clust: $pos->{$clust}{beg}\t$pos->{$clust}{end}\n";
		}
		foreach my $namez (sort keys %{$temp{$clust}{name}}) {
			$pos->{$clust}{name}{$namez} = 1;
		}
		$clustCount ++;
	}
	my $maxclust = (keys %{$pos}) + 1; # not in any cluster
	return ($pos, $maxclust);
}

sub find_random {
	my ($coorBeg, $coorLen, $temp, $pos) = @_;
	#print "#Find random\n";
	$pos = -100 if not defined $pos;
	my %temp = %{$temp};
	for (my $i = 0; $i < 1000; $i++) {
		my $bad = 0;
		my $rand = int(rand($coorLen) + $coorBeg);
		foreach my $clust (sort keys %temp) {
			my $beg = int($temp{$clust}{beg}/$temp{$clust}{total}+0.5);
			my $end = int($temp{$clust}{end}/$temp{$clust}{total}+0.5);
			if (intersect($rand-$window, $rand+$window, $beg-250,$beg+250) == 1) {
				$bad ++;
				#print "i=$i: clust=$clust, beg = $beg, end=$end, rand=$rand, pos=$pos, bad=$bad\n";
			}
			if (intersect($rand-$window, $rand+$window, $end-250,$end+250) == 1) {
				$bad ++;
				#print "i=$i: clust=$clust, beg = $beg, end=$end, rand=$rand, pos=$pos, bad=$bad\n";
			}
			#if ($i == 0)  {print "clust=$clust, beg=$beg, end=$end, pos=$pos, bad=$bad, rand=$rand\n";}
		}
		$bad ++ if intersect($rand-$window, $rand+$window, $pos-$window,$pos+$window) == 1;
		return $rand if $bad == 0;
	}
	for (my $i = 0; $i < 1000; $i++) {
		my $bad = 0;
		my $rand = int(rand($coorLen) + $coorBeg);
		foreach my $clust (sort keys %temp) {
			my $beg = int($temp{$clust}{beg}/$temp{$clust}{total}+0.5);
			my $end = int($temp{$clust}{end}/$temp{$clust}{total}+0.5);
			if (intersect($rand-$window, $rand+$window, $beg-150,$beg+150) == 1) {
				$bad ++;
				#print "i=$i: clust=$clust, beg = $beg, end=$end, rand=$rand, pos=$pos, bad=$bad\n";
			}
			if (intersect($rand-$window, $rand+$window, $end-150,$end+150) == 1) {
				$bad ++;
				#print "i=$i: clust=$clust, beg = $beg, end=$end, rand=$rand, pos=$pos, bad=$bad\n";
			}
			#if ($i == 0)  {print "clust=$clust, beg=$beg, end=$end, pos=$pos, bad=$bad, rand=$rand\n";}
		}
		$bad ++ if intersect($rand-$window, $rand+$window, $pos-$window,$pos+$window) == 1;
		return $rand if $bad == 0;
	}
	for (my $i = 0; $i < 1000; $i++) {
		my $bad = 0;
		my $rand = int(rand($coorLen) + $coorBeg);
		foreach my $clust (sort keys %temp) {
			my $beg = int($temp{$clust}{beg}/$temp{$clust}{total}+0.5);
			my $end = int($temp{$clust}{end}/$temp{$clust}{total}+0.5);
			if (intersect($rand-$window, $rand+$window, $beg-100,$beg+100) == 1) {
				$bad ++;
				#print "i=$i: clust=$clust, beg = $beg, end=$end, rand=$rand, pos=$pos, bad=$bad\n";
			}
			if (intersect($rand-$window, $rand+$window, $end-100,$end+100) == 1) {
				$bad ++;
				#print "i=$i: clust=$clust, beg = $beg, end=$end, rand=$rand, pos=$pos, bad=$bad\n";
			}
			#if ($i == 0)  {print "clust=$clust, beg=$beg, end=$end, pos=$pos, bad=$bad, rand=$rand\n";}
		}
		$bad ++ if intersect($rand-$window, $rand+$window, $pos-$window,$pos+$window) == 1;
		return $rand if $bad == 0;
	}
	for (my $i = 0; $i < 1000; $i++) {
		my $bad = 0;
		my $rand = int(rand($coorLen) + $coorBeg);
		foreach my $clust (sort keys %temp) {
			my $beg = int($temp{$clust}{beg}/$temp{$clust}{total}+0.5);
			my $end = int($temp{$clust}{end}/$temp{$clust}{total}+0.5);
			if (intersect($rand-$window, $rand+$window, $beg-50,$beg+50) == 1) {
				$bad ++;
				#print "i=$i: clust=$clust, beg = $beg, end=$end, rand=$rand, pos=$pos, bad=$bad\n";
			}
			if (intersect($rand-$window, $rand+$window, $end-50,$end+50) == 1) {
				$bad ++;
				#print "i=$i: clust=$clust, beg = $beg, end=$end, rand=$rand, pos=$pos, bad=$bad\n";
			}
			#if ($i == 0)  {print "clust=$clust, beg=$beg, end=$end, pos=$pos, bad=$bad, rand=$rand\n";}
		}
		$bad ++ if intersect($rand-$window, $rand+$window, $pos-$window,$pos+$window) == 1;
		return $rand if $bad == 0;
	}
	die "CAnnot find random for $coorBeg $coorLen\n";
}

sub printclust {
	my ($pos) = @_;
	(print "UNDEF POS\n" and return )if not defined $pos;
	foreach my $clust (sort keys %{$pos}) {
		print "$clust = $pos->{$clust}{beg}\n";
	}
	print "DONE\n";
}

sub process_PEAKNOPKFILE {
	my ($peakFileNTEMP, $pos, $data, $outClust, $info, $rep) = @_;	
	my $undef = 0;
	printclust($pos);
	my $maxclust = (keys %{$pos}) + 1; # not in any cluster
	#$data = $dataz 
	my $dataz;
	if (defined $info and $info =~ /^(INFO|ORIG)$/) {
		$dataz = $data;
	}
#	if (not defined $info or (defined $info and $info ne "SHUF")) {
#		$dataz = $data;
#	}
	print "Doing $peakFileNTEMP (maxclust=$YW$maxclust$N)\n";
	my $lc = 0;
	my @head;
	my %seq;
	my %conv;
	open (my $in1, "<", $peakFileNTEMP) or die "Cannot read from $peakFileNTEMP: $!\n";
	my %clustz;
	while (my $line = <$in1>) {
		$lc ++;
		chomp($line);
		$line =~ s/\t+/\t/g;
		my @arr = split("\t", $line);
		if ($lc == 1) {
			@head = @arr;
			for (my $i = 5; $i < @arr; $i++) {
				my $pos = $i - 5;
				if ($i < @arr-1 and $arr[$i] eq "C" and $arr[$i+1] eq "G") {
					$seq{C}{$i} = 2;
				}
				elsif ($i > 0 and $arr[$i] eq "G" and $arr[$i-1] eq "C") {
					$seq{G}{$i} = 2;
				}
				else {
					$seq{C}{$i} = 1 if $arr[$i] eq "C";
					$seq{G}{$i} = 1 if $arr[$i] eq "G";
				}
			}
			next;
		}
		my $readName = $arr[0];
		my ($date, $num1, $num2) = $readName =~ /.+m(\d{3,6})_(\d{6}).+\/(\d+)\//;
		die if not defined $num2;
		for (my $k = 0; $k < @arr; $k++) {
			die "arr k=$k is empty\n" if not defined $arr[$k];
		}
#		print "readname = $readName\n";
		$dataz->{$readName}{value} = join("\t", @arr[5..@arr-1]);
# -------[0  beg  1]----------------[0   end   1] ---------
		my $goodclust;# = 0;
		if (defined $info and $info eq "INFO") {
			foreach my $clust (sort keys %{$pos}) {
				foreach my $name (sort keys %{$pos->{$clust}{name}}) {
					if ($name =~ /$date.*$num1.*$num2/) {
						$goodclust = $clust;
						$dataz->{$readName}{clust}{$clust} = 999;
#						$dataz->{$readName}{clustbegtot}{$clust} = 999;
#						$dataz->{$readName}{clustendtot}{$clust} = 999;
						$dataz->{$readName}{clustbeg}{$clust} = 1;
						$dataz->{$readName}{clustend}{$clust} = 1;						
						last;
					}
				}
				last if defined $goodclust;
			}
			#for (my $i = $end0; $i <= $end1; $i++) {
			if (not defined $goodclust) {
#				print "Read = $readName, date=$date, num=$num1 2=$num2, clust=UNK\n" if $undef < 10;
#				print "More undef detected, more than 10\n" if $undef == 10;
				$undef ++;
				$dataz->{$readName}{clust}{$maxclust} = 1;
				next;
			}
			$clustz{$goodclust} ++;
		}
		else {
			foreach my $clust (sort keys %{$pos}) {
			$dataz->{$readName}{clust}{$clust} = $maxclust if not defined $dataz->{$readName}{clust}{$clust};
			my $beg0 = $pos->{$clust}{beg} - $window; $beg0 = 0 if $beg0 <= 0;
			my $beg1 = $pos->{$clust}{beg} + $window; $beg1 = @arr - 5 if $beg1 > @arr - 5;
			my $end0 = $pos->{$clust}{end} - $window; $end0 = 0 if $end0 <= 0;
			my $end1 = $pos->{$clust}{end} + $window; $end1 = @arr - 5 if $end1 > @arr - 5;
			my ($convertedBeg, $convertedEnd, $totalBeg, $totalEnd) = (0,0,0,0);
			for (my $i = $beg0; $i <= $beg1; $i++) {
				$convertedBeg ++ if $arr[$i] =~ /^[6789]$/;
				$totalBeg ++ if $arr[$i] =~ /^[456789]$/;
			}
			for (my $i = $end0; $i <= $end1; $i++) {
				$convertedEnd ++ if $arr[$i] =~ /^[6789]$/;
				$totalEnd ++ if $arr[$i] =~ /^[456789]$/;
			}
			my $convertedBest =  $convertedBeg > $convertedEnd ? $convertedBeg : $convertedEnd;
			my $totalBest =  $convertedBeg > $convertedEnd ? $totalBeg : $totalEnd;
			$convertedBest = $totalBest == 0 ? 0 : int(1000*$convertedBest/$totalBest+0.5)/1000;
#			$best = int(1000*$best+0.5)/1000;
	#		$dataz->{$readName}{clust}{$clust} = $best;
			$dataz->{$readName}{clust}{$clust} = $convertedBeg + $convertedEnd;
#			$dataz->{$readName}{clustbegtot}{$clust} = $totalBeg;
#			$dataz->{$readName}{clustendtot}{$clust} = $totalEnd;
			$dataz->{$readName}{clustbeg}{$clust} = $convertedBeg == 0 ? 0 : $totalBeg == 0 ? $convertedBeg : int(100*$convertedBeg/$totalBeg+0.5);
			$dataz->{$readName}{clustend}{$clust} = $convertedEnd == 0 ? 0 : $totalEnd == 0 ? $convertedEnd : int(100*$convertedEnd/$totalEnd+0.5);
			#for (my $i = $end0; $i <= $end1; $i++) {
			#	$dataz->{$readName}{clust}{$clust} ++ if $arr[$i] =~ /^[6789]$/;
			#}
			}
			
		}
	}
	print "lc=$lc, undef=$undef\n";
	my $totalshuf = $lc - 1;
	my $totalorig = 0;
	close $in1;
#		foreach my $read (sort {$dataz->{$strandz}{$b}{sort}{1} <=> $dataz->{$strandz}{$a}{sort}{1}} keys %{$dataz->{$strandz}}) {
#			my $clust = $dataz->{$strandz}{$read}{sort}{1};
	my %totalz;
	foreach my $clust (sort keys %{$pos}) {
		$totalz{$clust} = 0;
		$clustz{$clust} = 0 if not defined $clustz{$clust};
	}
	if (defined $info and $info eq "INFO") {
		foreach my $clust (sort keys %clustz) {
			print "ORIG GOODCLUST cluster $clust = $clustz{$clust}\n" if defined $info and $info eq "INFO";
			print $outClust "$clust\tALL\t$clustz{$clust}\t$totalshuf\tMAIN\t0\n";
		}
	}
	foreach my $read (sort keys %{$dataz}) {
		my $count = 0;
		foreach my $clust (sort {$dataz->{$read}{clust}{$b} <=> $dataz->{$read}{clust}{$a}} keys %{$dataz->{$read}{clust}}) {
			$count ++;
			if (defined $info and $info eq "INFO") {
				$dataz->{$read}{sort}{$count} = $clust;
				$totalz{$clust} ++;
	}
			elsif ($count == 1 and $dataz->{$read}{clustbeg}{$clust} > 5 and $dataz->{$read}{clustend}{$clust} > 5) {
				$totalz{$clust} ++;
				$dataz->{$read}{sort}{$count} = $clust;
				if (not defined $opt_R or (defined $opt_R and defined $info and $info ne "ORIG")) {
					print $outClust "$clust\tBEG\t$dataz->{$read}{clustbeg}{$clust}\t$totalshuf\t$info\t$rep\n";
					print $outClust "$clust\tEND\t$dataz->{$read}{clustend}{$clust}\t$totalshuf\t$info\t$rep\n";
				}
			}
			elsif ($count == 1) {
				$dataz->{$read}{sort}{1} = 99;
				$totalz{99} ++;
				if (not defined $opt_R or (defined $opt_R and defined $info and $info ne "ORIG")) {
					print $outClust "99\tBEG\t$dataz->{$read}{clustbeg}{$clust}\t$totalshuf\t$info\t$rep\n";
					print $outClust "99\tEND\t$dataz->{$read}{clustend}{$clust}\t$totalshuf\t$info\t$rep\n";
				}
			}
			last;
		}
		foreach my $clust (sort {$a <=> $b} keys %{$dataz->{$read}{clust}}) {
		}
	}
	my $sumz = 0;
	if (defined $info and $info ne "INFO") {# "SHUF") {
		foreach my $clust2 (sort {$a <=> $b} keys %totalz) {
			$totalz{$clust2} = 0 if not defined $totalz{$clust2};
			$sumz += $totalz{$clust2};
		}
		foreach my $clust2 (sort {$a <=> $b} keys %totalz) {
			my $perc = int($totalz{$clust2} / $sumz * 10000)/100;
			print $outClust "$clust2\tALL\t$totalz{$clust2}\t$sumz\t$info\t$rep\n";
			print "SHUF cluster $clust2 = $totalz{$clust2}\n";
		}
		print "Total = $sumz\n";
		return;
	}
	else {
		return($dataz);
	}
}

__END__
	$pos->{2}{beg} = 649;
	$pos->{3}{beg} = 829;
	$pos->{4}{beg} = 1084;
	$pos->{5}{beg} = 1056;
	$pos->{6}{beg} = 1232;
	$pos->{7}{beg} = 1328;
	$pos->{2}{end} = 862;
	$pos->{3}{end} = 1102;
	$pos->{4}{end} = 1246;
	$pos->{5}{end} = 1444;
	$pos->{6}{end} = 1404;
	$pos->{7}{end} = 1444;
#pFC19FIXED  649   862   PFC8_SNRPN_REVERSE.2 1000  -
#pFC19FIXED  829   1102  PFC8_SNRPN_REVERSE.3 1000  -
#pFC19FIXED  1084  1246  PFC8_SNRPN_REVERSE.4 1000  -
#pFC19FIXED  1056  1444  PFC8_SNRPN_REVERSE.5 590   -
#pFC19FIXED  1232  1404  PFC8_SNRPN_REVERSE.6 760   -
#pFC19FIXED  1328  1444  PFC8_SNRPN_REVERSE.7 1000  -
	print "Total Read = $lc\n";
	return \%{$pos};

}
__END__



open (my $in2, "<", $nopkFileNTEMP) or die "Cannot read from $nopkFileNTEMP: $!\n";
while (my $line = <$in2>) {
	$lc ++;
	chomp($line);
	$line =~ s/[\t]+/\t/g;
	my @arr = split("\t", $line);
	if ($lc == 1) {
		@head = @arr;
		next;
	}
	for (my $k = 5; $k < @arr; $k++) {
		die "died at $nopkFileNTEMP: k=$k value=$arr[$k]\n" if $arr[$k] !~ /^\d+$/;
	}
	my $readName = $arr[0];
	$data->{$readName}{value} = join("\t", @arr[5..@arr-1]);
#	for (my $i = 0; $i < @arr; $i++) {
#		my $head = $head[$i]; $head = "$LPR UNKNOWN$N" if not defined $head;
#	print "$head=$arr[$i]\n";
	foreach my $clust (sort keys %{$pos}) {
		$data->{$readName}{clust}{$clust} = 0 if not defined $data->{$readName}{clust}{$clust};
		my $beg0 = $pos->{$clust}{beg} - $window; $beg0 = 0 if $beg0 <= 0;
		my $beg1 = $pos->{$clust}{beg} + $window; $beg1 = @arr - 5 if $beg1 > @arr - 5;
		my $end0 = $pos->{$clust}{end} - $window; $end0 = 0 if $end0 <= 0;
		my $end1 = $pos->{$clust}{end} + $window; $end1 = @arr - 5 if $end1 > @arr - 5;
		my $value = 0; my $total1 = 0;
		for (my $i = $beg0; $i <= $beg1; $i++) {
			$value ++ if $arr[$i] =~ /^[6789]$/;
			$total1 ++ if $arr[$i] =~ /^[456789]$/;
		}
		my $beg2 = $value;#$total == 0 ? 0 : $value / $total;
		$value = 0; my $total2 = 0;
		for (my $i = $end0; $i <= $end1; $i++) {
			$value ++ if $arr[$i] =~ /^[6789]$/;
			$total2 ++ if $arr[$i] =~ /^[456789]$/;
		}
		my $end2 = $value;#$total == 0 ? 0 : $value / $total;
		my $best =  $beg2 > $end2 ? $beg2 : $end2;#/$total : $end2/$total;#total == 0 ? 0 : int(1000*$value / $total+0.5)/1000;
		my $besttot =  $beg2 > $end2 ? $total1 : $total2;
		$best = $besttot == 0 ? 0 : int(1000*$best/$besttot+0.5)/1000;
		$best = int(1000*$best+0.5)/1000;
#		$data->{$readName}{clust}{$clust} = $best;
		$data->{$readName}{clust}{$clust} = $beg2 + $end2;
#		$data->{$readName}{clustbegtot}{$clust} = $total1;
#		$data->{$readName}{clustendtot}{$clust} = $total2;
		$data->{$readName}{clustbeg}{$clust} = $beg2 == 0 ? 0 : $total1 == 0 ? $beg2 : int(100*$beg2/$total1+0.5);
		$data->{$readName}{clustend}{$clust} = $end2 == 0 ? 0 : $total1 == 0 ? $end2 : int(100*$end2/$total1+0.5);
#		$data->{$readName}{clustbeg}{$clust} = $beg2 == 0 ? 0 : int(100*$beg2/$total1+0.5);
#		$data->{$readName}{clustend}{$clust} = $end2 == 0 ? 0 : int(100*$end2/$total1+0.5);
	}
}
close $in2;


OUT0

__END__
	mycolorz = c(
\"blue4\", #-2
\"black\", #-1
\"white\", #0, 1
\"white\", #4
\"red2\", #6
\"red2\", #7
\"red2\", #8
\"red2\", #9
\"white\", #10
\"wheat2\", #11
\"wheat2\") #12
	mybreakz = c(
-2.5
-1.5, #1
-0.5, #2
2.5, #3
4.5, #4
5.5, #5
6.5, #6
7.5, #7
8.5, #8
9.5, #9
98.5, #10
99.5, #11
100.5) #12


"blue4", #-2
"black", #-1
"white", #0, 1
"red2", #6
"white", #10
"wheat2") #12
__END__
	ggplot(mybg0,aes(variable,bgvalue)) + geom_line(aes(color=type),size=1) +
	theme_bw() +
	coord_cartesian(ylim=c(-0.5,1)) +
	theme(panel.grid=element_blank()) + ylab(\"Fraction Conversion\") + xlab(\"bp from start of Amplicon\") +
	ggtitle(\"Background Non-Peak Conversion\")
	dev.off()
	

__END__
			for (my $i = @FOR; $i < @REV; $i++) {
				if (@FOR == 0) {
					$FOR[$i] = $grey;
					$print{$clust}{FOR}{name}[$i] = "grey$greycount";
					$greycount ++;
				}
				else {
					my $random = int(rand(@FOR));
					#$FOR[$i] = $grey;
					$FOR[$i] =$FOR[$random];
					my $readName = $print{$clust}{FOR}{name}[$random];
					$print{$clust}{FOR}{name}[$i] = $readName;
				}
			}

