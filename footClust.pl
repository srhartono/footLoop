#!/usr/bin/perl

use strict; use warnings; use mitochy; use Getopt::Std; use FAlite;
use vars qw($opt_v $opt_d);
getopts("vd:");


foreach my $input1 (sort @folder) {
my ($input1) = @;
die "\nusage: $YW$0$N $CY<footPeak bed file>$N\n\n" unless @ARGV == 1;
my $dist = defined $opt_d ? $opt_d : 50;
print "Dist = $dist\n";

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");
# get gene name
my ($gene) = $fileName1 =~ /^(.+)\_(Pos|Neg|Unk)/;
die "Cannot parse gene name from file=$LCY$input1$N\n" unless defined $gene;
$gene = uc($gene);

#get fasta file
my ($faFile) = <$folder1/../*.fa>;
my %data; my $linecount = 0;
my $total_read_all = 0; my $total_original_read = 0;
open (my $in1, "sort -k1,1 -k2,2n -k3,3n $input1|") or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /^#/;
	$linecount ++;
	my ($read, $beg, $end) = split("\t", $line);
	my ($num) = $read =~ /^.+\/(\d+)\/ccs/; 
	die "Read must be in this format: <anything>/<hole number>/ccs\n\n$read\n\n" if not defined $num;
	my $check = 0;
	$total_original_read ++;
	if (defined $data{$num}) {
		for (my $i = 0; $i < @{$data{$num}}; $i++) {
			my $beg2 = $data{$num}[$i][0];
			my $end2 = $data{$num}[$i][1];
			if ($beg < $end2 + $dist) {
				$data{$num}[$i][1] = $end;
				$check = 1;
				#print "$linecount: $num: MERGED $beg2-$end2 with $beg-$end into: $beg2-$end\n";
				last;
			}
		}
	}
	if ($check == 0) {
		push(@{$data{$num}}, [$beg,$end,$read,$gene]);
		$total_read_all ++;
#		print "$linecount: $num: $beg-$end\n";
	}
#	die if $linecount == 25;
}
close $in1;

my $total_read = (keys %data);
print "Unique read = $LGN$total_read$N, total read used = $LGN$total_read_all$N, original all read = $LGN$total_original_read$N\n";

open (my $out1, ">", "$input1.temp") or die "Cannot write to $fileName1.temp: $!\n";
print $out1 "id\tbeg\tend\n";
foreach my $num (sort keys %data) {
	for (my $i = 0; $i < @{$data{$num}}; $i++) {
		my $beg = $data{$num}[$i][0];
		my $end = $data{$num}[$i][1];
		my $read = $data{$num}[$i][2];
		print $out1 "$num.$i\t$beg\t$end\n";
	}
}
close $out1;

my $Rscript = "

set.seed(420)
library(ggplot2)
df = read.table(\"$input1.temp\",row.names=1,header=T,sep=\"\\t\")
dm = kmeans(df,5,nstart=20)
df\$cluster = dm\$cluster
df = df[order(df\$cluster, df[,1], df[,2]),]
df\$y = seq(1,dim(df)[1])
df\$ymax = seq(2,dim(df)[1]+1)
colnames(df) = c(\"x\",\"xmax\",\"clust\",\"y\",\"ymax\")
df2 = as.data.frame(aggregate(df[,c(1,2)],by=list(df\$clust),function(x)mean(x,trim=0.05)))
colnames(df2) = c(\"clust\",\"x2\",\"y2\")
df2 = df2[order(df2\$x2, df2\$y2),]
df2\$clust2 = seq(1,dim(df2)[1])
df\$id = rownames(df)
df = as.data.frame(merge(df,df2[,c(1,4)]))
df\$clust = df\$clust2
df = df[order(df\$clust, df\$x, df\$xmax),]
df\$y = seq(1,dim(df)[1])
df\$ymax = seq(2,dim(df)[1]+1)
df[,c(1,6)] = df[,c(6,1)]
colnames(df)[c(1,6)] = colnames(df)[c(6,1)]
df = df[,-7]
png(\"$input1.clust.png\",height=(dim(df)[1])*5,width=max(df\$xmax)+10)
ggplot(df, aes(x,y)) + geom_rect(aes(xmin=x,ymin=y,xmax=xmax,ymax=ymax,fill=as.factor(clust))) + theme_bw() + 
theme(panel.grid=element_blank()) + coord_fixed(ratio=10)
dev.off()
write.table(df,\"$input1.clust\",quote=F,row.names=F,col.names=T,sep=\"\\t\")
";
open (my $outR, ">", "$input1.temp.R") or die;
print $outR $Rscript;
system("R --vanilla --no-save < $input1.temp.R > $input1.temp.R.log 2>&1");
close $outR;


# process clust
my %cl;
open (my $in2, "<", "$input1.clust") or die;
while (my $line = <$in2>) {
	chomp($line);
	my ($num, $beg, $end, $y, $y2, $clust) = split("\t", $line);
	my ($ind) = "";
	($num, $ind) = $num =~ /^(\d+)\.(\d+)$/ if $num =~ /^\d+\.\d+$/;
	die "Cannot parse num from line=$line\n" if not defined $num;
	my ($mid) = int(($end + $beg)/2+0.5);
	push(@{$cl{$clust}{beg}}, $beg);
	push(@{$cl{$clust}{mid}}, $mid);
	push(@{$cl{$clust}{end}}, $end);
}
close $in2;

# process clust: get average beg/end point
open (my $out2, ">", "$input1.clust.bed") or die;
foreach my $clust (sort {$a <=> $b} keys %cl) {
	my $beg = tmm(@{$cl{$clust}{beg}});
	my $mid = tmm(@{$cl{$clust}{mid}});
	my $end = tmm(@{$cl{$clust}{end}});
	print $out2 "
}










__END__
close $out1;

#df = read.table("CALM3_Pos_20_0.65_CH.PEAK.local.bed.temp",row.names=1,header=F,sep="\t")

set.seed(420)
library(ggplot2)
df = read.table("CALM3_Pos_20_0.65_CH.PEAK.local.bed.temp",row.names=1,header=T,sep="\t")
dm = kmeans(df,6,nstart=20)
df$cluster = dm$cluster
df = df[order(df$cluster, df[,1], df[,2]),]
df$y = seq(1,dim(df)[1])
df$ymax = seq(2,dim(df)[1]+1)
colnames(df) = c("x","xmax","clust","y","ymax")
df2 = as.data.frame(aggregate(df[,c(1,2)],by=list(df$clust),function(x)mean(x,trim=0.05)))
colnames(df2) = c("clust","x2","y2")
df2 = df2[order(df2$x2, df2$y2),]
df2$clust2 = seq(1,dim(df2)[1])
df$id = rownames(df)
df = as.data.frame(merge(df,df2[,c(1,4)]))
df$clust = df$clust2
df = df[order(df$clust, df$x, df$xmax),]
df$y = seq(1,dim(df)[1])
df$ymax = seq(2,dim(df)[1]+1)
df[,c(1,6)] = df[,c(6,1)]
colnames(df)[c(1,6)] = colnames(df)[c(6,1)]
df = df[,-7]
png("CALM3_Pos_20_0.65_CH.PEAK.local.bed.clust.png",height=(dim(df)[1])*5,width=max(df$xmax)+10)
ggplot(df, aes(x,y)) + geom_rect(aes(xmin=x,ymin=y,xmax=xmax,ymax=ymax,fill=as.factor(clust))) + theme_bw() + 
coord_fixed(ratio=10)
dev.off()

