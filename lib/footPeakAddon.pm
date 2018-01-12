package footPeakAddon;

use strict; use warnings; use Getopt::Std; use FAlite; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_x $opt_R $opt_c $opt_t);

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/footLoop/lib';
   push(@INC, $libPath);
}
use myFootLib; use FAlite;
my $homedir = $ENV{"HOME"};
my $footLoopDir = dirname(dirname abs_path $0) . "/footLoop";

sub main {
	my ($input1, $faFile, $mygene, $minDis, $SEQ) = @_;
	my @coor = split("\t", $SEQ->{$mygene});
	my (%pk, %Rscripts);
	print date() . "\nusage: $YW$0$N [-c to use cpg] $CY<CALM3_Pos_20_0.65_CG.PEAK>$N $CY<location with lots of C>$N\n\n" and exit 1 unless @_ == 5;
	print date() . "Input cannot be directry!\n" and exit 1 if -d $input1;
	($input1) = getFullpath($input1);
	my ($folder, $fileName) = getFilename($input1, "folderfull");
	my ($gene, $strand, $window, $thres, $type, $isPeak) = $fileName =~ /^(\w+)_(Unk|Pos|Neg)_(\d+)_(\d+\.?\d*)_(\w+)\.(PEAK|NOPK)$/;
	open (my $outLog, ">", "$folder/footLoop_addition_logFile.txt") or die;

	if (defined $mygene) {
		LOG($outLog, date() . "footPeak.pl mygene=$mygene, input genez=$gene are not the same!\n") and exit 1 if uc($mygene) ne uc($gene);;
	}
	$mygene = $gene if not defined $mygene;
	
	my $lotsOfC = find_lots_of_C($faFile, $mygene, $outLog) if defined $faFile;
	my %bad;
	if (defined $lotsOfC) {
		my @lotsOfC = split(",", $lotsOfC);
		foreach my $coor (@lotsOfC) {
			my ($nuc, $beg, $end) = split(";", $coor);
			my $strands = $nuc eq "C" ? 0 : $nuc eq "G" ? 16 : 255;
			$bad{$strands}{$beg} = $end-1;
			LOG($outLog, date() . "$strands: coor=$coor, nuc=$nuc, beg=$beg, end=$end, beg-end=$beg-$bad{$strands}{$beg}\n");
		}
	}

	LOG($outLog, date() . "$input1; Undefined mygene=$mygene=, strand=$strand=, window=$window=, thres=$thres=, type=$type=, isPeak=$isPeak=\n") and exit 1 if not defined $isPeak or not defined $window;
	my %total; 
	$total{Pos}{peak} = 0; $total{Pos}{nopeak} = 0; $total{Pos}{total} = 0;
	$total{Neg}{peak} = 0; $total{Neg}{nopeak} = 0; $total{Neg}{total} = 0;
	$total{Unk}{peak} = 0; $total{Unk}{nopeak} = 0; $total{Unk}{total} = 0;
	my $log2 = "";
	LOG($outLog, date . "\n\nFolder $YW$folder$N: Processing files related to $LCY$input1$N\n");
	open (my $outLGENE, ">", "$folder/.0_RESULTS\_$mygene\_$strand\_$window\_$thres.TXT") if not defined $opt_x;
	open (my $outLEXTRA, ">", "$folder/.1_RESULTS_EXTRA\_$mygene\_$strand\_$window\_$thres.TXT") if not defined $opt_x;
	my %files;
	for (my $h = 0; $h < 4; $h++) {
		$type = $h % 4 == 0 ? 'CH' : $h % 4 == 1 ? 'CG' : $h % 4 == 2 ? 'GH' : 'GC';
		my $peakFile   = "$folder/$mygene\_$strand\_$window\_$thres\_$type.PEAK";
		my $nopkFile   = "$folder/$mygene\_$strand\_$window\_$thres\_$type.NOPK";
		$files{$peakFile} = 1;
		LOG($outLog, date . "h=$LGN$h\t$YW$peakFile\t$LCY$nopkFile\n$N");
	
		my ($folder1, $peakfileName) = getFilename($peakFile, "folder");
		my ($folder2, $nopeakfileName) = getFilename($nopkFile, "folder");
		
		my %data; my $totalnopeak = 0;
		my $linecount = 0;
		my $totalline = 0;
		if (-e $nopkFile) {
			($totalline) = `wc -l $nopkFile` =~ /^(\d+)/;
			open (my $in2, "<", $nopkFile) or LOG($outLog, date() . "Cannot read from $nopkFile: $!\n") and exit 1 ;
			LOG($outLog, date . "\tProcessing $LPR$nopkFile$N ($LGN$totalline$N lines)\n");
			while (my $line = <$in2>) {
				chomp($line);
				$linecount ++;
				next if $linecount == 1;
				LOG($outLog, date . "\tDone $totalnopeak / $totalline\n") if $totalnopeak % 500 == 0;
				my ($name, $val, $totalPeak, $peaks) = parse_peak($line, \%bad, $minDis, $outLog);
				$val = "$name\t" . join("\t", @{$val});
				push(@{$data{peak}}, $val) if $totalPeak > 0;
				push(@{$data{nopeak}}, $val) if $totalPeak == 0;
				$totalnopeak ++;
				$pk{$peakFile}{$name} = $peaks if defined $peaks;
			}
			close $in2;
		}
	
		my $peakCount = defined $data{peak} ? @{$data{peak}} : 0;
		my $nopeakCount = defined $data{nopeak} ? @{$data{nopeak}} : 0;
		my $nopeakPrint ="$folder2\t$nopeakfileName\t$peakCount\t$nopeakCount\t$totalnopeak\t$totalline";
		LOG($outLog, date() . "$nopeakPrint\n");
		$total{$type}{peak} += $peakCount;
		$total{$type}{nopeak} += $nopeakCount;
		$total{$type}{total} += $totalnopeak;
		
		my $totalpeak = 0;
		if (-e $peakFile) {
			open (my $in1, "<", $peakFile) or LOG($outLog, date() . "Cannot read from $peakFile: $!\n") and exit 1 ;
			($totalline) = `wc -l $peakFile` =~ /^(\d+)/;
			LOG($outLog, date . "\tProcessing $LPR$peakFile$N ($LGN$totalline$N lines)\n");
			$linecount = 0;
			while (my $line = <$in1>) {
				chomp($line);
				$linecount ++;
				next if $linecount == 1;
				LOG($outLog, date . "\tDone $totalpeak / $totalline\n") if $totalpeak % 500 == 0;
				my ($name, $val, $totalPeak, $peaks) = parse_peak($line, \%bad, $minDis, $outLog);
				$val = "$name\t" . join("\t", @{$val});
				push(@{$data{peak}}, $val) if $totalPeak > 0;
				push(@{$data{nopeak}}, $val) if $totalPeak == 0;
				$totalpeak ++;
				$pk{$peakFile}{$name} = $peaks if defined $peaks;
			}
			close $in1;
		}
		$peakCount = defined $data{peak} ? @{$data{peak}} - $peakCount : 0;
		$nopeakCount = defined $data{nopeak} ? @{$data{nopeak}} - $nopeakCount : 0;
		my $peakPrint ="$folder1\t$peakfileName\t$peakCount\t$nopeakCount\t$totalpeak\t$totalline";
		LOG($outLog, date() . "$peakPrint\n");
		$total{$type}{peak} += $peakCount;
		$total{$type}{nopeak} += $nopeakCount;
		$total{$type}{total} += $totalpeak;
	
	
	
		open (my $out1, ">", "$folder1/$peakfileName.out") or LOG($outLog, date() . "Cannot write to $peakfileName.out: $!\n") and exit 1 ;
		open (my $out2, ">", "$folder1/$nopeakfileName.out") or LOG($outLog, date() . "Cannot write to $nopeakfileName.out: $!\n") and exit 1 ;
		if (defined $data{peak}) {
			foreach my $val (sort @{$data{peak}}) {
				print $out1 "$val\n";
			}
		}
		if (defined $data{nopeak}) {
			foreach my $val (sort @{$data{nopeak}}) {
				print $out2 "$val\n";
			}
		}
		
		close $out1;
		close $out2;
	
		$log2 .= "\#Folder\tFile\tPeak\tNoPeak\tTotalRead\tTotalLineInFile\n" if $h == 0;
		LOG($outLog, date . "$peakPrint\n$nopeakPrint\n") if defined $opt_x;
		$log2 .= "$peakPrint\n$nopeakPrint\n" if not defined $opt_x;
	
		mkdir "$folder1/remove" if not -d "$folder1/remove/";
		
		my $peakFileBackup = "$folder1/remove/$peakfileName";
		my $peakFileTemp   = $peakFileBackup;
		my $count = 0;
		while (-e $peakFileTemp) {
			$peakFileTemp = $peakFileBackup . $count;
			$count ++;
		}
		if (-e $peakFile) {
	#		LOG($outLog, date . "\tmv $peakFile $peakFileTemp\n");
	#		system("/bin/mv $peakFile $peakFileTemp") if not defined $opt_x;
	#		LOG($outLog, date . "\tmv $folder1/$peakfileName.out $folder1/$peakfileName\n");
	#		system("mv $folder1/$peakfileName.out $folder1/$peakfileName") if not defined $opt_x;
		}
		
		my $nopkFileBackup = "$folder1/remove/$nopeakfileName";
		my $nopkFileTemp   = $nopkFileBackup;
		$count = 0;
		while (-e $nopkFileTemp) {
			$nopkFileTemp = $nopkFileBackup . $count;
			$count ++;
		}
		if (-e $nopkFile) {
	#		LOG($outLog, date . "\t/bin/mv $nopkFile $nopkFileTemp\n");
	#		system("/bin/mv $nopkFile $nopkFileTemp") if not defined $opt_x;
	#		LOG($outLog, date . "\tmv $folder1/$nopeakfileName.out $folder1/$nopeakfileName\n");
	#		system("mv $folder1/$nopeakfileName.out $folder1/$nopeakfileName") if not defined $opt_x;
		}
	}
	
	foreach my $file (sort keys %pk) {
		next if not defined $pk{$file};
		next if defined $pk{$file} and keys %{$pk{$file}} == 0;
		open (my $outPEAKS, ">", "$file.PEAKS") or die;
		open (my $outRPEAKS, ">", "$file.RPEAKS") or die;
		my $currtype = $file =~ /_CH/ ? "CH" : $file =~ /_CG/ ? "CG" : $file =~ /_GH/ ? "GH" : $file =~ /_GC/ ? "GC" : "TYPE_UNK";
		foreach my $name (sort keys %{$pk{$file}}) {
			next if not defined $pk{$file}{$name};
			my @arr = @{$pk{$file}{$name}};
			foreach my $peakz (sort @{$pk{$file}{$name}}) {
				my ($beg, $end) = split("-", $peakz);
				my ($chr0, $beg0, $end0, $name0, $val0, $strand0) = @coor;
				$end0 = $end + $beg0;
				$beg0 = $beg + $beg0;
				my ($junk, $readname) = $name =~ /^(\w+)\.(.+)$/;
				print $outPEAKS  "$chr0\t$beg0\t$end0\t$name\t0\t$strand0\t$file\n";
				print $outRPEAKS "$name\t$beg\t$end\n";
			}
		}
		close $outPEAKS;
		close $outRPEAKS;
	}
	
	if (not defined $opt_x) {
		my @typez = qw(CH CG GH GC);
		foreach my $typez (@typez[0..@typez-1]) {
			$total{$typez}{total} = 0 if not defined $total{$typez}{total};
			$total{$typez}{peak} = 0 if not defined $total{$typez}{peak};
			$total{$typez}{nopeak} = 0 if not defined $total{$typez}{nopeak};
			$total{$typez}{peak} = $total{$typez}{total} == 0 ? 0 : int(1000 * $total{$typez}{peak} / $total{$typez}{total}+0.5)/10;
			$total{$typez}{nopeak} = $total{$typez}{total} == 0 ? 0 : int(1000 * $total{$typez}{nopeak} / $total{$typez}{total}+0.5)/10;
			my @folder = split("/", $folder);
			my $foldershort = $folder[@folder-1];
			   $foldershort = $folder[@folder-2] if not defined ($foldershort) or (defined $foldershort and $foldershort =~ /^[\s]*$/);
			print $outLGENE "#folder\tfilename\tGene\tStrand\ttotal\tpeak.perc\n" if $typez eq "CH";
			print $outLGENE "$foldershort\t$fileName\t$mygene\t$typez\ttotal=$total{$typez}{total}\tpeak=$total{$typez}{peak}\n";
			print $outLEXTRA "$log2";
		}
		close $outLGENE;
		close $outLEXTRA;
	}
	system("cat $folder/.0_RESULTS\_$mygene\_$strand\_$window\_$thres.TXT") if not defined $opt_x;
#	LOG($outLog, date . "\tcd $folder && run_Rscript.pl *MakeHeatmap.R\n");
#	system("cd $folder && run_Rscript.pl *MakeHeatmap.R") if not defined $opt_x and defined $opt_R;

	my ($sampleName) = $folder =~ /\/\d+_(m\d+_\d+)_\d+_\w+/;
	if (not defined $sampleName) {
		$sampleName = $folder;
		$sampleName =~ s/\//_/g;
	}
	foreach my $file (sort keys %files) {
		next if not defined $files{$file};
		my $peakFile = $file . ".out";
		LOG($outLog, date() . "Doing $peakFile\n");
		my $nopkFile = $file . ".out"; $nopkFile =~ s/\.PEAK.out/.NOPK.out/;
		my $totpeak = linecount($peakFile);
		my $totnopk = linecount($nopkFile);
		my $bedFile = $peakFile . ".RPEAKS"; $bedFile =~ s/.out.RPEAKS$/.RPEAKS/;
		for (my $p = 0; $p < 2; $p ++) {
			my $currFile = $p == 0 ? $peakFile : $nopkFile;
			my ($currFolder, $currFilename) = getFilename($currFile, "folderfull");
			my $pngout = "$currFolder/$sampleName\_$currFilename.png";
			my $pdfout = "$currFolder/$sampleName\_$currFilename.pdf";
			next if linecount($currFile) <= 1;
			my $Rscript = "library(ggplot2);library(reshape2)\nlibrary(grid)\nlibrary(gridExtra)\n";
			$Rscript .= "
				df = read.table(\"$currFile\",skip=1,sep=\"\\t\")
				colnames(df) = c(\"V1\",seq(1,dim(df)[2]-1))
				h = hclust(dist(df[,-1]))
				df = df[h\$order,]
				df\$y = seq(1,dim(df)[1])
			";
			if (-e $bedFile and linecount($bedFile) > 0) {
				$Rscript .= "bed = read.table(\"$bedFile\",sep=\"\\t\");\n";# bed = bed[,c(4,2,3)];colnames(bed) = c(\"V1\",\"V2\",\"V3\");\n";
				$Rscript .= "bed = merge(subset(df,select=c(\"V1\",\"y\")),bed,by=\"V1\")\n";
			}
			$Rscript .= "
				dm = melt(df,id.vars=c(\"V1\",\"y\"))
				dm\$variable = as.numeric(as.character(dm\$variable))
				p = ggplot(dm,aes(variable,y)) +  
					geom_tile(aes(fill=as.factor(value))) + 
					theme_bw() + theme(legend.position=\"none\") + 
#					coord_fixed() +
					scale_fill_manual(values=c(\"0\"=\"grey\",\"1\"=\"white\",\"4\"=\"cornsilk\",\"5\"=\"cornsilk\",
								  					   \"6\"=\"green4\",\"7\"=\"seagreen4\",\"8\"=\"red4\",\"9\"=\"maroon4\")) +
					scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
               theme(line = element_blank(),
                     axis.text = element_blank(),
                     axis.title = element_blank()
                    ) + ggtitle(paste(\"(peak=\",$totpeak,\"; nopk=\",$totnopk,\")\",sep=\"\"))
				";
			if (-e $bedFile and linecount($bedFile) > 0) {
				$Rscript .= "
				if (length(bed) != 0 & dim(bed)[1] > 0) {
					bed\$variable=1
					p = p + geom_rect(data=bed,aes(xmin=V2,xmax=V3,ymin=y-0.5,ymax=y+0.5),
											size=0.5,fill=rgb(0,0,0,alpha=0),color=rgb(1,0,0,alpha=0.25))
				}
				";
			}
				my $change = $currFile =~ /PEAK.out$/ ? "df2[df2 < 8] = 0; df2[df2 >= 8] = 1" : "df2[df2 < 6] = 0; df2[df2 >= 6] = 1";

				$Rscript .= "

            df2 = subset(df,select=c(-V1,-y));
				$change
            df2 = data.frame(x=seq(1,dim(df2)[2]), y=apply(df2,2,mean));
            if (dim(df2[df2\$y > 0,])[1] > 15) {
               df2 = df2[df2\$y > 0,]
               df2\$x = as.numeric(as.character(df2\$x));
               df2\$y = as.numeric(as.character(df2\$y));
               df2\$x2 = df2\$x
               df2\$y2 = df2\$y
               for (i in 1:(dim(df2)[1]-10)) {
                  a = df2[df2\$x >= df2[i,]\$x & df2\$x <= df2[i+10-1,]\$x,]
                  if (length(a) != 0 & dim(a)[1] != 0) {
                        df2[i,]\$y2 = mean(a\$y)
                        df2[i,]\$x2 = mean(a\$x)
                  }
               }
#              df2 = rbind(df2[seq(1,5),],df2); df2\$x2[1:5] = df2\$x[1:5]
               mins = seq(1,as.integer(df2[1,]\$x2)-1,10)
               df2 = rbind(data.frame(x=mins,y=0,x2=mins,y2=0),df2)
            } else {
               df2 = data.frame(x=seq(1,dim(df)[2]), y=0, x2=seq(1,dim(df)[2]), y2=0);
            }
            p2 = ggplot(df2,aes(x2,y2)) + geom_point(aes(x=x,y=y),size=1) + geom_line(color=rgb(1,0,0,alpha=1)) + theme_bw()+
               scale_x_continuous(expand = c(0,0)) +
               scale_y_continuous(expand = c(0,0)) +
               theme(line = element_blank(),axis.text = element_blank(),axis.title = element_blank()) +
               annotate(geom='text',x=10,y=1,label=\"- 100 \%\",size=5,hjust=0) +
               annotate(geom='text',x=10,y=0,label=\"- 0   \%\",size=5,hjust=0) +
               coord_cartesian(ylim=c(-0.05,1.05))

				png(\"$pngout\",width=dim(df)[2]*2,height=dim(df)[1]*16)
            grid.arrange(p,p2,ncol=1,nrow=2);
				dev.off()
			";

			open (my $outRscript, ">", "$currFile.R") or die;
			print $outRscript $Rscript;
			$Rscripts{"$currFile.R"} = 1;
			close $outRscript;
		}
	}
	foreach my $outRscript (sort keys %Rscripts) {
		LOG($outLog, date() . "run_Rscript.pl $outRscript > $outRscript.LOG 2>&1\n");
		system("run_Rscript.pl $outRscript > $outRscript.LOG 2>&1") == 0 or LOG($outLog, date() . "Failed to run_Rscript.pl $outRscript: $!\n");
	}
}
###############
# Subroutines #
###############


sub make_heatmap {
	

}





sub parse_peak {
	my ($ARG, $bad, $minDis, $outLog) = @_;
	my ($name, $isPeak, $mygene, $type, $strand, @val) = split("\t", $ARG);
	my %bad; %bad = %{$bad} if defined $bad;
	my $name_want = "CALM3.m160130_030742_42145_c100934342550000001823210305251633_s1_p0/9477/ccs";
	shift(@val) if $val[0] eq "";
#	for (my $i = 0; $i < @val; $i++) {
#		if ($val[$i] !~ /^[456789]$/) {print "."} else {print "$val[$i]";}
#		LOG($outLog, date() . "\n") if $i != 0 and ($i+1) % 100 == 0;
#	}
#	LOG($outLog, date() . "\n");
#	0001234000
#	0123456789
#	len=10, e1=len(e1), e2=10-len(e2)
	my $peaks;
	my %peak; $peak{curr} = 0; #my $edge = 0; my $edge2 = 0; my $zero = 0; my $edge1 = 0;
	my $Length = @val; 
	my $print = "Total length = $Length\n";
	my ($edge1) = join("", @val) =~ /^(0+)[\.1-9A-Za-z]/;
	$edge1 = defined $edge1 ? length($edge1) : 0;
	my ($edge2) = join("", @val) =~ /[\.1-9A-Za-z](0+)$/;
	$edge2 = defined $edge2 ? @val-length($edge2) : @val;
	for (my $i = 0; $i < @val; $i++) {
		my $val = $val[$i];
		if ($i % 100 == 0) {$print .= "\n$YW" . $i . "$N:\t";}
		if ($val[$i] =~ /[89]/) {
			$peak{beg} = $i if $peak{curr} == 0;
			$print .= "${LPR}$val[$i]$N" if $peak{curr} == 0;
			$print .= "${LRD}$val[$i]$N" if $peak{curr} == 1;
			$peak{curr} = 1;
		}
		elsif ($val[$i] =~ /[23]/) {
			$peak{end} = $i+1;
			push(@{$peak{peak}}, "$peak{beg}-$peak{end}");
			undef $peak{beg}; undef $peak{end};
			$peak{curr} = 0;
			$val[$i] =~ tr/23/89/;
			$print .= "${LPR}$val$N";
		}
		else {
			$print .= "EDGE1" if $i == $edge1;
			$print .= "${LGN}$val[$i]$N" if $val =~ /^[46]$/;
			$print .= "${LGN}$val[$i]$N" if $val =~ /^[57]$/;
			$print .= "." if $val[$i] eq 1;
			$print .= "x" if $val[$i] eq 0;# and $i < $edge1;
			$print .= "EDGE2" if $i == $edge2 - 1;
		}
	}
	my (%nopeak, @peak);
	$strand = $type =~ /^C/ ? 0 : $type =~ /^G/ ? 16 : 255;
	my %peak2;
	$print .= "\n";
#	LOG($outLog, date() . "\nDoing $YW$name$N\n" if $name eq $name_want;#"SEQ_76074" or $name eq "SEQ_34096" or $name eq "SEQ_62746");
	if (defined $peak{peak}) {
		foreach my $peak (sort @{$peak{peak}}) {
			my ($beg, $end) = split("-", $peak);
#			LOG($outLog, date() . "$name: $beg to $end\n" if $name eq 77011 or $name eq "$name_want");
			my $checkBad = 0;
#			if (not defined $bad{$strand}) {
#				push(@peak, "$beg-$end");
#				push(@{$peak2{peak}}, $peak);
#			}
#			next if not defined $bad{$strand};
			foreach my $begBad (sort keys %{$bad{$strand}}) {
				my $endBad = $bad{$strand}{$begBad};
#			foreach my $begBad (sort keys %bad) {
#				my $endBad = $bad{$strand}{$begBad};
#				LOG($outLog, date() . "\t$beg-$end in begBad=$begBad to endBad=$endBad?\n" if $name eq "$name_want");
				if ($beg >= $begBad and $beg <= $endBad and $end >= $begBad and $end <= $endBad) {
#					for (my $m = $beg; $m <= $end; $m++) {
#						$nopeak{$m} = 1;
#					}
					LOG($outLog, date() . "\t\t$LGN YES$N peak=$beg-$end, bad=$begBad-$endBad\n") if $isPeak eq "PEAK";# if $name eq "$name_want");
					$checkBad = 1; last;
				}
			}
			if ($checkBad != 1) {
				next if not defined $bad{$strand};
				foreach my $begBad (sort keys %{$bad{$strand}}) {
					my $endBad = $bad{$strand}{$begBad};
#					LOG($outLog, date() . "$name: $beg-$end in begBad=$begBad to endBad=$endBad?\n" if $name eq "SEQ_76074" or $name eq "SEQ_34096" or $name eq "$name_want");
					#if (($beg >= $begBad and $beg <= $endBad) or ($end >= $begBad and $end <= $endBad)) {
						my @valz = @val;
						my ($goodC, $badC) = (0,0);
						for (my $m = $beg; $m < $end; $m++) {
							if ($m >= $begBad and $m <= $endBad) {
								$badC ++ if $valz[$m] =~ /[2389]/;
							}
							else {
								$goodC ++ if $valz[$m] =~ /[2389]/;
							}
						}
						if ($goodC < 5 and $badC >= 9) {
#							LOG($outLog, date() . "\t$LRD NO!$N beg=$beg, end=$end, begBad=$begBad, endBad=$endBad, badC = $badC, goodC = $goodC\n" if $name eq "SEQ_76074" or $name eq "SEQ_34096" or $name eq "SEQ_62746");
#							LOG($outLog, date() . "\t$YW$name$N $LRD NO!$N beg=$beg, end=$end, begBad=$begBad, endBad=$endBad, badC = $badC, goodC = $goodC\n");
							$checkBad = 1; last;
						}
#						else {
#							#LOG($outLog, date() . "\t$LGN OKAY!$N beg=$beg, end=$end, begBad=$begBad, endBad=$endBad, badC = $badC, goodC = $goodC\n" if $name eq "SEQ_76074" or $name eq "SEQ_34096" or $name eq "SEQ_62746");
#						}
#					LOG($outLog, date() . "\t\t$LGN YES$N\n" if $name eq "$name_want");
					#}
				}
			}
#			LOG($outLog, date() . "\t$name checkbad = $checkBad\n" if $name eq "SEQ_76074" or $name eq "SEQ_34096" or $name eq "SEQ_62746");
			
			if ($checkBad == 1) {
				$print .= "\tend=$end > 100 + edge1=$edge1 OR beg=$beg < edge2=$edge2-100; Peak Bad : $LRD$peak$N\n";
#				LOG($outLog, date() . "\tend=$end > 100 + edge1=$edge1 OR beg=$beg < edge2=$edge2-100; Peak bad : $LRD$peak$N\n" if $name eq "$name_want");
				for (my $j = $beg; $j <= $end; $j++) {
					$nopeak{$j} = 1;
				}
			}
			elsif ($end > 10 + $edge1 or $beg < $edge2 - 10) {
#				LOG($outLog, date() . "something wrong\n";# if $name eq "$name_want");
				$print .= "\tend=$end > 100 + edge1=$edge1 OR beg=$beg < edge2=$edge2-100; Peak Used: $LGN$peak$N\n";
				push(@peak, "$beg-$end");
				push(@{$peak2{peak}}, $peak);
			}
			else {
				$print .= "\tend=$end > 10 + edge1=$edge1 OR beg=$beg < edge2=$edge2-10; Peak Not : $LRD$peak$N\n";
				for (my $j = $beg; $j <= $end; $j++) {
					$nopeak{$j} = 1;
				}
			}
		}
	}
	my $totalpeak = scalar(@peak);
#	my @val2;
#	for (my $i = 0; $i < @val; $i++) {
#		my $val = $val[$i];
#		$val2[$i] = $val;
	#	if ($val =~ /^(8|9)$/ and defined $nopeak{$i}) { # Peak Converted CpG or CH
	#		#$val2[$i] = 7 if $val eq 9;
	#		#$val2[$i] = 6 if $val eq 8;
	#	}
	#}
	#die $print if $totalpeak > 1;
	$print .= "$name\t$totalpeak\n" if $isPeak eq "PEAK";
#	LOG($outLog, date() . "$print\n" if $isPeak eq "PEAK";# and $print =~ /; Peak Not/;# if $totalpeak == 1;# or $name eq "SEQ_100022") and exit 1 ;
#	exit 0 if $isPeak eq "PEAK";
	return ($name, \@val, $totalpeak, $peak2{peak});
}

sub find_lots_of_C {
	my ($seqFile, $mygene, $outLog) = @_;#, $geneIndex, $box) = @_; #$geneIndexesFa;
	my %seq;
	LOG($outLog, date . "\n${YW}2. Parsing in sequence for genes from sequence file $CY$seqFile$N\n");
	open(my $SEQIN, "<", $seqFile) or LOG($outLog, date() . "\n$LRD!!!$N\tFATAL ERROR: Could not open $CY$seqFile$N: $!") and exit 1 ;
	my $fasta = new FAlite($SEQIN);
	my %lotsOfC;

	while (my $entry = $fasta->nextEntry()) {
	   my $gene = uc($entry->def); $gene =~ s/^>//;
	   my $seqz = uc($entry->seq);
		next if $gene ne $mygene;
	   LOG($outLog, date . "\t\tgenez=$gene ($gene)\n");
	   
	
		my $minlen = 6;
	   my $seqz2 = $seqz;#join("", @{$seq{$gene}{seq}});
	   while ($seqz2 =~ /(C){$minlen,99}/g) {
	      my ($prev, $curr, $next) = ($`, $&, $');
	      my ($curr_C) = length($curr);
	      my ($next_C) = $next =~ /^(C+)[AGTN]*$/;
	      $next_C = defined $next_C ? length($next_C) : 0;
	      my ($beg_C) = defined $prev ? length($prev) : 0;
	      my ($end_C) = $curr_C + $next_C + $beg_C;
	      my $length = $curr_C + $next_C;
	      ($prev) = $prev =~ /^.*(\w{$minlen})$/ if length($prev) > $minlen; $prev = "NA" if not defined $prev;
	      ($next) = $next =~ /^(\w{$minlen}).*$/ if length($next) > $minlen; $next = "NA" if not defined $next;
	      LOG($outLog, date() . "$gene: $beg_C to $end_C ($length)\n\tPREV=$prev\n\tCURR=$curr\n\tNEXT=$next\n");
	      $lotsOfC{$gene} .= "C;$beg_C;$end_C,";
	   }
		$seqz2 = "";
	   $seqz2 =$seqz;# join("", @{$seq{$gene}{seq}});
	   while ($seqz2 =~ /(G){$minlen,99}/g) {
	      my ($prev, $curr, $next) = ($`, $&, $');
	      my ($curr_G) = length($curr);
	      my ($next_G) = $next =~ /^(G+)[ACTN]*$/;
	      $next_G = defined $next_G ? length($next_G) : 0;
	      my ($beg_G) = defined $prev ? length($prev) : 0;
	      my ($end_G) = $curr_G + $next_G + $beg_G;
	      my $length = $curr_G + $next_G;
	      ($prev) = $prev =~ /^.*(\w{$minlen})$/ if length($prev) > $minlen; $prev = "NA" if not defined $prev;
	      ($next) = $next =~ /^(\w{$minlen}).*$/ if length($next) > $minlen; $next = "NA" if not defined $next;
	      LOG($outLog, date() . "$gene: $beg_G to $end_G ($length)\n\tPREV=$prev\n\tCURR=$curr\n\tNEXT=$next\n");
	      $lotsOfC{$gene} .= "G;$beg_G;$end_G,";
	   }
		$lotsOfC{$gene} =~ s/,$//;
		return $lotsOfC{$gene} if defined $lotsOfC{$gene};
		return;
	}
}
#foreach my $gene (keys %lotsOfC) {
#   $gene = uc($gene);
#   $lotsOfC{$gene} =~ s/;$//;
#   LOG($outLog, date() . "$gene\t$lotsOfC{$gene}\n");
#   my $beg2 = $geneIndex{$gene};
#   foreach my $lines (@{$box->{$gene}}) {
#      LOG($outLog, date() . "GENEZ = $gene, lines = $lines\n");
#   }
#   LOG($outLog, date() . "genez=$gene,beg=$beg2\n");
#}
#push


1;
__END__
# 0 = not converted
# 1 = converted C
# 2 = A T or G (non C)
# 3 = Non converted CpG
# 4 = Converted CpG
# 5 = PEAK Converted CpG
# 6 = No data
# 9 = PEAK Converted C

   # For nucleotide
# 10 = Nucleotide A
# 11 = Nucleotide C
# 12 = Nucleotide T
# 13 = Nucleotide G






__END__
END1
#my (@inputs, $input1);
#if ($input1 =~ /.orig$/) {
#}
#else {
#	(@inputs) = <$folders/*Pos50.orig>;
#	LOG($outLog, date() . "Must have 1 input only! (" . scalar(@inputs) . "):\n" . join("\n-", @inputs) . "\n") and exit 1 if @inputs != 1 and not defined $opt_x;
#	$input1 = defined $opt_x ? $folders : $inputs[0];
#	LOG($outLog, date() . "INPUT1=$input1\n");
#}

