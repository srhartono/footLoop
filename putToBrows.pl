#!/usr/bin/perl

use warnings; use strict; use Getopt::Std;
use vars qw($opt_f $opt_b $opt_o $opt_p $opt_t);
getopts("f:b:o:p:t:");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/jeep/lib';
   push(@INC, $libPath);
}
use myFootLib; use FAlite;

my $usage = "Usage: -f <BED file name> -b <Buffer length(or zero if none)> -o <Output file name> -p <Conversion threshold used> -t <track name>\n";

die $usage if not defined ($opt_f) or not defined ($opt_b) or not defined ($opt_o) or not defined ($opt_p) or not defined ($opt_t) or not -e ($opt_f);

my $output = $opt_o; $output =~ s/\/+/\//g;
my @folders = split("/", $output);
$output = "";
print "Output0 = $output\n";
for (my $i = 0; $i < @folders-1; $i++) {
	$output .= $folders[$i] . "/";
#	$output =~ s/^\/+/\// if $output =~ /^\//;
	print localtime() . ": Creating Directory: $output\n" if not -d $output;
	mkdir $output if not -d $output;
}
#die "Output folder $output does not exist somehow (permission issue?)\n" unless -d $output;
$output .= $folders[@folders-1];
my $GTFPos = $output . "_pos.gtf";
open(my $GTFPosfile, ">", $GTFPos) or die "Cannot write to $GTFPos: $!\n";
my $pname = $opt_t . "Pos";
my $nname = $opt_t . "Neg";
print $GTFPosfile "track name=$pname color=255,0,0\n";

my $GTFNeg = $output . "_neg.gtf";
open(my $GTFNegfile, ">", $GTFNeg) or die "Cannot write to $GTFNeg: $!\n";
print $GTFNegfile "track name=$nname color=0,0,255\n";

my $bed = $opt_f;
open(my $bedFile, "<", $bed) or die "No BED file found.\n";

my $i;
my $start;
my $end;

while(my $bedLine = <$bedFile>)
{
  my @bedLine = split ("\t", $bedLine);
  chomp $bedLine[3];
  my $Posreads = $bedLine[3] . "_Pos" . $opt_p . ".txt";
  my $Negreads = $bedLine[3] . "_Neg" . $opt_p . ".txt";
  my $z = 1;
  my $readsFile1 = $Posreads . "Clust";
  if(open(my $readsFile, "<", $Posreads))
  {
    open(my $rF1, ">", $readsFile1) or die "cannot open file \n";
    my $Rscript = "Poscluster.R";
    open(my $out, ">", $Rscript);
    print $out "
    library(\"GMD\")
    df = read.table(\"./$Posreads\", sep=\"\t\", row.names=1)
    df2 = df
    df2[df2 != 9] = 0
    df\$sum = apply(df2,1,sum)
    df = df[order(-df\$sum),]
    dimz = dim(df)[1]
    df = df[1:dimz,]
    df = subset(df,select=-sum)
    df2 = df
    df2[df2 !=9] = 0
    h = heatmap.3(
            df2,
            dendrogram=\"row\",
            Rowv=TRUE, Colv=FALSE,
            labRow=FALSE, labCol=FALSE
              )
    id = h\$rowInd
    id = rev(id)
    df = df[id,]
    dev.off()
    write.table(df, file = \"./$readsFile1\", sep = \"\t\")
    ";
    close $out;
    system("R --vanilla --no-save < $Rscript");
    open(my $readsFile3, "<", $readsFile1);
    my $trash = <$readsFile3>;
    while(my $readsLine = <$readsFile3>)
    {
      my @readsLine = split ("\t", $readsLine);
      $readsLine[0] =~ s/SEQ_//;
      $readsLine[0] =~ s/"//g;
      print $GTFPosfile "$bedLine[0]\tSOURCE\texon\t$bedLine[1]\t" , $bedLine[1]+1 , "\t0\t+\t0\tgene_id \"$readsLine[0]$z\"; transcript_id \"$readsLine[0]$z\"\n"; 
      print $GTFPosfile "$bedLine[0]\tSOURCE\texon\t" , $bedLine[2]-1 , "\t$bedLine[2]\t0\t+\t0\tgene_id \"$readsLine[0]$z\"; transcript_id \"$readsLine[0]$z\"\n";
      $i = 1;
      for($i; $i<@readsLine; $i++)
      {
        if($readsLine[$i] == 9 || $readsLine[$i] == 5)
        {
          $start = $i + $bedLine[1] - 1 - $opt_b + 1;
          $i++;
          while($readsLine[$i] == 9 || $readsLine[$i] == 5)
          {
            $i++;
          }
          $end = $i + $bedLine[1] - $opt_b - 1;
          print $GTFPosfile "$bedLine[0]\tSOURCE\texon\t$start\t$end\t0\t+\t0\tgene_id \"$readsLine[0]$z\"; transcript_id \"$readsLine[0]$z\"\n";
        }
      }
      $z++;
    }
  }
  else 
  {
    print "No positive reads file for $bedLine[3] detected. Omitting positive $bedLine[3] track from Genome Broswer.\n";
  }
  my $y = 1;
  my $readsFile1N = $Negreads . "Clust";
  if(open(my $readsFileN, "<", $Negreads))
  {
    open(my $rF1N, ">", $readsFile1N);
    my $RscriptN = "Negcluster.R";
    open(my $outN, ">", $RscriptN);
    print $outN "
    library(\"GMD\")
    df = read.table(\"./$Negreads\", sep=\"\t\", row.names=1)
    df2 = df
    df2[df2 != 9] = 0
    df\$sum = apply(df2,1,sum)
    df = df[order(-df\$sum),]
    dimz = dim(df)[1]
    df = df[1:dimz,]
    df = subset(df,select=-sum)
    df2 = df
    df2[df2 !=9] = 0
    h = heatmap.3(
            df2,
            dendrogram=\"row\",
            Rowv=TRUE, Colv=FALSE,
            labRow=FALSE, labCol=FALSE
              )
    id = h\$rowInd
    id = rev(id)
    df = df[id,]
    dev.off()
    write.table(df, file = \"./$readsFile1N\", sep = \"\t\")
    ";
    close $outN;
    system("R --vanilla --no-save < $RscriptN");
    open(my $readsFile3N, "<", $readsFile1N);
    my $trash = <$readsFile3N>;
    while(my $readsLineN = <$readsFile3N>)
    {
      my @readsLineN = split ("\t", $readsLineN);
      $readsLineN[0] =~ s/SEQ_//;
      $readsLineN[0] =~ s/"//g;
      print $GTFNegfile "$bedLine[0]\tSOURCE\texon\t$bedLine[1]\t" , $bedLine[1]+1 , "\t0\t+\t0\tgene_id \"$readsLineN[0]$y\"; transcript_id \"$readsLineN[0]$y\"\n";
      print $GTFNegfile "$bedLine[0]\tSOURCE\texon\t" , $bedLine[2]-1 , "\t$bedLine[2]\t0\t+\t0\tgene_id \"$readsLineN[0]$y\"; transcript_id \"$readsLineN[0]$y\"\n";
      $i = 1;
      for($i; $i<@readsLineN; $i++)
      {
        if($readsLineN[$i] == 9 || $readsLineN[$i] == 5)
        {
          $start = $i + $bedLine[1] - 1 - $opt_b + 1;
          $i++;
          while($readsLineN[$i] == 9 || $readsLineN[$i] == 5)
          {
            $i++;
          }
          $end = $i + $bedLine[1] - $opt_b - 1;
          print $GTFNegfile "$bedLine[0]\tSOURCE\texon\t$start\t$end\t0\t+\t0\tgene_id \"$readsLineN[0]$y\"; transcript_id \"$readsLineN[0]$y\"\n";
          $y++;
        }
      }
      $y++;
    }
  }
  else
  {
    print "No negative reads file for $bedLine[3] detected. Omitting negative $bedLine[3] track from Genome Broswer.\n";
  }
system("rm $readsFile1");
system("rm $readsFile1N");
}
system("rm Rplots.pdf");
system("rm Poscluster.R");
system("rm Negcluster.R");
#system("/bin/mv $GTFPos $opt_o");
#system("/bin/mv $GTFNeg $opt_o");
