# footLoop Pipeline v1.6

# 0. SYNOPSIS

```
# 1. Install:
git clone https://github.com/srhartono/footLoop

# 2. Check if required software exists:
## Make sure you have bedtools (v2.25.0), bowtie2 (v2.2.6), bismark2 (v0.20.0), R (v3.6.1), samtools (v0.1.19)
## NEWER VERSIONS of these softwares MIGHT NOT WORK due to output file differences. If you get errors it's very likely due to this.
## Make sure you have required R packages, which are NOT checked: ggplot2, reshape2, grid, gridExtra, GMD, labeling, and RColorBrewer.

cd footLoop
./check_software.pl


## Download hg19.fa.gz from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz, gunzip it, then put into footLoop folder.

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz

# 3. Example run:
# Download example.tar.gz and move everything in it into footLoop folder.
tar zxvf example.tar.gz
mv example/* footLoop/
cd footLoop
./footLoop.pl -r PCB190425.fq.gz -n PCB190425_MAP -g hg19.fa -l PCB190425 -x -10 -y 10 -i geneIndexes.bed
./footPeak.pl -n PCB190425_MAP -o PCB190425_PEAK
./footClust.pl -n PCB190425_PEAK
./footPeak_graph.pl -n PCB190425_PEAK
./footPeak_GCprofile.pl -n PCB190425_PEAK -i geneIndexes.bed
./footPeak_GTF.pl -n PCB190425_PEAK
./footStats.pl -n PCB190425_PEAK
./footRepro.pl -n PCB190425_PEAK # Will not work, as example fastq only contains 1 biorep for each gene.

#4. Results folders:

## BED files of peaks: PCB190425_PEAK/PEAKS_GENOME/
## GTF files of peaks: PCB190425_PEAK/GTF/PEAK/
## PNG footprint of peaks: PCB190425_PEAK/PNG/PEAK/
## PDF of GC profile: PCB190425_PEAK/GCPROFILE/PDF/PEAK/
## PDF of Reproducibility Correlation: PCB190425_PEAK/FOOTREPRO/*.pdf
## Stats summary: PCB190425_PEAK/99_FOOTSTATS/1_PEAKSTATS.TXT


```

# 1. USAGE

## 1a. Download/Install

`git clone https://github.com/srhartono/footLoop`

### Other required packages/softwares (and versions used for this pipeline)

Softwares (also in softwares/ folder):

**NEWER VERSIONS of these softwares MIGHT NOT WORK due to output file differences. If you get errors it's very likely due to this.**

**These are also avaiable in softwares folder. Make sure to install, then put their folder location into $PATH after unzipping!**


- [bedtools2 (v2.25.0)](https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz)
- [bowtie2 (v2.2.6)](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.6/bowtie2-2.2.6-source.zip/download)
- [bismark2 (v0.20.0)](https://www.bioinformatics.babraham.ac.uk/projects/bismark/bismark_v0.20.0.tar.gz)
- [R (v3.6.1)](https://cloud.r-project.org/src/base/R-3/R-3.6.1.tar.gz)
- [samtools (v0.1.19-96b5f2294a)](https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2/download)

R packages:

- RColorBrewer v1.1-2 
- gridExtra v2.3
- labeling v0.3
- reshape2 v1.4.3    
- ggplot2 v3.1.0
- GMD v0.3.3         

To Install R packages, type these in an R session:

```
install.packages("RColorBrewer")
install.packages("gridExtra")
install.packages("labeling")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("gplots")
```

To install GMD:

```
cd footLoop/softwares
tar zxvf GMD_0.3.1.1.tar.gz
R CMD INSTALL GMD
```

## 1b. Mapping (footLoop.pl)

```
footLoop.pl -r <read.fastq> -n <output directory> -l <label> -i <index.bed6> -g <genome fasta>

Options (default are in [ ] brackets)
-x: add x bp on the start of amplicon (strand specific). E.g. add 100bp upstream of start site: -x -100. [-x 0]
-y: [0] Same as –x, but add x bp to the end of amplicon. [-y 0]
-x and -y are useful as sometimes read start/end exceed amplicon start/end by 10-100 basepairs, and Bismark will not map these reads even though the rest of the read maps perfectly to the amplicon.

-l: alphanumeric label
-L: minimum read length (e.g. 500bp: –L 500). Add “p” to make it based on amplicon length (e.g.85%: -L 85p). [-L 85p]
-q: minimum map quality. [-q 0]
-Z: toggle to use non-stringent mapping (-H for more explanation). [off]
-F toggle to redo bismark mapping even if a .sam/.bam file is present in output_dir. [off]
```

## 1c. Peak calling (footPeak.pl)


```
footPeak.pl -n <footLoop output directory> -o <output directory>

Options:
-l: minimum peak length in basepair. [-l 100]
-w: cytosines used to calculate conversion (other bases are ignored)/ [-w 20]
-t: conversion threshold in fraction [-t 0.55]
-G: only process this gene, e.g. -G CALM3. [off]
```

## 1d. Clustering (footClust.pl)


```
footClust.pl -n <footPeak output directory>

Options:
-D: how tight (basepairs) reads in a cluster should be placed relative to each other. [-D 200]
-R: toggle to order reads in reverse if the gene’s strand is negative. [off]
-G: only process this gene, e.g. -G CALM3. [off]
```

## 1e. Graphing (footPeak_graph.pl)


```
footPeak_graph.pl -n <footPeak_output directory>

Options:
-c: Include Cytosines in CpG context. [off]
-r: Option to create PNG graphs. [-r 1]
   -r 0: do not run any R scripts
   -r 1: run only relevant R scripts (default)
   -r 2: run ALL R scripts
 -R Option to create PDF graphs. [-R 0]
   -R 0: do not run any R scripts (default)
   -R 1: run only relevant R scripts
   -R 2: run ALL R scripts
-G: only process this gene, e.g. -G CALM3. [off]
```

## 1f. Creating GTF for UCSC (footPeakGTF.pl)

```
footPeakGTF.pl -n <footPeak output directory>

Options:
-G: only process this gene, e.g. -G CALM3. [off]
```

## 1g. Calculating GC Profile (footPeak_GCprofile.pl)


```
footPeak_GCprofile.pl -n <footPeak output directory> -i <gene index file>

Options:
-G: only process this gene, e.g. -G CALM3. [off]
```

## 1h. Stats summary (footStats.pl)


```
footStats.pl -n <footPeak output directory>

Options:
-G: only process this gene, e.g. -G CALM3. [off]

Will produce 0_SUMMARY.TXT and 1_PEAKSTATS.TXT, which contain statistics of the run, e.g. how many peaks in a read.
```


## 1i. Reproducibility (footRepro.pl)

```
footRepro.pl -n <footPeak output directory>

Options:
-c: Include Cytosines in CpG context. [off]
-G: only process this gene, e.g. -G CALM3. [off]
```

`PS: This will only work if the gene have more than 1 bioreps/pacbio runs.`

# 2. Detailed Explanation

## 2A. Mapping

Reads were mapped using Bismark v0.2.0 (Krueger and Andrews, 2011) to their respective amplicon regions, buffered by 10bp off their beginning and ends. We used bismark default setting with slightly relaxed minimum score threshold (--score_min of L,0,-0.3 instead of L,0,-0.2), as our data contains much more cytosine conversions in all contexts (CHH, CHG, and CpG, especially in non-promoter regions), whereas most Cytosine conversions in an average methylation sequencing experiment will only come from unmethylated CpGs. Nevertheless, even with relaxed mapping threshold, there is virtually nonexistent risk of misalignment as not only the read length is high quality and very long (~2kb), these reads were mapped to their own amplicon regions, and there were very little sequences in common between these regions. Furthermore, to ensure high quality read, we only keep reads which length are 95% of their respective amplicon lengths.

## 2B. Strand Assignment

For each read, we assign strand based on their conversion types as follows:
1.	We excluded regions that are very GC poor (C < 6 and G < 6) and regions with very few conversions (C->T < 6 and G->A < 6).
2.	If the number of C->T conversions is within +/- 10% of G->A conversions, then we assigned these as unknown. These regions were analyzed separately (Supplementary Figure SX), but in general there are very few of these regions (Supplementary Figure SX).
3.	Otherwise, reads with more C to T conversions were assigned as non-template/R-loop forming strand, and regions with more G to A conversions were assigned as template strand.
Furthermore, regions with indels often create false positives, as it contain higher amount of mismatched nucleotides adjacent to it. Therefore, to be more conservative, we masked up to 5bp around these regions and we treated these regions as “no data” in downstream calculation.

## 2C. Peak Calling

For peak calling, we performed threshold-based sliding window method. Specifically, we calculated conversion percentage every 20 consecutive cytosines, and only those with at least 55% conversion were called as peak of at least 100 bp lengths. As a control, we measured the level of cytosine conversion in background and non-R-loop forming regions, as well as measuring the G conversions of each R-loop forming region.

## 2D. Clustering

For each gene, R-loop peaks were clustered using their start and stop coordinates using k-means clustering. K was determined automatically by minimizing intra-cluster differences, iterating until minimum within-group distance of each clusters is at most 3.


## 2E. Reproduciblity

For each gene replicate, we quantified distribution of reads in each cluster. These distributions were used to measure Pearson correlation coefficient between each replicates. As control, we shuffled the position of each read around the amplicon region and used its position to determine its cluster. A shuffled read is considered to be inside a specific cluster if their start and end positions fall between +/- 100bp of mean start and end of all reads in that cluster. All reads that weren’t inside any specific cluster were put in an extra cluster. Then we calculated Pearson correlation coefficient in the same manner as described above.

# 3. LICENSE

>footLoop pipeline Version 1.2 
>
>Copyright (C) 2019 Stella Regina Hartono
>
>This program is free software: you can redistribute it and/or modify
>it under the terms of the GNU General Public License as published by
>the Free Software Foundation, either version 3 of the License, or
>(at your option) any later version.
>
>This program is distributed in the hope that it will be useful,
>but WITHOUT ANY WARRANTY; without even the implied warranty of
>MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
>GNU General Public License for more details.
>
>The license can be found at https://www.gnu.org/licenses/gpl-3.0.en.html. 
>By downloading or using this software, you agree to the terms and conditions of the license. 
