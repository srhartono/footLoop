FootLoop method

##0. Download/Install

`git clone https://github.com/srhartono`

##1. Mapping (footLoop.pl)

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

##2. Peak calling (footPeak.pl)


```
footPeak.pl -n <footLoop output directory> -o <output directory>

Options:
-l: minimum peak length in basepair. [-l 100]
-w: cytosines used to calculate conversion (other bases are ignored)/ [-w 20]
-t: conversion threshold in fraction [-t 0.55]
-G: only process this gene, e.g. -G CALM3. [off]
```

##3. Clustering (footClust.pl)


```
footClust.pl -n <footPeak output directory>

Options:
-D: how tight (basepairs) reads in a cluster should be placed relative to each other. [-D 200]
-R: toggle to order reads in reverse if the gene’s strand is negative. [off]
-G: only process this gene, e.g. -G CALM3. [off]
```

##4. Graphing (footPeak_graph.pl)


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
#5. Creating GTF for UCSC (footPeakGTF.pl)

```
footPeakGTF.pl -n <footPeak output directory>

Options:
-G: only process this gene, e.g. -G CALM3. [off]
```

##EXTRA 1: Calculating GC Profile (footPeak_GCprofile.pl)


```
footPeak_GCprofile.pl -n <footPeak output directory> -i <gene index file>

Options:
-G: only process this gene, e.g. -G CALM3. [off]
```


##EXTRA 2: Reproducibility (footRepro.pl)


```
footRepro.pl -n <footPeak output directory>

Options:
-c: Include Cytosines in CpG context. [off]
-G: only process this gene, e.g. -G CALM3. [off]


Extra:
Creating % Conversion for UCSC (footPeak_conv.pl)
footPeak_conv.pl -n <footPeak output directory>

Options:
-c: Include Cytosines in CpG context. [off]
-G: only process this gene, e.g. -G CALM3. [off]
```


# Detailed Explanation

##1. Mapping

Reads were mapped using Bismark v0.2.0 (Krueger and Andrews, 2011) to their respective amplicon regions, buffered by 10bp off their beginning and ends. We used bismark default setting with slightly relaxed minimum score threshold (--score_min of L,0,-0.3 instead of L,0,-0.2), as our data contains much more cytosine conversions in all contexts (CHH, CHG, and CpG, especially in non-promoter regions), whereas most Cytosine conversions in an average methylation sequencing experiment will only come from unmethylated CpGs. Nevertheless, even with relaxed mapping threshold, there is virtually nonexistent risk of misalignment as not only the read length is high quality and very long (~2kb), these reads were mapped to their own amplicon regions, and there were very little sequences in common between these regions. Furthermore, to ensure high quality read, we only keep reads which length are 95% of their respective amplicon lengths.

##2. Strand Assignment

For each read, we assign strand based on their conversion types as follows:
1.	We excluded regions that are very GC poor (C < 6 and G < 6) and regions with very few conversions (C->T < 6 and G->A < 6).
2.	If the number of C->T conversions is within +/- 10% of G->A conversions, then we assigned these as unknown. These regions were analyzed separately (Supplementary Figure SX), but in general there are very few of these regions (Supplementary Figure SX).
3.	Otherwise, reads with more C to T conversions were assigned as non-template/R-loop forming strand, and regions with more G to A conversions were assigned as template strand.
Furthermore, regions with indels often create false positives, as it contain higher amount of mismatched nucleotides adjacent to it. Therefore, to be more conservative, we masked up to 5bp around these regions and we treated these regions as “no data” in downstream calculation.

##3. Peak Calling

For peak calling, we performed threshold-based sliding window method. Specifically, we calculated conversion percentage every 20 consecutive cytosines, and only those with at least 55% conversion were called as peak of at least 100 bp lengths. As a control, we measured the level of cytosine conversion in background and non-R-loop forming regions, as well as measuring the G conversions of each R-loop forming region.

##4. Clustering

For each gene, R-loop peaks were clustered using their start and stop coordinates using k-means clustering. K was determined automatically by minimizing intra-cluster differences, iterating until minimum within-group distance of each clusters is at most 3.


##5. Reproduciblity

For each gene replicate, we quantified distribution of reads in each cluster. These distributions were used to measure Pearson correlation coefficient between each replicates. As control, we shuffled the position of each read around the amplicon region and used its position to determine its cluster. A shuffled read is considered to be inside a specific cluster if their start and end positions fall between +/- 100bp of mean start and end of all reads in that cluster. All reads that weren’t inside any specific cluster were put in an extra cluster. Then we calculated Pearson correlation coefficient in the same manner as described above.

