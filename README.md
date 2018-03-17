### V2.96
- Fixed bugs on peak calling, naming, pathing, labeling (dashed)
- Added -v on footLoop.pl and footPeak.pl to do version checking

### V2.95
- Created footPeakGTF.pl which will turn footPeak peaks (in .CALL/*.out) into gtf format in GTF/ folder
- Fixed LOTS of bugs in all scripts (pathing bugs, naming bugs)
- Beautified scripts to make them more easily maintained
- Reworked folders so now all in/out folders and files are inside footPeak's -o output folder. E.g. footPeakGTF.pl will put output file in <footPeakFolder>/GTF. Cluster fils will be put in <footPeakFolder>/FOOTCLUST

Working on:
- 75% done on HMM based peak caller. The framework is done, just need fine tuning and bug testing (bug testing might take longer time depending on how many bugs)
- Make track hub that can sync with official UCSC sessions (e.g hg19)

To do:
- Create rDNA track hub if rDNA can't be synced with official UCSC sessions above.
- Create better classes for each footprint "types" (e.g. trailing peaks)

Known issues:
- Clustering is still hardcoded into 5 clusters.
- Sometimes in a cluster, if a minority of reads have 2 peaks that are separated by long distance, the aggregate cluster peak will be long as well. This will be fixed by making a smarter clustering.
- There might be unexpected bugs with footLoop -l function (label) as this haven't been tested very thoroughly.


### How to clone:

`git clone ssh://git@github.com/srhartono/footLoop.git`

**IMPORTANT! When it ask for this, type "yes"**:

`Are you sure you want to continue connecting (yes/no)?` **yes**


### 2.0 - 2.72

Massive rework on scrips

- footLoop.pl deals with read map/filter
- footPeak.pl deals with peak calling


Improved graphs and added stats


### Initial Commit 4/28/2017

git checkout ID: 1735a7ec4a8e11b76e588a31dc19fdc0f8a4d83c

Command to checkout: `git checkout -b <NEWBRANCH> 1735a7ec4a8e11b76e588a31dc19fdc0f8a4d83c`

md5sum for files:

|md5sum | file |verrsion|
|a638c3b8baee302ade667a3744a83ac3 | footLoop.pl| z4|
|7f6271300053b429099043bc0f8bdfb4 | footPeak.pl| z2|
|0acc108024a11a688ef3bb74ec49f5af | foot.pm| |
|76fb117bb0005dac1f154ef5f5a8a42c | putToBrows.pl| z2|


