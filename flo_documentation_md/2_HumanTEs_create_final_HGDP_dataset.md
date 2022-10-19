Final steps to creating the proper dataset for HGDP analyses
================

# Unbiased Population Variation of Human Transposable Elements - Script 2

This is the second out of eight scripts describing the creation and
analysis of the dataset of human TE abundance. This script is the first
script to work with the full mapping dataset of the HGDP. To understand
how the dataset was created, please have a look at the sourceforge
documentation <https://sourceforge.net/p/human-te-dynamics/wiki/Home/> .

The main purpose of this script is to create the final dataset that will
be used for all subsequent analyses.

Also, during my initial analysis, the versions and terminologies used
were not always concordant. I tried my best to make everything coherent
and based on the same dataset for these scripts. If ynathing is unclear,
just contact me or Robert.

I disabled the automatic evaluation of the code, as there is no point in
recreating the dataset each time using absolute paths. But running all
of the code adjusting for the paths of the input file should create the
files necessary for the subsequent analyses (unless indicated
differently, e.g. for the analyses including samples with ancient DNA).

For many analyses and processing steps, we will require the tidyverse
package:

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
    ## ✔ ggplot2 3.3.6      ✔ purrr   0.3.4 
    ## ✔ tibble  3.1.8      ✔ dplyr   1.0.10
    ## ✔ tidyr   1.2.1      ✔ stringr 1.4.1 
    ## ✔ readr   2.1.2      ✔ forcats 0.5.2 
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

Next, you need to set your working directory, e.g. the working directory
to the ‘data’ sourceforge repository. This will have to be adjusted at
each computer. Please also make sure if all other paths used in the
various Markdown scripts match up to your personal setup and adjsut them
if necessary.

``` r
setwd('/Users/rpianezza/TE')
```

To read in the initial dataset, run read_delim on the summary dataset,
likely named
‘USEME_HGDP_complete_reflib6.2_mq10_batchinfo_cutoff0.01copies.txt’:

``` r
HGDP <- read_delim("/Users/rpianezza/TE/summary-HGDP/nocutoff/HGDP_complete_reflib6.2_mq10_batchinfo.txt",delim="\t",col_names=FALSE,comment="#")
```

We assign names for the different contents of the file.

``` r
names(HGDP)<-c("ID","Pop","sex","Country","type","familyname","length","reads","copynumber","batch")
```

**ID**: Referring to the unique HGDP identifier of each sample,
e.g. HGDP00001  
**Pop**: Defining the population were the sample was taken,
e.g. ‘Brahui’  
**Sex**: Sex of the sample, either ‘male’ or ‘female’  
**Country**: Referring to the general area of origin, either a continent
or part of a continent (elements: Central_South_Asia Africa, Oceania,
Europe, Middle_East, America, East_Asia). Originally, there were two
alternative names for samples from Oceania, either ‘Oceania’ or
‘Oceania\_(SGDP),Oceania’. For the sake of simplicity, occurences of
‘Oceania\_(SGDP),Oceania’, were changed to ‘Oceania’.

``` r
for (i in 1:length(HGDP$Country)){
  if (HGDP$Country[i]=="Oceania_(SGDP),Oceania"){
    HGDP$Country[i]<-"Oceania"}}
```

**type**: Defining the type of sequence. Either ‘scg’, ‘scgx’, ‘te’ or
‘krab’.  
**familyname**: Defining the name of the corresponding consensus
sequence.  
**length**: Length of the corresponding consensus sequence.  
**reads**: THE NUMBER OF READS (CHECK AGAIN SOURCEFORGE
DOCUMENTATION).  
**copynumber**: The normalized copynumber of the sequence in the
repsective sample.  
**batch**: Initially defining the computer on which it was conducted,
now hopefully meaningless.

More detailed explanation for ‘batch’: The creation of the different
mapstat files for all HGDP datasets was a long process that was
initially conducted on two different computers. Despite trying to use
the exact same versions of all programs, we discovered that the version
from the two computers had inherent biases regarding abundance. So all
analyses originally conducted on vetlinux04 (marked ‘flo’ in the initial
dataset) were rerun on the vetgrid27 computer. An overview of the
datasets including who ran them and a acheck of the md5 sum of the
initial cram-file including the NGS reads can be found in the
Excel-Sheet

## Outlier removal after mapping

There are various reasons why we would want to modify our data after
mapping. For example, one or various single-copy genes could turn out to
vary significantly more than previously expected. Or specific TE
sequences are basically nonexistent in any of the samples and are thus
i) not interesting to analyze and ii) cause problems in the pairwise
comparisons (e.g. when you then have to divide by 0). Of course, these
issues are inherently different, and thus also need to be addressed
differently. If we remove a TE from the analysis, not much changes. We
can just exclude a bunch of sequences, perform all analyses with the
subset excluding these sequences, and pretend like they were never in
the reference library v6.2 to begin with. On the other hand, excluding
only a single SCG means that we have to recalculate the copy number
estimates for all TEs in all samples. There is a script on Sourceforge
that does exactly that (CHECK AGAIN WHAT THE SCRIPT IS CALLED AND WHERE
IT IS!).

The idea behind the outlier removal conducted here was the above
described simplification of graphical display of pairwise relative
comparisons. In other words, the log-version of plots with very low copy
number estimates did not work properly.

``` r
cutoff <- HGDP %>% group_by(familyname) %>% filter(copynumber <0.01)
uniquecutofffamnames<-unique(cutoff$familyname) 
removetes<-c()
for (i in uniquecutofffamnames){
  if (max(HGDP$copynumber[HGDP$familyname==i])<0.01){
    removetes<-c(removetes,i)}}
`%notin%` <- Negate(`%in%`) # This defines the function %notin% as doing the reverse of the function %in%
HGDPcutoff<-HGDP %>% filter(familyname %notin% removetes)
write.csv(HGDPcutoff,"/Users/rpianezza/TE/summary-HGDP/rp_USEME_HGDP_complete_reflib6.2_mq10_batchinfo_cutoff0.01.txt",row.names = FALSE)
```

The idea in this code is to remove all TEs from the HGDP dataset that do
not possess a single instance out of all individual samples where the
normalized coverage exceeds 0.01 copies. This means that all these TEs
have incredibly low coverage and do never have more than 0.01 copies in
each of the sample. The reason why this approach was chosen over using
the mean of the TE as a cutoff is since with the mean, it might still be
that one or very few samples have a significantly higher abundance of
the TE, but they are just overshadowed with a very low abundance in most
other samples. By using the maximum value as a cutoff, such a scenario
is avoided. The value 0.01 was chosen arbitrarily and a higher value
could probably be chosen, as it is debatable if anything below at least
1 copy is of any interest for these analyses. The TEs removed in this
step but still present in the original reflibrary v62. are listed below:
“X5a_DNA”,“LTR88a”,“LTR81B”,“LTR89”,“LTR81C”,“LTR86B1”,“Eutr3”“,”LTR88b”,“MER104A”,“BSRd”,“MamGypLTR2”,“X11_LINE”,“LTR91”,
“CHARLIE8A”,“CR1_HS”,“LTR86C”,“UCON39”,“LTR88c”,“HERVR”

The dataset HGDPcutoff will be the base dataset for all analyses, unless
specified afterwards. It is stored in the sourceforge repository
