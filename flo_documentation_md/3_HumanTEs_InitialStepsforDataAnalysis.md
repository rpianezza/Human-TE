Preparing Data analysis of the HGDP dataset
================

# Unbiased Population Variation of Human Transposable Elements - Script 3

This is the third out of eight scripts describing the creation and
analysis of the dataset of human TE abundance. This script is used to
describe the loading of the main HGDP dataset and the most commonly used
subsets. The original dataset creation is described in sourceforge at
<https://sourceforge.net/p/human-te-dynamics/wiki/Home/> . The dataset
was further subsetted to exclude TEs with very low copy numbers. This
process is further explained in script 2.

For all analyses, you would likely want to operate in tidyverse. Thus
include:

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

If you are doing analyses with R or RStudio, you could consider changing
your working directory to your working copy of the svn data repository,
e.g. with a command like this:

``` r
setwd('/Users/rpianezza/TE')
```

However, if you want to change your working directory for Rmarkdown, you
have to include the following code chung in the ‘setup’ at the very
beginning of a markdown document:

``` r
knitr::opts_knit$set(root.dir = '/Users/rpianezza/TE')
```

To read in the dataset and name the columns, run the following commands
(adjust the path if necessary):

``` r
HGDPcutoff<-read_delim("/Users/rpianezza/TE/summary-HGDP/USEME_HGDP_complete_reflib6.2_mq10_batchinfo_cutoff0.01.txt",comment="#")
```

    ## Rows: 1394352 Columns: 10
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (7): ID, Pop, sex, Country, type, familyname, batch
    ## dbl (3): length, reads, copynumber
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
names(HGDPcutoff)<-c("ID","Pop","sex","Country","type","familyname","length","reads","copynumber","batch")
```

To create plots that can be further adjusted in softwares like Adobe
Illustrator, you can use the postscript command:

``` r
postscript(file="plot.eps")
#ggplot(plotobject)
dev.off()
```

To analyse subsets in plots, you can use tidyverse’s ‘group_by’ method,
e.g.

``` r
HGDPcutoff %>% group_by(sex,familyname)
```

Alternatively, but less ideally (though you will see this method in some
of my scripts, it should not be preferred over the tidyverse method),
you can create subsets from your initial dataset, e.g.:

``` r
ftecutoff<-subset(HGDPcutoff, sex=="female" & type=="te")
mtecutoff<-subset(HGDPcutoff, sex=="male" & type=="te")
```
