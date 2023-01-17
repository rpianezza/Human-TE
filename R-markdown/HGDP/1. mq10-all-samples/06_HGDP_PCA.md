PCA for different subsets of the dataset
================

In this script I create different PCAs for different subsets of the HGDP
summary dataset. The focus is the `copynumber` of every transposon for
each individual.

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
    ## ✔ ggplot2 3.4.0      ✔ purrr   0.3.4 
    ## ✔ tibble  3.1.8      ✔ dplyr   1.0.10
    ## ✔ tidyr   1.2.1      ✔ stringr 1.4.1 
    ## ✔ readr   2.1.2      ✔ forcats 0.5.2 
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library("ggpubr")
HGDPcutoff <- read_delim("/Volumes/Temp1/rpianezza/TE/summary-HGDP/USEME_HGDP_mq0_cutoff0.01.txt")
```

    ## Rows: 1396835 Columns: 10
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (7): ID, Pop, sex, Country, type, familyname, batch
    ## dbl (3): length, reads, copynumber
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
names(HGDPcutoff) <- c("ID","pop","sex","country","type","familyname","length","reads","copynumber","batch")

HGDPcutoff <- HGDPcutoff %>% mutate(country = recode(country, "Oceania_(SGDP),Oceania"="Oceania"))
```

# Function for PCA plotting

First, I create a function that takes a subset of the dataset (`data`)
and a string as `title` for the plot title and gives back the complete
PCA plot. Note that `ellipses` are shown in all the plots, to remove
them you should modify the function here.

``` r
PCA <- function(data, title){
  m <- filter(data, sex=='male')
  f <- filter(data, sex=='female')
  len <- length(unique(data$familyname))
  males <- length(unique(m$ID))
  females <- length(unique(f$ID))

  f_matrix<-matrix(as.vector(f$copynumber),nrow=females,ncol=len,byrow=T)
  f_fram<-data.frame(f_matrix)

  names(f_fram)<-unique(f$familyname)
  f_matrixcont<-matrix(as.vector(f$country),nrow=females,ncol=len,byrow=T)
  f_framcont<-data.frame(f_matrixcont)
  f_contcol<-c(f_framcont$X1)

  fHGDP.pca <- prcomp(f_fram, center = TRUE, scale = TRUE)
  
  m_matrix<-matrix(as.vector(m$copynumber),nrow=males,ncol=len,byrow=T)
  m_fram<-data.frame(m_matrix)

  names(m_fram)<-unique(m$familyname)
  m_matrixcont<-matrix(as.vector(m$country),nrow=males,ncol=len,byrow=T)
  m_framcont<-data.frame(m_matrixcont)
  m_contcol<-c(m_framcont$X1)

  mHGDP.pca <- prcomp(m_fram, center = TRUE, scale = TRUE)
  
  
  library(ggbiplot)
  f_PCA <- ggbiplot(fHGDP.pca, var.axes=FALSE, groups = f_contcol, ellipse = TRUE)+ ggtitle("Females")+ theme(plot.title = element_text(size = 8, hjust = 0.5)) 
  m_PCA <- ggbiplot(mHGDP.pca, var.axes=FALSE, groups = m_contcol, ellipse = TRUE)+ ggtitle("Males")+ theme(plot.title = element_text(size = 8, hjust = 0.5)) 
  
figure <- ggarrange(f_PCA, m_PCA, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom", align = "hv", font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

annotate_figure(figure, top = text_grob(title, color = "black", size = 20), fig.lab = "")
}
```

# PCA for the most variable TEs

My first idea is to compare the PCA for the `copynumber` of all the
**repetitive sequences** (RepSeq) in the dataset with a more specific
PCA, the one considering only the most **variable** RepSeq. To create
this subset I used the same **filters** applied in scripts 1 and 2.

``` r
TE_cutoff <- filter(HGDPcutoff, type=='te') %>% group_by(familyname, sex) %>% mutate(max=max(copynumber), min=min(copynumber)) %>% mutate(diff = max-min, ratio = max/min)

out_abs <- filter(TE_cutoff, diff>200 & diff<Inf)
out_rel <- filter(TE_cutoff, ratio>5 & ratio<Inf & max>1.5)
out_abs_names <- c(unique(out_abs$familyname))
out_rel_names <- c(unique(out_rel$familyname))
out_names <- c(out_abs_names, out_rel_names[!(out_rel_names %in% out_abs_names)]) %>% sort()

most_variable<-subset(HGDPcutoff, type=="te") %>% filter(familyname %in% out_names)

PCA(HGDPcutoff, "All the repetitive sequences")
```

    ## Warning in matrix(as.vector(m$copynumber), nrow = males, ncol = len, byrow =
    ## T): data length [932910] is not a sub-multiple or multiple of the number of rows
    ## [553]

    ## Warning in matrix(as.vector(m$country), nrow = males, ncol = len, byrow = T):
    ## data length [932910] is not a sub-multiple or multiple of the number of rows
    ## [553]

    ## Loading required package: plyr

    ## ------------------------------------------------------------------------------

    ## You have loaded plyr after dplyr - this is likely to cause problems.
    ## If you need functions from both plyr and dplyr, please load plyr first, then dplyr:
    ## library(plyr); library(dplyr)

    ## ------------------------------------------------------------------------------

    ## 
    ## Attaching package: 'plyr'

    ## The following object is masked from 'package:ggpubr':
    ## 
    ##     mutate

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     arrange, count, desc, failwith, id, mutate, rename, summarise,
    ##     summarize

    ## The following object is masked from 'package:purrr':
    ## 
    ##     compact

    ## Loading required package: scales

    ## 
    ## Attaching package: 'scales'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     discard

    ## The following object is masked from 'package:readr':
    ## 
    ##     col_factor

    ## Loading required package: grid

![](06_HGDP_PCA_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

- **Africans** are clearly separated from the others.
- **Eurasians**, which comprehend `Europe`, `Middle_East`,
  `Central_South_Asia` and `East_Asia`, are always clustered together.
  This makes sense because there is no strong geographic barrier among
  these populations, thus there was a strong gene flow across thousands
  of years which led to a more homogeneous genetic background.
- **Oceanians** are also a separated group, with some individuals closer
  to Africans and others to Eurasians.
- **Americans** are clustered together, slightly separated from
  Eurasians.

All the described patterns are more evident in **females**. In the
males, there are some individuals from America and Eurasia that are
clustered together with the africans. This was previously noted in
Scripts 1-2.

``` r
PCA(most_variable, "Most variable repetitive sequences")
```

![](06_HGDP_PCA_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

Filtering for the most variable sequences does not change the figure so
much. I only notice that the percentage of **explained variability** has
increased for both the plots and both axes.

# PCA for different RepSeq families

``` r
classification_tot <- read_tsv("/Users/rpianezza/TE/ric-documentation-Rmd/other-files/repbase_classification.txt", col_names = c("familyname", "superfamily", "shared_with"))
```

    ## Rows: 1386 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): familyname, superfamily, shared_with
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
families_in_HGDP <- filter(HGDPcutoff, type=="te") %>% distinct(familyname) %>% pull()

classification <- filter(classification_tot, familyname %in% families_in_HGDP)

HGDP_class <- full_join(HGDPcutoff, classification, by="familyname")

write_tsv(HGDP_class, "/Volumes/Temp1/rpianezza/TE/summary-HGDP/HGDP_cutoff_classified.tsv", col_names = TRUE)

DNA_names <- c("Crypton", "hAT", "Helitron", "Kolobok", "Mariner/Tc1", "Merlin", "MuDR", "piggyBac", "DNA transposon")
LINE_names <- c("L1", "CR1", "L2", "Crack", "RTE", "RTEX", "R4", "Vingi", "Tx1", "Penelope")
SINE_names <- c("SINE1/7SL", "SINE2/tRNA", "SINE3/5S", "SINE")
LTR_names <- c("ERV1", "ERV2", "ERV3", "Gypsy", "Endogenous Retrovirus", "LTR Retrotransposon", "Long terminal repeat", "Non-LTR Retrotransposon")
satellites_names <- c("Satellite", "satellite", "SAT")

(classification <- HGDP_class %>% mutate(class = case_when(superfamily %in% DNA_names ~ "DNA", superfamily %in% LINE_names ~ "LINE", superfamily %in% SINE_names ~ "SINE", superfamily %in% LTR_names ~ "LTR", superfamily %in% satellites_names ~ "satellite")))
```

    ## # A tibble: 1,396,835 × 13
    ##    ID       pop   sex   country type  famil…¹ length reads copyn…² batch super…³
    ##    <chr>    <chr> <chr> <chr>   <chr> <chr>    <dbl> <dbl>   <dbl> <chr> <chr>  
    ##  1 HGDP000… Brah… male  Centra… scg   chr1:9…   5136 1105.   0.861 ro    <NA>   
    ##  2 HGDP000… Brah… male  Centra… scg   chr1:1…   3064  832.   1.09  ro    <NA>   
    ##  3 HGDP000… Brah… male  Centra… scg   chr1:1…   3239  901.   1.11  ro    <NA>   
    ##  4 HGDP000… Brah… male  Centra… scg   chr1:1…   4035 1102.   1.09  ro    <NA>   
    ##  5 HGDP000… Brah… male  Centra… scg   chr1:1…   2500  733.   1.17  ro    <NA>   
    ##  6 HGDP000… Brah… male  Centra… scg   chr1:1…   2599  580.   0.894 ro    <NA>   
    ##  7 HGDP000… Brah… male  Centra… scg   chr1:1…   2124  477.   0.899 ro    <NA>   
    ##  8 HGDP000… Brah… male  Centra… scg   chr1:2…   6284 1527.   0.973 ro    <NA>   
    ##  9 HGDP000… Brah… male  Centra… scg   chr1:2…   3222  889.   1.10  ro    <NA>   
    ## 10 HGDP000… Brah… male  Centra… scg   chr1:3…   3698  868.   0.940 ro    <NA>   
    ## # … with 1,396,825 more rows, 2 more variables: shared_with <chr>, class <chr>,
    ## #   and abbreviated variable names ¹​familyname, ²​copynumber, ³​superfamily

This part of the code is work in progress. My idea was to give a bit of
context to each RepSeq in the dataset, adding a column with the family
of the sequence. Anyway, I did not find a way to do that apart from
manually annotate the family of each TE, using the info available in
RepBase. I know that this approach is not only super slow but also error
prone, so I will work on that to improve its reliability.

``` r
LINE <- filter(classification, type=="te", class=="LINE")
PCA(LINE, "LINEs")
```

![](06_HGDP_PCA_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
DNA <- filter(classification, type=="te", class=="DNA")
PCA(DNA, "DNA transposons")
```

![](06_HGDP_PCA_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
LTR <- filter(classification, type=="te", class=="LTR")
PCA(LTR, "LTR retrotransposons")
```

![](06_HGDP_PCA_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

``` r
simple_repeats <- filter(classification, type=="te", class=="satellite")
PCA(simple_repeats, "Satellites")
```

![](06_HGDP_PCA_files/figure-gfm/unnamed-chunk-6-4.png)<!-- -->

I notice that the pattern previously described is evident in **non-LTR
retrotransposons** as well as in **DNA transposons**, but not in
**simple repeats**. We do not expect simple repeats to rapidly spread
into different populations as TEs, so this is an expected result, but
still nice to see.

``` r
L1 <- filter(classification, type=="te", superfamily=="L1")
PCA(L1, "LINE-1 retrotransposons")
```

![](06_HGDP_PCA_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
L2 <- subset(classification, type=="te") %>% filter(grepl("L2|L3|L4",familyname))
PCA(L2, "LINE-2/3/4 retrotransposons")
```

![](06_HGDP_PCA_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

The same pattern is well shown in **LINE-1** retrotransposons, known to
be the most active TE family in humans, but is totally absent in
**LINE-2/3-4** retrotransposons, which are expected to be extint in
humans.
