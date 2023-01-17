PCA for different subsets of the dataset
================

In this script I create different PCAs for different subsets of the HGDP
summary dataset. The focus is the `copynumber` of every transposon for
each individual.

First, I did a PCA for all the TEs and all the samples. Subsequently, I
excluded the “low quality samples”, identified as samples in which the
scg copynumber variation was high. Since all the outliers in the 1st PCA
were poor quality samples, I did all the other PCAs just for the
high-quality subset.

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
library(ggpubr)

pcr_free <- read_tsv("/Volumes/Temp1/rpianezza/SGDP/ric-documentation/SGDP-no-PCR.tsv")
```

    ## Rows: 261 Columns: 1
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): ID
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
SGDP <- read_tsv("/Volumes/Temp2/rpianezza/SGDP/summary/USEME_SGDP_cutoff") %>% filter(biosample %in% pcr_free$ID)
```

    ## Rows: 470028 Columns: 10
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (7): biosample, sex, pop, country, type, familyname, batch
    ## dbl (3): length, reads, copynumber
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
SGDP_stable_names <- SGDP %>% filter(type=="scg") %>% group_by(biosample, pop, country, sex) %>% dplyr::summarise(var = var(copynumber)) %>% arrange(desc(var)) %>% filter(var<0.01) %>% ungroup() %>% distinct(biosample) %>% pull()
```

    ## `summarise()` has grouped output by 'biosample', 'pop', 'country'. You can
    ## override using the `.groups` argument.

``` r
SGDP_stable <- SGDP %>% filter(biosample %in% SGDP_stable_names)
```

# Function for PCA plotting

``` r
PCA <- function(data, title){
  
  m_variable <- data %>% filter(sex=="male") %>% group_by(familyname) %>% dplyr::summarise(sd = sd(copynumber)) %>% filter(sd>0) %>% select(familyname) %>% pull
  f_variable <- data %>% filter(sex=="female") %>% group_by(familyname) %>% dplyr::summarise(sd = sd(copynumber)) %>% filter(sd>0) %>% select(familyname) %>% pull
  
  m <- filter(data, sex=='male', familyname %in% m_variable)
  f <- filter(data, sex=='female', familyname %in% f_variable)
  m_len <- length(unique(m$familyname))
  f_len <- length(unique(f$familyname))
  males <- length(unique(m$biosample))
  females <- length(unique(f$biosample))

  f_matrix<-matrix(as.vector(f$copynumber),nrow=females,ncol=f_len,byrow=T)
  f_fram<-data.frame(f_matrix)
  names(f_fram)<-unique(f$familyname)
  f_matrixcont<-matrix(as.vector(f$country),nrow=females,ncol=f_len,byrow=T)
  f_framcont<-data.frame(f_matrixcont)
  f_contcol<-c(f_framcont$X1)
  fHGDP.pca <- prcomp(f_fram, center = TRUE, scale = TRUE)
  
  m_matrix<-matrix(as.vector(m$copynumber),nrow=males,ncol=m_len,byrow=T)
  m_fram<-data.frame(m_matrix)
  names(m_fram)<-unique(m$familyname)
  m_matrixcont<-matrix(as.vector(m$country),nrow=males,ncol=m_len,byrow=T)
  m_framcont<-data.frame(m_matrixcont)
  m_contcol<-c(m_framcont$X1)
  mHGDP.pca <- prcomp(m_fram, center = TRUE, scale = TRUE)
  
 m_contcol <- factor(m_contcol, levels = c("Africa", "America", "Central Asia and Siberia", "East Asia", "South Asia", "West Eurasia", "Oceania"))
  
  f_contcol <- factor(f_contcol, levels = c("Africa", "America", "Central Asia and Siberia", "East Asia", "South Asia", "West Eurasia", "Oceania"))
  
  library(ggbiplot)
  f_PCA <- ggbiplot(fHGDP.pca, var.axes=FALSE, groups = f_contcol, ellipse = FALSE)+ ggtitle("Females")+ theme(plot.title = element_text(size = 8, hjust = 0.5)) 
m_PCA <- ggbiplot(mHGDP.pca, var.axes=FALSE, groups = m_contcol, ellipse = FALSE)+ ggtitle("Males")+ theme(plot.title = element_text(size = 8, hjust = 0.5)) 
  
figure <- ggarrange(f_PCA, m_PCA, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom", align = "hv", font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

annotate_figure(figure, top = text_grob(title, color = "black", size = 20), fig.lab = "")
}
```

# PCA for the most variable TEs

``` r
TE_cutoff <- filter(SGDP, type=='te') %>% group_by(familyname, sex) %>% mutate(max=max(copynumber), min=min(copynumber)) %>% mutate(diff = max-min, ratio = max/min)

out_abs <- filter(TE_cutoff, diff>200 & diff<Inf)
out_rel <- filter(TE_cutoff, ratio>5 & ratio<Inf & max>1.5)
out_abs_names <- c(unique(out_abs$familyname))
out_rel_names <- c(unique(out_rel$familyname))
out_names <- c(out_abs_names, out_rel_names[!(out_rel_names %in% out_abs_names)]) %>% sort()

most_variable<-subset(SGDP, type=="te") %>% filter(familyname %in% out_names)

TE_stable <- filter(SGDP_stable, type=="te") %>% group_by(familyname, sex) %>% mutate(max=max(copynumber), min=min(copynumber)) %>% mutate(diff = max-min, ratio = max/min)

out_abs_s <- filter(TE_stable, diff>200 & diff<Inf)
out_rel_s <- filter(TE_stable, ratio>5 & ratio<Inf & max>1.5)
out_abs_names_s <- c(unique(out_abs_s$familyname))
out_rel_names_s <- c(unique(out_rel_s$familyname))
out_names_s <- c(out_abs_names_s, out_rel_names_s[!(out_rel_names_s %in% out_abs_names_s)]) %>% sort()

most_variable_s <- subset(SGDP_stable, type=="te") %>% filter(familyname %in% out_names_s)

PCA(TE_cutoff, "All the repetitive sequences")
```

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

![](03_SGDP_PCA_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
PCA(TE_stable, "All the repetitive sequences (only high quality samples)")
```

![](03_SGDP_PCA_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
PCA(most_variable_s, "Most variable repetitive sequences")
```

![](03_SGDP_PCA_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

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
families_in_SGDP <- filter(SGDP_stable, type=="te") %>% distinct(familyname) %>% pull()

classification <- filter(classification_tot, familyname %in% families_in_SGDP)

SGDP_class <- full_join(SGDP_stable, classification, by="familyname")

write_tsv(SGDP_class, "/Volumes/Temp1/rpianezza/SGDP/summary/SGDP_classified.tsv", col_names = TRUE)

DNA_names <- c("Crypton", "hAT", "Helitron", "Kolobok", "Mariner/Tc1", "Merlin", "MuDR", "piggyBac", "DNA transposon")
LINE_names <- c("L1", "CR1", "L2", "Crack", "RTE", "RTEX", "R4", "Vingi", "Tx1", "Penelope")
SINE_names <- c("SINE1/7SL", "SINE2/tRNA", "SINE3/5S", "SINE")
LTR_names <- c("ERV1", "ERV2", "ERV3", "Gypsy", "Endogenous Retrovirus", "LTR Retrotransposon", "Long terminal repeat", "Non-LTR Retrotransposon")
satellites_names <- c("Satellite", "satellite", "SAT")

(classification <- SGDP_class %>% mutate(class = case_when(superfamily %in% DNA_names ~ "DNA", superfamily %in% LINE_names ~ "LINE", superfamily %in% SINE_names ~ "SINE", superfamily %in% LTR_names ~ "LTR", superfamily %in% satellites_names ~ "satellite")))
```

    ## # A tibble: 236,717 × 13
    ##    biosam…¹ sex   pop   country type  famil…² length reads copyn…³ batch super…⁴
    ##    <chr>    <chr> <chr> <chr>   <chr> <chr>    <dbl> <dbl>   <dbl> <chr> <chr>  
    ##  1 SAMEA33… fema… Kinh  East A… scg   chr1:9…   4152 2824.   0.962 flo   <NA>   
    ##  2 SAMEA33… fema… Kinh  East A… scg   chr1:9…   5136 2964.   0.816 flo   <NA>   
    ##  3 SAMEA33… fema… Kinh  East A… scg   chr1:1…   3064 1971.   0.910 flo   <NA>   
    ##  4 SAMEA33… fema… Kinh  East A… scg   chr1:1…   3239 2256.   0.985 flo   <NA>   
    ##  5 SAMEA33… fema… Kinh  East A… scg   chr1:1…   4035 2736.   0.959 flo   <NA>   
    ##  6 SAMEA33… fema… Kinh  East A… scg   chr1:1…   2500 1810.   1.02  flo   <NA>   
    ##  7 SAMEA33… fema… Kinh  East A… scg   chr1:1…   2599 1826.   0.993 flo   <NA>   
    ##  8 SAMEA33… fema… Kinh  East A… scg   chr1:1…   2124 1297.   0.863 flo   <NA>   
    ##  9 SAMEA33… fema… Kinh  East A… scg   chr1:2…   6284 3834.   0.862 flo   <NA>   
    ## 10 SAMEA33… fema… Kinh  East A… scg   chr1:2…   3222 2167.   0.951 flo   <NA>   
    ## # … with 236,707 more rows, 2 more variables: shared_with <chr>, class <chr>,
    ## #   and abbreviated variable names ¹​biosample, ²​familyname, ³​copynumber,
    ## #   ⁴​superfamily

``` r
LINE <- filter(classification, type=="te", class=="LINE")
PCA(LINE, "LINEs")
```

![](03_SGDP_PCA_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
DNA <- filter(classification, type=="te", class=="DNA")
PCA(DNA, "DNA transposons")
```

![](03_SGDP_PCA_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
LTR <- filter(classification, type=="te", class=="LTR")
PCA(LTR, "LTR retrotransposons")
```

![](03_SGDP_PCA_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

``` r
simple_repeats <- filter(classification, type=="te", class=="satellite")
PCA(simple_repeats, "Satellites")
```

![](03_SGDP_PCA_files/figure-gfm/unnamed-chunk-6-4.png)<!-- -->

I notice that the pattern previously described is evident in **non-LTR
retrotransposons** as well as in **DNA transposons**, but not in
**simple repeats**. We do not expect simple repeats to rapidly spread
into different populations as TEs, so this is an expected result, but
still nice to see.

``` r
L1 <- filter(classification, type=="te", superfamily=="L1")
PCA(L1, "LINE-1 retrotransposons")
```

![](03_SGDP_PCA_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
L2 <- subset(classification, type=="te") %>% filter(grepl("L2|L3|L4",familyname))
PCA(L2, "LINE-2/3/4 retrotransposons")
```

![](03_SGDP_PCA_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->
