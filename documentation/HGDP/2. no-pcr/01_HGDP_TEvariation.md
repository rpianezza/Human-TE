HGDP - Analyzing the TEs with the highest variance
================

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

no_pcr_samples <- read_tsv("/Volumes/Temp1/rpianezza/investigation/HGDP-no-PCR/HGDP-only-pcr-free-samples.tsv", col_names = ("ID"))
```

    ## Rows: 676 Columns: 1
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): ID
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
HGDP <- read_csv("/Volumes/Temp1/rpianezza/TE/summary-HGDP/USEME_HGDP_complete_reflib6.2_mq10_batchinfo_cutoff0.01.txt", col_names = c("ID","pop","sex","country","type","familyname","length","reads","copynumber","batch"), skip=1) %>% type_convert() %>% filter(ID %in% no_pcr_samples$ID)
```

    ## Rows: 1394352 Columns: 10
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (7): ID, pop, sex, country, type, familyname, batch
    ## dbl (3): length, reads, copynumber
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   ID = col_character(),
    ##   pop = col_character(),
    ##   sex = col_character(),
    ##   country = col_character(),
    ##   type = col_character(),
    ##   familyname = col_character(),
    ##   batch = col_character()
    ## )

## Absolute variation

``` r
f_MMM <- filter(HGDP, type == "te", sex == "female") %>% group_by(familyname) %>% dplyr::summarise(min = min(copynumber), mean = mean(copynumber), max = max(copynumber))

m_MMM <- filter(HGDP, type == "te", sex == "male") %>% group_by(familyname) %>% dplyr::summarise(min = min(copynumber), mean = mean(copynumber), max = max(copynumber))
```

### Females

``` r
f_outliers_names <- dplyr::mutate(f_MMM, diff = max-min) %>% filter(diff>200 & diff<Inf)

f_outliers <- filter(HGDP, familyname %in% f_outliers_names$familyname, type == "te", sex == "female")
f_outliers <- f_outliers[order(f_outliers$copynumber,decreasing=T),]
f_outliers$familyname<-factor(f_outliers$familyname,levels=unique(f_outliers$familyname))

ggplot(f_outliers, aes(x=familyname, y=log(copynumber))) + geom_boxplot(notch=F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Most variable repetitive sequences - Females") + theme(plot.title = element_text(hjust = 0.5))
```

![](01_HGDP_TEvariation_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

### Males

``` r
m_outliers_names <- mutate(m_MMM, diff = max-min) %>% filter(diff>200 & diff<Inf)

m_outliers <- filter(HGDP, familyname %in% m_outliers_names$familyname, type == "te", sex == "male")
m_outliers <- m_outliers[order(m_outliers$copynumber,decreasing=T),]
m_outliers$familyname<-factor(m_outliers$familyname,levels=unique(m_outliers$familyname))

ggplot(m_outliers, aes(x=familyname, y=log(copynumber))) + geom_boxplot(notch=F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
ggtitle("Most variable repetitive sequences - Males") + theme(plot.title = element_text(hjust = 0.5))
```

![](01_HGDP_TEvariation_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->
\## Relative comparison

### Females

``` r
f_outliers_names <- mutate(f_MMM, ratio = max/min) %>% filter(ratio>2 & ratio<Inf & max>1.5)

f_outliers <- filter(HGDP, familyname %in% f_outliers_names$familyname, type == "te", sex == "female")
f_outliers <- f_outliers[order(f_outliers$copynumber,decreasing=T),]
f_outliers$familyname<-factor(f_outliers$familyname,levels=unique(f_outliers$familyname))

ggplot(f_outliers, aes(x=familyname, y=log(copynumber))) + geom_boxplot(notch=F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Most variable repetitive sequences - Relative comparison - Females") + theme(plot.title = element_text(hjust = 0.5))
```

![](01_HGDP_TEvariation_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Now we have a dataset containing only the TEs with the highest
differences in abundance between minimum and maximum value.

### Males

``` r
m_outliers_names <- mutate(m_MMM, ratio = max/min) %>% filter(ratio>2 & ratio<Inf & max>1.5)

m_outliers <- filter(HGDP, familyname %in% m_outliers_names$familyname, type == "te", sex == "male")
m_outliers <- m_outliers[order(m_outliers$copynumber,decreasing=T),]
m_outliers$familyname<-factor(m_outliers$familyname,levels=unique(m_outliers$familyname))

ggplot(m_outliers, aes(x=familyname, y=log(copynumber))) + geom_boxplot(notch=F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Most variable repetitive sequences - Relative comparison - Males") + theme(plot.title = element_text(hjust = 0.5))
```

![](01_HGDP_TEvariation_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Note that the L1ME5 TE, which was always present in the analysis with
all the samples, is completely disappeared from the absolute variant TE
plots. This TE was previously identified as having an uneven coverage,
with just one small region that was present at very high frequency that
caused its high variability in copynumber estimate. Furthermore, this
sequence was identified as AT-microsatellite, thus confirming the fact
that PCR is having a lot of problem in dealing with satellites.
