Analyzing the TEs with the highest relative difference between their
minimum and maximum abundance estimate
================

This is the second script written by me, Riccardo. This script works
with the version of the HGDP dataset created in Script 2 by Florian. We
perform all analyses separately for males and females, as we established
that there are significant differences between the sexes in Script 5
from Florian.

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

``` r
library("ggpubr")
#setwd(~/human-data)
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

females <- filter(HGDPcutoff, sex == 'female')
males <- filter(HGDPcutoff, sex == 'male')
```

The idea of the first part of the code is to plot the **general
distribution of TE abundance** estimates in a boxplot. First I will
perform the analysis on the female subset, later on the male subset
using the same settings.

Looking at this plot for all 965 TEs seems excessive (and not
informative), so I created a subset. First, I identify the **minimum**,
**mean** and **maximum** value for each TE sequence. This code chunk
creates a data frame containing respectively the name of the TE and its
minimum, mean and maximum value in the female subset of the data set.
Everything containing minimum, mean and maximum value for each TE is
named from now `MMM`.

``` r
f_MMM <- filter(HGDPcutoff, type == "te", sex == "female") %>% group_by(familyname) %>% summarise(min = min(copynumber), mean = mean(copynumber), max = max(copynumber))

m_MMM <- filter(HGDPcutoff, type == "te", sex == "male") %>% group_by(familyname) %>% summarise(min = min(copynumber), mean = mean(copynumber), max = max(copynumber))
```

## Females

From the whole dataset, I want to create a subset containing only the
TEs with the highest relative differences between minimum and maximum
value. I only kept TEs which satisfy these conditions:

- The **Fold difference** between the minimum and the maximum copy
  number value should be greater than a chosen values (2, in this case)
- The **Fold difference** must be lower than `Infinite` as otherwise
  this caused some issues
- The **maximum** should be greater than a specific value (1.5, in this
  case). This condition has been added to avoid TEs with high
  fold-difference but very low copynumber. For example, a TE with `max`
  copynumber of 0.4 and `min` copynumber of 0.1 met the condition of
  “fold difference \> 2”, but it’s not informative at all.

For the production of the plot we do want all values of the respective
TE, and not just its minimum and maxmimum value. Thus, we need to create
a subset of the full dataset containing only the TEs with the **highest
fold difference** between minimum and maximum value. As we already have
a vector containing the names of the TEs we want for the female dataset,
we just need to subset the dataset. Then we order the dataset and make
sure that the order is correctly displayed in the follwoing plot.

``` r
f_outliers_names <- mutate(f_MMM, ratio = max/min) %>% filter(ratio>2 & ratio<Inf & max>1.5)

f_outliers <- filter(HGDPcutoff, familyname %in% f_outliers_names$familyname, type == "te", sex == "female")
f_outliers <- f_outliers[order(f_outliers$copynumber,decreasing=T),]
f_outliers$familyname<-factor(f_outliers$familyname,levels=unique(f_outliers$familyname))

ggplot(f_outliers, aes(x=familyname, y=log(copynumber))) + geom_boxplot(notch=F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
```

![](02_HGDP_Relative-TEvariation_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

Now we have a dataset containing only the TEs with the highest
differences in abundance between minimum and maximum value.

## Males

``` r
m_outliers_names <- mutate(m_MMM, ratio = max/min) %>% filter(ratio>2 & ratio<Inf & max>1.5)

m_outliers <- filter(HGDPcutoff, familyname %in% m_outliers_names$familyname, type == "te", sex == "male")
m_outliers <- m_outliers[order(m_outliers$copynumber,decreasing=T),]
m_outliers$familyname<-factor(m_outliers$familyname,levels=unique(m_outliers$familyname))

ggplot(m_outliers, aes(x=familyname, y=log(copynumber))) + geom_boxplot(notch=F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
```

![](02_HGDP_Relative-TEvariation_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

## Details for each relevant TE

In this section I go in details analyzing the **copynumber geographical
distribution** for the TE families previously selected in the subsets. I
decided to perform the analysis only for the sequences not previously
found in script 1.

First, I analysed the difference between the two sex-related subsets. I
found that the **females** dataset (`fTEoutlier_order`) has 19 different
`familyname`, while the **males** has 20. There are 3 sequence families
present in the females and not in the males:

- `LSAU`
- `GSATII`
- `Eutr8`

And 4 sequence families present in the males and not in the females:

- `EUTREP13`
- `AmnSINE1_HS`
- `EuthAT-N1`
- `KER`

``` r
f_familynames <- f_outliers$familyname
m_familynames <- m_outliers$familyname

# TEs present in the female subset and not in the males
length(unique(f_familynames))
```

    ## [1] 19

``` r
setdiff(f_familynames,m_familynames)
```

    ## [1] "LSAU"   "GSATII" "Eutr8"

``` r
# TEs present in the male subset and not in the females
length(unique(m_familynames))
```

    ## [1] 20

``` r
setdiff(m_familynames,f_familynames)
```

    ## [1] "EUTREP13"    "AmnSINE1_HS" "EuthAT-N1"   "KER"

I used the same function written in the Script 1:

``` r
plotTEfamily <- function(data, sex, famname, binwidht, x_title, y_title, x_numbers, y_numbers){
filtered <- filter(data, familyname==famname)
ggplot(data = filtered, mapping = aes(x = copynumber, fill = Country)) +
  geom_histogram(binwidth = binwidht) + 
  ggtitle(paste0(sex,' - ',famname)) + theme(plot.title = element_text(size = 8, hjust = 0.5)) +
  {if(x_title=='n'){
  theme(axis.title.x=element_blank())}} +
  {if(y_title=='n'){
  theme(axis.title.y=element_blank())}} +
  {if(x_numbers=='n'){
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())}} +
{if(y_numbers=='n'){
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())}}
}
```

### Non-LTR retrotransposons

#### LINE-1

The LINE1 variant **L1ME5** appears in both female and male plots as one
of the TEs with the highest variance. Here I look in details at its
distribution in the two datasets.

``` r
(mL1ME5_plot<-plotTEfamily(males, 'Males', 'L1ME5', 0.5, 'y', 'y', 'y', 'y'))
```

![](02_HGDP_Relative-TEvariation_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
(fL1ME5_plot<-plotTEfamily(females, 'Females', 'L1ME5', 0.5, 'y', 'y', 'y', 'y'))
```

![](02_HGDP_Relative-TEvariation_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

**Sex distribution**: we notice a clear separation in both sexes between
low copynumber and high copynumber individuals. Relatively to the
females, there are much more male with low copy number.

In script 5, Florian shows that females has higher copynumber relative
to males for this TE.

**Geographic distribution**: we can see a very interesting pattern.
Almost every **african** individual has very low copy number, while the
high copynumber individuals are almost all **out of Africa** (OOA).

I think this is consistent with a model in which there was a strong
bottleneck due to a migration OOU by a small population around 60.000
years ago. It is known that TEs are usually spreading faster in
bottlenecked populations. Thus, we can speculate that at least one big
L1ME5 invasion was triggered after the OOU.

In some way, it seems that this invasion ignored some Y chromosomes but
hit almost every X chromosomes. This would explain the “mixed”
geographic distribution that we see in males with low copynumber. In
other words, we still have many males with few insertions because the Y
chromosome is almost “free”, while only the african females has low
copynumber, and all the OOU females have high copynumber due to the
invasion.

#### AmnSINE_Hs

``` r
(mAmnSINE1_HS_plot<-plotTEfamily(males, 'Males', 'AmnSINE1_HS', 0.1, 'y', 'y', 'y', 'y'))
```

![](02_HGDP_Relative-TEvariation_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
(fAmnSINE1_HS_plot<-plotTEfamily(females, 'Females', 'AmnSINE1_HS', 0.1, 'y', 'y', 'y', 'y'))
```

![](02_HGDP_Relative-TEvariation_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

This sequence entered the dataset only for 1 outliers, seen in the male
plot from `Central South Asia`. I think this is not interesting.

### Endogenous retroviruses

`EUTREP13` and `HERVI`.

``` r
# EUTREP13
mEUTREP13_plot<-plotTEfamily(males, 'Males', 'EUTREP13', 0.2, 'n', 'n', 'y', 'y')
fEUTREP13_plot<-plotTEfamily(females, 'Females', 'EUTREP13', 0.2, 'n', 'n', 'y', 'y')
```

# HERVI

``` r
mHERVI_plot<-plotTEfamily(males, 'Males', 'HERVI', 0.2, 'n', 'n', 'y', 'y')
fHERVI_plot<-plotTEfamily(females, 'Females', 'HERVI', 0.2, 'n', 'n', 'y', 'y')
```

``` r
m_ERV_figure <- ggarrange(mHERVI_plot, mEUTREP13_plot, ncol = 1, nrow = 2, common.legend = TRUE, legend = "right", align = "hv", font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

m_ERV_final <- annotate_figure(m_ERV_figure, left = text_grob("Count", color = "black", rot = 90), bottom = text_grob("Copynumber", color = "black"), top = text_grob("Males", color = "black", size = 20), fig.lab = "")


f_ERV_figure <- ggarrange(fHERVI_plot, fEUTREP13_plot, ncol = 1, nrow = 2, common.legend = TRUE, legend = "right", align = "hv", font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f_ERV_final <- annotate_figure(f_ERV_figure, left = text_grob("Count", color = "black", rot = 90), bottom = text_grob("Copynumber", color = "black"), top = text_grob("Females", color = "black", size = 20), fig.lab = "")

f_ERV_final
```

![](02_HGDP_Relative-TEvariation_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
m_ERV_final
```

![](02_HGDP_Relative-TEvariation_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

For **HERVI**, we see something against the tide. The african
individuals are not in the lower-copynumber part of the distribution,
but more in the middle, in both sexes. We also see possible invasion(s)
in african as well as middle-east populations, since there are
individuals with high copynumber in these groups for both sexes.

It could be interesting to identify the exact origin of these
individuals with high copynumber (maybe they are part of small or
bottlenecked populations as well?). HERV-K are considered the most
active human endogenous retroviruses.

I don’t think that for `EUTREP13` we can deduce somthing interesting
from this plot.

### DNA transposons

`Eutr8`, `EuthAT-N1` and `MER69A`.

``` r
# Eutr8
mEutr8_plot<-plotTEfamily(males, 'Males', 'Eutr8', 0.2, 'n', 'n', 'y', 'y')
fEutr8_plot<-plotTEfamily(females, 'Females', 'Eutr8', 0.2, 'n', 'n', 'y', 'y')

# EuthAT-N1
mEuthAT_N1_plot<-plotTEfamily(males, 'Males', 'EuthAT-N1', 0.2, 'n', 'n', 'y', 'y')
fEuthAT_N1_plot<-plotTEfamily(females, 'Females', 'EuthAT-N1', 0.2, 'n', 'n', 'y', 'y')

# MER69A
mMER69A_plot<-plotTEfamily(males, 'Males', 'MER69A', 0.2, 'n', 'n', 'y', 'y')
fMER69A_plot<-plotTEfamily(females, 'Females', 'MER69A', 0.2, 'n', 'n', 'y', 'y')
```

``` r
m_DNA_figure <- ggarrange(mEutr8_plot, mEuthAT_N1_plot, mMER69A_plot, ncol = 1, nrow = 3, common.legend = TRUE, legend = "right", align = "hv", font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

m_DNA_final <- annotate_figure(m_DNA_figure, left = text_grob("Count", color = "black", rot = 90), bottom = text_grob("Copynumber", color = "black"), top = text_grob("Males", color = "black", size = 20), fig.lab = "")


f_DNA_figure <- ggarrange(fEutr8_plot, fEuthAT_N1_plot, fMER69A_plot, ncol = 1, nrow = 3, common.legend = TRUE, legend = "right", align = "hv", font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f_DNA_final <- annotate_figure(f_DNA_figure, left = text_grob("Count", color = "black", rot = 90), bottom = text_grob("Copynumber", color = "black"), top = text_grob("Females", color = "black", size = 20), fig.lab = "")

f_DNA_final
```

![](02_HGDP_Relative-TEvariation_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
m_DNA_final
```

![](02_HGDP_Relative-TEvariation_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->

The copynumber for these 3 DNA transposons is not changing that much. I
can only note that the **africans** are always less abundant in the
right part of the plots. In other words, they have lower copynumber in
general.

### Not classified repseq

There are some sequences which are included in RepBase but that are not
characterized yet. In our subset we have: `Eutr11`, `UCON75`, `UCON54`,
`MARE8`. See <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5753468/> for
more information.

``` r
# Eutr11
mEutr11_plot<-plotTEfamily(males, 'Males', 'Eutr11', 0.2, 'n', 'n', 'y', 'y')
fEutr11_plot<-plotTEfamily(females, 'Females', 'Eutr11', 0.2, 'n', 'n', 'y', 'y')

# UCON75
mUCON75_plot<-plotTEfamily(males, 'Males', 'UCON75', 0.2, 'n', 'n', 'y', 'y')
fUCON75_plot<-plotTEfamily(females, 'Females', 'UCON75', 0.2, 'n', 'n', 'y', 'y')

# UCON54
mUCON54_plot<-plotTEfamily(males, 'Males', 'UCON54', 0.2, 'n', 'n', 'y', 'y')
fUCON54_plot<-plotTEfamily(females, 'Females', 'UCON54', 0.2, 'n', 'n', 'y', 'y')

# MARE8
mMARE8_plot<-plotTEfamily(males, 'Males', 'MARE8', 0.2, 'n', 'n', 'y', 'y')
fMARE8_plot<-plotTEfamily(females, 'Females', 'MARE8', 0.2, 'n', 'n', 'y', 'y')
```

``` r
m_NI_figure <- ggarrange(mEutr11_plot, mUCON75_plot, mUCON54_plot, mMARE8_plot, ncol = 2, nrow = 2, common.legend = TRUE, legend = "right", align = "hv", font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

m_NI_final <- annotate_figure(m_NI_figure, left = text_grob("Count", color = "black", rot = 90), bottom = text_grob("Copynumber", color = "black"), top = text_grob("Males", color = "black", size = 20), fig.lab = "")


f_NI_figure <- ggarrange(fEutr11_plot, fUCON75_plot, fUCON54_plot, fMARE8_plot, ncol = 2, nrow = 2, common.legend = TRUE, legend = "right", align = "hv", font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f_NI_final <- annotate_figure(f_NI_figure, fUCON75_plot, left = text_grob("Count", color = "black", rot = 90), bottom = text_grob("Copynumber", color = "black"), top = text_grob("Females", color = "black", size = 20), fig.lab = "")

f_NI_final
```

![](02_HGDP_Relative-TEvariation_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
m_NI_final
```

![](02_HGDP_Relative-TEvariation_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->
Their copynumber is not changing that much. I don’t think these are
relevant in our analysis.

### Satellites

I also performed the same analysis for the two **satellites** present in
the subdataset, `TAR1`and `GSATII`. Again, I’m not sure if satellites
are relevant in our analyses.

``` r
# TAR1
mTAR1_plot<-plotTEfamily(males, 'Males', 'TAR1', 1, 'n', 'n', 'y', 'y')
fTAR1_plot<-plotTEfamily(females, 'Females', 'TAR1', 1, 'n', 'n', 'y', 'y')

# GSATII
mGSATII_plot<-plotTEfamily(males, 'Males', 'GSATII', 10, 'n', 'n', 'y', 'y')
fGSATII_plot<-plotTEfamily(females, 'Females', 'GSATII', 10, 'n', 'n', 'y', 'y')
```

``` r
m_SAT_figure <- ggarrange(mTAR1_plot, mGSATII_plot, ncol = 1, nrow = 2, common.legend = TRUE, legend = "right", align = "hv", font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

m_SAT_final <- annotate_figure(m_SAT_figure, left = text_grob("Count", color = "black", rot = 90), bottom = text_grob("Copynumber", color = "black"), top = text_grob("Males", color = "black", size = 20), fig.lab = "")


f_SAT_figure <- ggarrange(fTAR1_plot, fGSATII_plot, ncol = 1, nrow = 2, common.legend = TRUE, legend = "right", align = "hv", font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f_SAT_final <- annotate_figure(f_SAT_figure, fUCON75_plot, left = text_grob("Count", color = "black", rot = 90), bottom = text_grob("Copynumber", color = "black"), top = text_grob("Females", color = "black", size = 20), fig.lab = "")

f_SAT_final
```

![](02_HGDP_Relative-TEvariation_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
m_SAT_final
```

![](02_HGDP_Relative-TEvariation_files/figure-gfm/unnamed-chunk-17-2.png)<!-- -->

- `TAR1`: **africans** are, as always, in the left part of the plot. We
  can also notice some individuals from `Middle East` with high
  copynumber in both sexes.
- `GSATII`: we see the **bimodal** distribution already noticed for many
  TEs. The geographical distribution pattern is also matched, as well as
  the sex-specific pattern. This is the only analysed satellite showing
  clearly this patterns.
