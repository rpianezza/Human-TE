Kruskal_Wallis test to determine the significance of various factors on
the estimates of copy numbers in the HGDP dataset
================

# Unbiased Population Variation of Human Transposable Elements - Script 6

This is the sixth out of eight scripts describing the creation and
analysis of the dataset of human TE abundance. This script works with
the version of the HGDP dataset created in Script 2. Here, we perform
all analyses on the whole dataset as well as seperately for males and
females (as we established there are differences betwene the sexes in
Script 5).

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

HGDPcutoff$ID<-as.factor(HGDPcutoff$ID)
HGDPcutoff$Pop<-as.factor(HGDPcutoff$Pop)
HGDPcutoff$sex<-as.factor(HGDPcutoff$sex)
HGDPcutoff$Country<-as.factor(HGDPcutoff$Country)
HGDPcutoff$type<-as.factor(HGDPcutoff$type)
HGDPcutoff$familyname<-as.factor(HGDPcutoff$familyname)
HGDPcutoff$copynumber<-as.numeric(HGDPcutoff$copynumber)
HGDPcutoff$reads<-as.numeric(HGDPcutoff$reads)
HGDPcutoff$batch<-as.factor(HGDPcutoff$batch)


ftecutoff<-subset(HGDPcutoff, sex=="female" & type=="te")
mtecutoff<-subset(HGDPcutoff, sex=="male" & type=="te")
```

The idea of this code is to test for the influence of specific
parameters on the obtained copy number estimates. Kruskal-Wallis is a
nonparametric test, thus we do not need any assumption about the
underlying distribution of the data. Simply put, it is a nonparametric
version of ANOVA. We use this test to determine if different inherent
properties of the data are inherently different in their corresponding
copy numbers. For example, when using the full dataset, we can confirm
previous results showing a strong discrepancy in copy number estimates
between males and females:

``` r
kruskal.test(copynumber ~ sex, data = HGDPcutoff)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  copynumber by sex
    ## Kruskal-Wallis chi-squared = 83.149, df = 1, p-value < 2.2e-16

The Kurskal-Wallis test additionally allows us to explore the potential
effects of other important data properties. For the additional tests, I
again look at the genders separately. When reporting these results, it
is likely that a correction for multiple testing would have to be
applied.

The most interesting test for me was to see the influence of the
continent/region from which the samples were from. This test is
significant with a 5% significance level for both seces, which nicely
illustrates that we do actually observe regional abundance differences.
Though if the significance can survive multiple testing correction I am
not sure about.

``` r
kruskal.test(copynumber ~ Country, data = mtecutoff)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  copynumber by Country
    ## Kruskal-Wallis chi-squared = 13.085, df = 6, p-value = 0.0417

``` r
kruskal.test(copynumber ~ Country, data = ftecutoff)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  copynumber by Country
    ## Kruskal-Wallis chi-squared = 15.345, df = 6, p-value = 0.01773

Contrarily to the Continent/Region-level, we do not observe significance
when looking at variance on the population level:

``` r
kruskal.test(copynumber ~ Pop, data = mtecutoff)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  copynumber by Pop
    ## Kruskal-Wallis chi-squared = 24.679, df = 53, p-value = 0.9997

``` r
kruskal.test(copynumber ~ Pop, data = ftecutoff)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  copynumber by Pop
    ## Kruskal-Wallis chi-squared = 29.022, df = 43, p-value = 0.9492

There is a very significant relation between copy number and length. To
think about if this is expected and if it makes sense is certainly a
good idea.

``` r
kruskal.test(copynumber ~ length, data = mtecutoff)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  copynumber by length
    ## Kruskal-Wallis chi-squared = 398209, df = 685, p-value < 2.2e-16

``` r
kruskal.test(copynumber ~ length, data = ftecutoff)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  copynumber by length
    ## Kruskal-Wallis chi-squared = 197803, df = 685, p-value < 2.2e-16

No significance is observed when looking at variation on the individual
level, which is probably not very surprising.

``` r
kruskal.test(copynumber ~ ID, data = mtecutoff)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  copynumber by ID
    ## Kruskal-Wallis chi-squared = 65.487, df = 552, p-value = 1

``` r
kruskal.test(copynumber ~ ID, data = ftecutoff)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  copynumber by ID
    ## Kruskal-Wallis chi-squared = 46.421, df = 274, p-value = 1

Lastly, we could look at the batch information. This information should
now be meaningless regarding technical differences. However, it seems
like it still has a significant effect. I believe this is solely caused
by the fact that very often whole populations are in one or the other
batch and it is not a random distribution.

``` r
kruskal.test(copynumber ~ batch, data = mtecutoff)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  copynumber by batch
    ## Kruskal-Wallis chi-squared = 0.049687, df = 1, p-value = 0.8236

``` r
kruskal.test(copynumber ~ batch, data = ftecutoff)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  copynumber by batch
    ## Kruskal-Wallis chi-squared = 1.9675, df = 1, p-value = 0.1607
