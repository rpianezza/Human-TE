HGDP - Mining for other interesting TEs
================

This is the third script written by me, Riccardo. This script works with
the version of the HGDP dataset created in Script 2 by Florian. We
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
library("dplyr")
library("reticulate")

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

HGDP_<-read_delim("/Users/rpianezza/TE/summary-HGDP/nocutoff/HGDP_complete_reflib6.2_mq10_batchinfo.txt",comment="#")
```

    ## Rows: 1410083 Columns: 10
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (7): HGDP00001, Brahui, male, Central_South_Asia, scg, chr1:916864-92101...
    ## dbl (3): 4152, 1051.8410596026488, 1.01525454783029
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
names(HGDP_)<-c("ID","Pop","sex","Country","type","familyname","length","reads","copynumber","batch")
HGDP <- HGDP_ %>% mutate(Country = replace(Country, Country=="Oceania_(SGDP),Oceania", "Oceania"))
```

# Mining the dataset looking for other interesting TEs

The idea of this script is to play a bit with the thresholds used in the
firsts two scripts, looking for other TEs with interesting
distributions.

## List of previosly analysed TEs

First, I create a list containing all the `familyname` of the TEs
previously taken into account in my analyses. In total, I already found
and analyzed `58` TE.

``` r
TE_cutoff <- filter(HGDPcutoff, type=='te') %>% group_by(familyname, sex) %>% mutate(max=max(copynumber), min=min(copynumber)) %>% mutate(diff = max-min, ratio = max/min)

out1_abs <- filter(TE_cutoff, diff>200 & diff<Inf)
out1_rel <- filter(TE_cutoff, ratio>2 & ratio<Inf & max>1.5)

out1_abs_names <- c(unique(out1_abs$familyname))
out1_rel_names <- c(unique(out1_rel$familyname))
(out1_names <- c(out1_abs_names, out1_rel_names[!(out1_rel_names %in% out1_abs_names)]) %>% sort())
```

    ##  [1] "6kbHsap"     "ALR"         "ALR_"        "ALR1"        "ALR2"       
    ##  [6] "ALRa_"       "ALRb"        "ALU"         "AmnSINE1_HS" "CER"        
    ## [11] "EuthAT-N1"   "Eutr11"      "Eutr8"       "EUTREP13"    "GSATII"     
    ## [16] "HERVI"       "HSATI"       "HSATII"      "KER"         "L1"         
    ## [21] "L1HS"        "L1MC1"       "L1ME5"       "L1P_MA2"     "L1PA10"     
    ## [26] "L1PA15"      "L1PA16"      "L1PA3"       "L1PA4"       "L1PA6"      
    ## [31] "L1PA7"       "L1PA7_5"     "L1PA8"       "L1PB1"       "L1PB2"      
    ## [36] "L1PB2c"      "L1PB4"       "L1PREC1"     "L1PREC2"     "LSAU"       
    ## [41] "MARE8"       "MER2"        "MER22"       "MER69A"      "MLT1B"      
    ## [46] "MLT2A1"      "MLT2A2"      "MSTA"        "SVA_A"       "TAR1"       
    ## [51] "THE1_I"      "THE1A"       "THE1B"       "THE1C"       "THE1D"      
    ## [56] "TIGGER1"     "UCON54"      "UCON75"

## Playing with the thresholds

To not miss any possible interesting TE, here I use the file without the
**cutoff** at 0.1 copynumber for TEs, the one we used until now.

``` r
TE <- filter(HGDP, type=='te') %>% group_by(familyname, sex) %>% mutate(max=max(copynumber), min=min(copynumber), mean=mean(copynumber)) %>% mutate(diff = max-min, ratio = max/min)
```

``` r
plotTEfamily <- function(data, Sex, famname, binwidht, x_title, y_title, x_numbers, y_numbers){
filtered <- filter(data, familyname==famname, sex==Sex)
ggplot(data = filtered, mapping = aes(x = copynumber, fill = Country)) +
  geom_histogram(binwidth = binwidht) +
  ggtitle(paste0(Sex,' - ',famname)) + theme(plot.title = element_text(size = 8, hjust = 0.5)) +
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

### Filter 1

``` r
all1 <- filter(TE, diff<200 & diff>150)
new1 <- filter(all1, !(familyname %in% out1_names))
(list1 <- unique(new1$familyname))
```

    ##  [1] "MER4A1" "L1MA4"  "LTR7A"  "L1PA13" "LTR12C" "L1MA1"  "L1MA3"  "ALRa"  
    ##  [9] "MER1B"  "MER9"

``` r
out2_names <- c(out1_names, list1) %>% sort()
```

Filtering for `diff<200` & `diff>150`, we exclude the TEs previously
found in script 1 but still with a significant variance. We found
`MER4A1`, `L1MA4`, `LTR7A`, `L1PA13`, `LTR12C`, `L1MA1`, `L1MA3`,
`ALRa`, `MER1B`, `MER9`. The LINEs among this list shows a pattern super
consistent to the one observed in previous scripts. The other sequences
do not look of particular interest.

``` r
plotTEfamily(TE, 'male', "L1MA4", 1, 'y', 'y', 'y', 'y') 
```

![](3_HGDP_Other_TEs_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
plotTEfamily(TE, 'female', "L1MA4", 1, 'y', 'y', 'y', 'y')
```

![](3_HGDP_Other_TEs_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
plotTEfamily(TE, 'male', "L1PA13", 1, 'y', 'y', 'y', 'y') 
```

![](3_HGDP_Other_TEs_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

``` r
plotTEfamily(TE, 'female', "L1PA13", 1, 'y', 'y', 'y', 'y')
```

![](3_HGDP_Other_TEs_files/figure-gfm/unnamed-chunk-6-4.png)<!-- -->

``` r
plotTEfamily(TE, 'male', "L1MA1", 1, 'y', 'y', 'y', 'y') 
```

![](3_HGDP_Other_TEs_files/figure-gfm/unnamed-chunk-6-5.png)<!-- -->

``` r
plotTEfamily(TE, 'female', "L1MA1", 1, 'y', 'y', 'y', 'y')
```

![](3_HGDP_Other_TEs_files/figure-gfm/unnamed-chunk-6-6.png)<!-- -->

``` r
plotTEfamily(TE, 'male', "L1MA3", 1, 'y', 'y', 'y', 'y') 
```

![](3_HGDP_Other_TEs_files/figure-gfm/unnamed-chunk-6-7.png)<!-- -->

``` r
plotTEfamily(TE, 'female', "L1MA3", 1, 'y', 'y', 'y', 'y')
```

![](3_HGDP_Other_TEs_files/figure-gfm/unnamed-chunk-6-8.png)<!-- --> It
may be interesting a comprehensive analysis of all **L1s**?

### Filter 2

This is a complicated filter, but it could be decisive. The idea is to
find the TEs in which **copynumber distributions** are present the
bigger **gaps**. The rationale is that a distribution with a significant
gap could represent an invasion of the TE among some populations, while
being in low copynumber in others. For example, in the 2nd script we
found the LINE1 variant `L1ME5` which was behaving exactly like that.

First, I order the dataset for `familyname`, `sex` and `copynumber`.
Then, I create a *CSV* file with the created dataset. The reason is that
I did not know how to perform the next coding steps in R, so I did that
in *Python 3*, in the script **find-gaps** which can be found in the
“ric-documentation/other-files/” folder.

``` r
(ordered <- arrange(TE, familyname, sex, copynumber))
```

    ## # A tibble: 814,752 × 15
    ## # Groups:   familyname, sex [1,968]
    ##    ID        Pop   sex   Country type  famil…¹ length  reads copyn…² batch   max
    ##    <chr>     <chr> <chr> <chr>   <chr> <chr>    <dbl>  <dbl>   <dbl> <chr> <dbl>
    ##  1 HGDP01373 Basq… fema… Europe  te    6kbHsap   6018 2.33e5    178. flo    429.
    ##  2 HGDP01366 Basq… fema… Europe  te    6kbHsap   6018 3.42e5    196. flo    429.
    ##  3 HGDP00323 Kala… fema… Centra… te    6kbHsap   6018 4.74e5    201. flo    429.
    ##  4 HGDP01276 Moza… fema… Middle… te    6kbHsap   6018 4.64e5    212. flo    429.
    ##  5 HGDP00797 Orca… fema… Europe  te    6kbHsap   6018 3.19e5    214. ro     429.
    ##  6 HGDP00581 Druze fema… Middle… te    6kbHsap   6018 3.20e5    218. flo    429.
    ##  7 HGDP01280 Moza… fema… Middle… te    6kbHsap   6018 4.16e5    221. flo    429.
    ##  8 HGDP00872 Maya  fema… America te    6kbHsap   6018 4.20e5    223. flo    429.
    ##  9 HGDP00244 Path… fema… Centra… te    6kbHsap   6018 3.86e5    227. ro     429.
    ## 10 HGDP00738 Pale… fema… Middle… te    6kbHsap   6018 3.57e5    228. ro     429.
    ## # … with 814,742 more rows, 4 more variables: min <dbl>, mean <dbl>,
    ## #   diff <dbl>, ratio <dbl>, and abbreviated variable names ¹​familyname,
    ## #   ²​copynumber

``` r
#write.csv(ordered, "/Users/rpianezza/TE/ric-documentation-Rmd/ordered")
```

The Python script allowed me to find the biggest gap in each TE
copynumber distribution, divided by sex. In the next code chunk I work a
bit on the dataset structure, resulting in the `data` dataset. It
contains two rows for each TE family (one for each sex). In every row is
reported the biggest gap in the copynumber distribution (`maxgap`) as
well as the ratio between the `maxgap` and the `mean`
(`gap_mean_ratio`).

I added the `gap_mean_ratio` to have a measure of the entity of the gap.
A TE with a mean copynumber of 1000 and a gap of 2 is not interesting,
while a TE with a mean of 5 and a gap of 2 **IS** interesting! Thus, I
ordered `data` for `gap_mean_ratio` in descendant order, so that the TEs
with higher ratio are showed first.

Lastly, I filtered for TEs with a `max > 1`, to remove all these TEs
with superlow copynumber.

``` r
gaps <- read_tsv("/Users/rpianezza/TE/ric-documentation-Rmd/gaps.csv")
```

    ## Rows: 1967 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): familyname, sex
    ## dbl (1): maxgap
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
gaps_sorted <- arrange(gaps, sex, desc(maxgap))
mean <- select(TE, familyname, sex, min, max, mean)
(data <- inner_join(x=gaps_sorted, y=mean, by=c('familyname','sex')) %>% unique() %>% mutate(gap_mean_ratio=maxgap/mean) %>% relocate(maxgap, .after = mean) %>% arrange(desc(gap_mean_ratio)) %>% filter(max>1))
```

    ## # A tibble: 1,507 × 7
    ##    familyname  sex        min     max    mean   maxgap gap_mean_ratio
    ##    <chr>       <chr>    <dbl>   <dbl>   <dbl>    <dbl>          <dbl>
    ##  1 HSATI       female 142.    2297.   334.    1716.             5.14 
    ##  2 AmnSINE1_HS male     0.242    2.71   0.447    1.83           4.10 
    ##  3 L1ME5       female   0.267    8.62   6.62     5.45           0.823
    ##  4 L1ME5       male     0.201    8.01   4.94     2.80           0.568
    ##  5 UCON53      female   0.312    1.16   0.566    0.220          0.388
    ##  6 MER22       male    52.6    364.   165.      56.5            0.341
    ##  7 UCON93      male     0.211    1.15   0.566    0.183          0.324
    ##  8 UCON46      female   0.433    1.32   0.750    0.241          0.321
    ##  9 UCON96      male     0.224    1.01   0.539    0.168          0.311
    ## 10 HERV23      male     0.690    1.33   0.865    0.264          0.306
    ## # … with 1,497 more rows

``` r
#all2 <- filter(data, gap_mean_ratio>0.3)
#new2 <- filter(all2, !(familyname %in% out2_names))
#(list2 <- unique(new2$familyname))
```

Starting from the top of the dataset, we have `HSATI` and `AmnSINE1_HS`,
already analysed in script 1/2. Their `gap_mean_ratio` is an outlier
among the rest of the TEs, but looking at their distribution we notice
that there is only 1 individual, probably an error, with high copynumber
that generate this high values.

The third TE, but the first **significant** TE, is the well known
`L1ME5`, which was analysed in script 2 and which was considered the
most interesting TE. Here it not disappoint us. Both for males and
females, the `gap_mean_ratio` is significantly higher than for the other
TEs.

I plotted the other TEs with high `gap_mean_ratio`, but none of them
showed a pattern similar to `L1ME5` or other interesting patterns.

I conclude that `L1ME5` remains the strongest candidate to be an active
TE in the human genome which invaded **non-africans** populations.
