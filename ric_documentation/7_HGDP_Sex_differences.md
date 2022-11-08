Comparing the mean TE copy numbers in males and females
================

In this script I re-wrote Florian’s script 5 and added some new
analysis. The topic is the comparison of TE abundances between the two
sexes.

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
library(dplyr)
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

## Calculating mean abundance for males and females

This chuck creates the dataset used for all the subsequent analyses.
Starting from the
*USEME_HGDP_complete_reflib6.2_mq10_batchinfo_cutoff0.01.txt* file, the
final dataset is composed by `965` rows, one for each TE family in the
dataset. For each TE family, we have:

- The **mean copynumber** among all the **females** in the dataset
  (`f_mean`), as well as in **males** (`m_mean`).
- The **log** of the two means (`f_mean_log`, `m_mean_log`).
- The **difference** between the copynumber in males and females
  (`diff`) and the absolute difference (`abs_diff`).
- The **ratio** between the copynumbers in the two sexes (`ratio`),
  calculated with the more abundant sex always at the numerator.
- Last, `more_in` indicates in which sex the TE is more abundant.

``` r
data <- filter(HGDPcutoff, type=="te") %>% group_by(familyname, sex) %>% summarise(mean=mean(copynumber))
```

    ## `summarise()` has grouped output by 'familyname'. You can override using the
    ## `.groups` argument.

``` r
f<-filter(data, sex=="female") %>% rename(f_mean = mean)
m<-filter(data, sex=="male") %>% rename(m_mean = mean)

(final_data <- inner_join(f, m, by = "familyname") %>% select(familyname, f_mean, m_mean) %>% mutate(f_mean_log=log(f_mean), m_mean_log=log(m_mean), diff=m_mean-f_mean, abs_diff=abs(diff), ratio=case_when(diff>=0 ~ m_mean/f_mean, diff<0 ~ f_mean/m_mean), more_in=case_when(diff>=0 ~ "male", diff<0 ~ "female")) %>% arrange(desc(ratio)) %>% mutate(familyname=fct_reorder(familyname,ratio)))
```

    ## # A tibble: 965 × 9
    ## # Groups:   familyname [965]
    ##    familyname     f_mean   m_mean f_mea…¹ m_mea…²     diff abs_d…³ ratio more_in
    ##    <fct>           <dbl>    <dbl>   <dbl>   <dbl>    <dbl>   <dbl> <dbl> <chr>  
    ##  1 HSATI      334.        2.00e+3   5.81    7.60   1.67e+3 1.67e+3  5.99 male   
    ##  2 X17_DNA      0.708     3.68e-1  -0.346  -1.00  -3.40e-1 3.40e-1  1.92 female 
    ##  3 UCON64       0.775     4.08e-1  -0.255  -0.897 -3.67e-1 3.67e-1  1.90 female 
    ##  4 Tigger16b    0.206     1.09e-1  -1.58   -2.22  -9.74e-2 9.74e-2  1.90 female 
    ##  5 X3_LINE      0.356     1.88e-1  -1.03   -1.67  -1.68e-1 1.68e-1  1.89 female 
    ##  6 UCON106      0.368     1.97e-1  -1.00   -1.62  -1.71e-1 1.71e-1  1.87 female 
    ##  7 X15_DNA      0.000994  5.47e-4  -6.91   -7.51  -4.47e-4 4.47e-4  1.82 female 
    ##  8 UCON55       0.000261  1.45e-4  -8.25   -8.84  -1.17e-4 1.17e-4  1.81 female 
    ##  9 UCON79       0.000387  6.95e-4  -7.86   -7.27   3.08e-4 3.08e-4  1.79 male   
    ## 10 CHARLIE4     0.718     4.08e-1  -0.331  -0.897 -3.11e-1 3.11e-1  1.76 female 
    ## # … with 955 more rows, and abbreviated variable names ¹​f_mean_log,
    ## #   ²​m_mean_log, ³​abs_diff

## Comparing mean abundance between males and females

We want to create a plot that shows a pairwise comparison of the mean
copy number estimates of males and females for each of the 965 TE
sequences. Additionally, we include a linear regression line (in grey)
as well as a dashed line representing the null expectation for the
regression, i.e. that males and females have the same mean abundance
values and the regression should thus have a slope of 1. The regression
line having a slope slightly lower than one is thus an indication
(though not a proof) that TEs in females (x-axis) might be on average
slightly more abundant than TEs in males (y-axis).

``` r
ggplot(final_data, aes(x=f_mean, y=m_mean, color=more_in))+ labs(color = "More abundant in:") +
  geom_point(size=1)+
  geom_smooth(method="lm",color="grey", se=F)+ ylab("Male mean")+ xlab("Female mean")+geom_abline(slope=1,linetype='dashed')+theme_bw()
```

    ## `geom_smooth()` using formula 'y ~ x'

![](7_HGDP_Sex_differences_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
ggplot(final_data, aes(f_mean_log, m_mean_log, color=more_in))+ labs(color = "More abundant in:") +
  geom_point(size=0.8)+
  geom_smooth(method="lm",color="grey",se=F)+theme_bw()+ ylab("Male mean (log)")+ xlab("Female mean (log)")
```

    ## `geom_smooth()` using formula 'y ~ x'

![](7_HGDP_Sex_differences_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

## Relative and absolute mean differences between the sexes

### Absolute differences

The following plots shows, respectively:

- The differences between sex abundances of all the 965 TEs.
- A subset of the first plot, with only TEs with `diff > 75`.

``` r
abs_subset<- filter(final_data, abs_diff>75) # This number is arbitrary, feel free to look at more/less TEs in the plot

ggplot(abs_subset, aes(reorder(familyname, -abs(diff)), abs_diff, fill=more_in)) + labs(fill = "More abundant in:") +
  geom_bar(stat="identity") + ylab("Difference between sex abundances") + xlab("Repetitive sequence families") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

![](7_HGDP_Sex_differences_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

### Relative differences

The following plots shows, respectively:

- The ratio between sex abundances of all the 965 TEs.
- A subset of the first plot, with only TEs with `ratio > 1.25`.

``` r
ggplot(final_data, aes(reorder(familyname, -ratio), ratio, fill=more_in)) + labs(fill = "More abundant in:") +
  geom_bar(stat="identity") + ylab("Ratio between sex abundances") + xlab("Repetitive sequence families") + theme(axis.text.x=element_blank())
```

![](7_HGDP_Sex_differences_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
rel_subset<- filter(final_data, ratio>1.25) # This number is arbitrary, feel free to look at more/less TEs in the plot

ggplot(rel_subset, aes(reorder(familyname, -ratio), ratio, fill=more_in)) + labs(fill = "More abundant in:") +
  geom_bar(stat="identity") + ylab("Ratio between sex abundances") + xlab("Repetitive sequence families") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

![](7_HGDP_Sex_differences_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

## L1 family differences between the two sexes

One of the most interesting TE groups found so far in the dataset are
for sure the `LINE-1` (**L1**) retrotransposon subfamilies. Here, I show
that almost every L1 is more abundant in the `females`. Only 2/113 L1
subfamilies are more abundant in `males`, and their ratio is very close
to 1. On the other hand, we have lot of L1 subfamilies with higher
copynumber in females, and for some of them the difference is very
consistent.

``` r
L1<-final_data %>% filter(str_detect(familyname, "^L1"))
ggplot(L1, aes(reorder(familyname, -abs(diff)), abs_diff, fill=more_in)) + labs(fill = "More abundant in:") +
  geom_bar(stat="identity") + ylab("Difference between sex abundances") + xlab("Repetitive sequence families") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=4))
```

![](7_HGDP_Sex_differences_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
(L1_male <- filter(L1, more_in == "male"))
```

    ## # A tibble: 2 × 9
    ## # Groups:   familyname [2]
    ##   familyname f_mean m_mean f_mean_log m_mean_log    diff abs_diff ratio more_in
    ##   <fct>       <dbl>  <dbl>      <dbl>      <dbl>   <dbl>    <dbl> <dbl> <chr>  
    ## 1 L1M2A_5    22.6   22.9         3.12       3.13 0.313    0.313    1.01 male   
    ## 2 L1M6B_5end  0.127  0.129      -2.06      -2.05 0.00151  0.00151  1.01 male

``` r
nrow(L1_male)
```

    ## [1] 2

``` r
nrow(L1)
```

    ## [1] 113

This plot suggests that the LINE-1s are in general more abundant on the
**X chromosome** than on the **Y**.

``` r
L1_rel <- filter(L1, ratio>1.08)

ggplot(L1_rel, aes(reorder(familyname, -abs(ratio)), ratio, fill=more_in)) + labs(fill = "More abundant in:") +
  geom_bar(stat="identity") + ylab("Difference between sex abundances") + xlab("Repetitive sequence families")
```

![](7_HGDP_Sex_differences_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

If we filter for the L1’s with the highest `ratio` between the two sexes
copynumber abundances, we find that the most relevant TE is again
`L1ME5`.
