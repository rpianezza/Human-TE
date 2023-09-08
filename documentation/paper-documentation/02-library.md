HGDP - Library insights
================

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.1     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ ggplot2   3.4.2     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.1     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(ggpubr)

theme_set(theme_bw())

HGDP <- read_delim("/Volumes/Temp1/rpianezza/0.old/summary-HGDP/HGDP_cutoff_classified.tsv")
```

    ## Rows: 1394352 Columns: 12
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (9): ID, pop, sex, country, type, familyname, batch, superfamily, shared...
    ## dbl (3): length, reads, copynumber
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
DNA_names <- c("Crypton", "hAT", "Helitron", "Kolobok", "Mariner/Tc1", "Merlin", "MuDR", "piggyBac", "DNA transposon")
LINE_names <- c("L1", "CR1", "L2", "Crack", "RTE", "RTEX", "R4", "Vingi", "Tx1", "Penelope")
SINE_names <- c("SINE1/7SL", "SINE2/tRNA", "SINE3/5S", "SINE")
LTR_names <- c("ERV1", "ERV2", "ERV3", "Gypsy", "Endogenous Retrovirus", "LTR Retrotransposon", "Long terminal repeat", "Non-LTR Retrotransposon")
satellites_names <- c("Satellite", "satellite", "SAT")

classification <- HGDP %>% mutate(class = case_when(superfamily %in% DNA_names ~ "DNA", superfamily %in% LINE_names ~ "LINE", superfamily %in% SINE_names ~ "SINE", superfamily %in% LTR_names ~ "LTR", superfamily %in% satellites_names ~ "Satellite"))
```

``` r
(class <- classification %>% filter(type == "te") %>% select(familyname, class) %>% distinct() %>% group_by(class) %>% summarise(count = n()))
```

    ## # A tibble: 6 × 2
    ##   class     count
    ##   <chr>     <int>
    ## 1 DNA         219
    ## 2 LINE        137
    ## 3 LTR         511
    ## 4 SINE         12
    ## 5 Satellite    22
    ## 6 <NA>         64

``` r
(pie <- ggplot(class, aes(x="", y=count, fill=class))+
  geom_col(width = 1) + labs(fill='') +
  coord_polar(theta="y")+
  theme_void()+
  #geom_text(aes(label = count), 
           #position = position_stack(vjust = 0.5), 
           # size = 3)+
  theme(legend.position = "bottom")+
 scale_fill_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", "#D89000")))
```

![](02-library_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
#ggsave("/Volumes/Temp1/rpianezza/paper/figures/approach/library.png", pie, dpi=600)
```

``` r
(scg <- classification %>% filter(type=="scg") %>% group_by(familyname) %>% summarise(copynumber = mean(copynumber)))
```

    ## # A tibble: 483 × 2
    ##    familyname                copynumber
    ##    <chr>                          <dbl>
    ##  1 chr10:100746263-100749916      0.945
    ##  2 chr10:100987794-100989884      1.02 
    ##  3 chr10:101034691-101037167      1.00 
    ##  4 chr10:101226194-101229794      0.979
    ##  5 chr10:102065389-102068038      1.03 
    ##  6 chr10:102394109-102396298      0.963
    ##  7 chr10:11316832-11319884        1.02 
    ##  8 chr10:116670033-116672635      0.948
    ##  9 chr10:117239599-117241923      0.883
    ## 10 chr10:120899837-120903175      0.978
    ## # ℹ 473 more rows

``` r
(scg_plot <- ggplot(scg, aes(familyname, copynumber))+
    geom_point() + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()))
```

![](02-library_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->
