HGDP - Analyzing the TE abundance in the two sexes
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
a_HGDP <- read_tsv("/Volumes/Temp1/rpianezza/PCA-copynumber-all-analysis/a_HGDP.tsv")
```

    ## Rows: 828 Columns: 2
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): ID
    ## dbl (1): a
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
HGDP_nobiased_samples <- filter(a_HGDP, (a > (-0.5)) & (a<0.5)) %>% select(ID) %>% pull()
HGDP_clean <- filter(HGDP, ID %in% HGDP_nobiased_samples)

DNA_names <- c("Crypton", "hAT", "Helitron", "Kolobok", "Mariner/Tc1", "Merlin", "MuDR", "piggyBac", "DNA transposon")
LINE_names <- c("L1", "CR1", "L2", "Crack", "RTE", "RTEX", "R4", "Vingi", "Tx1", "Penelope")
SINE_names <- c("SINE1/7SL", "SINE2/tRNA", "SINE3/5S", "SINE")
LTR_names <- c("ERV1", "ERV2", "ERV3", "Gypsy", "Endogenous Retrovirus", "LTR Retrotransposon", "Long terminal repeat", "Non-LTR Retrotransposon")
satellites_names <- c("Satellite", "satellite", "SAT")

classification <- HGDP_clean %>% mutate(class = case_when(superfamily %in% DNA_names ~ "DNA", superfamily %in% LINE_names ~ "LINE", superfamily %in% SINE_names ~ "SINE", superfamily %in% LTR_names ~ "LTR", superfamily %in% satellites_names ~ "satellite"))
```

# CV

``` r
cn <- classification %>% filter(type=="te") %>% group_by(familyname, class, superfamily) %>% summarise(mean = mean(copynumber), var = var(copynumber), cv = var/mean) %>% arrange(desc(cv))
```

    ## `summarise()` has grouped output by 'familyname', 'class'. You can override
    ## using the `.groups` argument.

``` r
no_sat <- cn %>% filter(class!="satellite", class!="NA")

top25 <- cn %>% filter(!(class %in% c("satellite", NA)), mean>1)%>% filter(cv > 2.2)

outliers <- HGDP_clean %>% filter(type=="te", familyname %in% top25$familyname) %>% inner_join(cn, by="familyname")
outliers <- outliers[order(outliers$cv,decreasing=T),]
outliers$familyname<-factor(outliers$familyname,levels=unique(outliers$familyname))

(box <- ggplot(outliers, aes(x=familyname, y=copynumber)) + geom_boxplot(notch=F, aes(color=class)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylim(0,7000)+ facet_grid(~sex) +
 scale_color_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")))
```

    ## Warning: Removed 660 rows containing non-finite values (`stat_boxplot()`).

![](cv-rethinking_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
active <- c("ALU","SVA_A","L1HS","HERV-K14CI","L1PA2")

rank <- no_sat %>% ungroup() %>% arrange(desc(mean)) %>% mutate(cn_rank = row_number()) %>% arrange(desc(cv)) %>% mutate(cv_rank = row_number())

(active <- rank %>% filter(familyname%in%active) %>% select(familyname, class, cn_rank, cv_rank))
```

    ## # A tibble: 5 × 4
    ##   familyname class cn_rank cv_rank
    ##   <chr>      <chr>   <int>   <int>
    ## 1 ALU        SINE        1       1
    ## 2 L1HS       LINE       17      15
    ## 3 SVA_A      SINE       14      48
    ## 4 L1PA2      LINE       80      53
    ## 5 HERV-K14CI LTR       482     183

``` r
(rank_diff <- rank %>% mutate(rank_increase = cn_rank-cv_rank) %>% filter(mean>10) %>% arrange(desc(rank_increase)))
```

    ## # A tibble: 482 × 9
    ##    familyname class superfamily            mean    var     cv cn_rank cv_rank
    ##    <chr>      <chr> <chr>                 <dbl>  <dbl>  <dbl>   <int>   <int>
    ##  1 HERV-K14CI LTR   ERV2                   10.0  1.38  0.137      482     183
    ##  2 LTR25      LTR   ERV1                   19.5  3.43  0.176      401     153
    ##  3 LTR21A     LTR   ERV1                   19.2  1.50  0.0781     403     261
    ##  4 MER4BI     LTR   Endogenous Retrovirus  31.2  4.12  0.132      329     192
    ##  5 MER6C      DNA   Mariner/Tc1            11.4  0.505 0.0444     470     342
    ##  6 L1MDB_5    LINE  L1                     16.3  0.914 0.0560     429     302
    ##  7 LTR22      LTR   ERV2                   51.2 12.9   0.253      254     128
    ##  8 HERV30I    LTR   Endogenous Retrovirus  10.4  0.435 0.0416     477     352
    ##  9 HERVFH21I  LTR   Endogenous Retrovirus  14.4  0.758 0.0525     440     316
    ## 10 LTR06      LTR   ERV1                   12.4  0.516 0.0418     460     350
    ## # ℹ 472 more rows
    ## # ℹ 1 more variable: rank_increase <int>

``` r
#write_csv(table_active, "/Volumes/Temp1/rpianezza/paper/tables/active.csv")
```

``` r
plot_cn_country <- function(data, fam){
  
cn <- data %>% filter(familyname==fam)

(plots <- ggplot(cn, aes(x = familyname, y = copynumber, color = country)) + geom_boxplot(notch = FALSE, width = 0.8, lwd = 0.2, outlier.size = 1) + facet_grid(~sex) + 
theme(strip.background = element_blank(), strip.text.x = element_blank()) + xlab(NULL))
}

alu <- classification %>% filter(familyname=="ALU")
(f1 <- ggplot(alu, aes(x = familyname, y = copynumber, color = country)) + geom_boxplot(notch = FALSE, width = 0.8, lwd = 0.2, outlier.size = 1) + facet_grid(~sex) + xlab(NULL))
```

![](cv-rethinking_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
f2 <- plot_cn_country(classification, "L1PA4")
f3 <- plot_cn_country(classification, "THE1B")
f4 <- plot_cn_country(classification, "L1")
f5 <- plot_cn_country(classification, "THE1_I")
f6 <- plot_cn_country(classification, "L1PREC1")
f7 <- plot_cn_country(classification, "THE1C")
f8 <- plot_cn_country(classification, "L1PA16")
f9 <- plot_cn_country(classification, "THE1D")
f10 <- plot_cn_country(classification, "L1PA6")

f11 <- plot_cn_country(classification, "THE1A")
f12 <- plot_cn_country(classification, "MSTA")
f13 <- plot_cn_country(classification, "L1PB2c")
f14 <- plot_cn_country(classification, "L1PA7")
f15 <- plot_cn_country(classification, "L1HS")
f16 <- plot_cn_country(classification, "L1PB1")
f17 <- plot_cn_country(classification, "L1PREC2")
f18 <- plot_cn_country(classification, "TIGGER1")
f19 <- plot_cn_country(classification, "MLT2A2")
f20 <- plot_cn_country(classification, "L1P_MA2")

f21 <- plot_cn_country(classification, "L1PA3")
f22 <- plot_cn_country(classification, "L1PA7_5")
f23 <- plot_cn_country(classification, "L1PA10")
f24 <- plot_cn_country(classification, "L1MA1")
(f25 <- plot_cn_country(classification, "L1PA15"))
```

![](cv-rethinking_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
cv_all <- ggarrange(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,f19,f20,f21,f22,f23,f24,f25, ncol = 1, common.legend = TRUE)
cv_10 <- ggarrange(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10, ncol = 1, common.legend = TRUE)


#ggsave(cv_all, file="/Volumes/Temp1/rpianezza/paper/cv_all.png", height = 25, width = 12, dpi=1000)
#ggsave(cv_10, file="/Volumes/Temp1/rpianezza/paper/cv_10.png", height = 15, width = 12, dpi=1000)
```

# Other stuff

``` r
library(scales)
```

    ## 
    ## Attaching package: 'scales'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     discard

    ## The following object is masked from 'package:readr':
    ## 
    ##     col_factor

``` r
#extract hex color codes for a plot with three elements in ggplot2 
hex <- hue_pal()(10)

#overlay hex color codes on actual colors
show_col(hex)
```

![](cv-rethinking_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
