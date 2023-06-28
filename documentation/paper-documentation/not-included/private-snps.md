HGDP - Detecting private SNPs of single TEs
================

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.0     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ ggplot2   3.4.1     ✔ tibble    3.1.8
    ## ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.1     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the ]8;;http://conflicted.r-lib.org/conflicted package]8;; to force all conflicts to become errors

``` r
library(ggpubr)
library(umap)

theme_set(theme_bw())

HGDP <- read_delim("/home/riccardo/Desktop/human-te/HGDP_cutoff_classified.tsv")
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

classification <- HGDP %>% mutate(class = case_when(superfamily %in% DNA_names ~ "DNA", superfamily %in% LINE_names ~ "LINE", superfamily %in% SINE_names ~ "SINE", superfamily %in% LTR_names ~ "LTR", superfamily %in% satellites_names ~ "satellite"))

HGDP_meta <- HGDP %>% select(ID, pop, country, sex) %>% distinct()

females <- HGDP %>% filter(sex=="female") %>% select(ID)
males <- HGDP %>% filter(sex=="male") %>% select(ID)
```

## Functions

Function which finds private SNPs in each population. Takes as input a
merged file (PCAable + metadata) and a frequency threshold (es. 0.5). It
counts how many SNPs are present in at least 50% of the individuals in a
population and not present in more than 15 other samples.

``` r
private <- function(data, threshold, percentage_pop_with_snp){
snps <- colnames(data)[3:length(colnames(data))]
populations <- data %>% select(location) %>% distinct() %>% pull()
private_snps <- tibble(location=NULL, SNP=NULL, thr=NULL, perc_shared_in_pop=NULL, number_external=NULL)

for (pop in populations){
  for (snp in snps){
    p <- data %>% select(ID, location, !!snp) %>% filter(location==pop)
    other <- filter(data, location!=pop)
    shared <- sum(p[[snp]] < threshold)
    count_pop <- p %>% summarise(count = n()) %>% pull()
    other_shared <- sum(other[[snp]] < threshold)
    if ((shared > (count_pop/100*percentage_pop_with_snp)) & other_shared<15){
      row <- tibble(location=pop, SNP=snp, thr=threshold, perc_shared_in_pop=shared/count_pop, number_external=other_shared)
      private_snps <- bind_rows(private_snps, row)
    }
  }
}
private_snps
}
```

Function to loop the function “private()” over multiple thresholds (from
0.5 to 0.9). A private SNP with 0.5 frequency is a very strong proof of
a bottleneck during the invasion, while a SNP with only 0.9 is not that
significant.

``` r
different_thresholds <- function(data, percentage_pop_with_snp){
result_tibble <- tibble(location=NULL, SNP=NULL, thr=NULL, perc_shared_in_pop=NULL, number_external=NULL)
for (t in seq(0.5, 0.95, +0.05)){
  tib <- private(data, t, percentage_pop_with_snp)
  result_tibble <- bind_rows(result_tibble, tib)
}

if (!is_empty(result_tibble)){
pops <- result_tibble %>% select(location) %>% distinct() %>% pull()
filtered_tibble <- tibble(location=NULL, SNP=NULL, thr=NULL, perc_shared_in_pop=NULL, number_external=NULL)
for (pop in pops){
  only_pop <- filter(result_tibble, location==pop)
  filtered <- only_pop %>% filter(thr==0.5)
for (t in seq(0.5, 0.95, +0.05)){
  th <- filter(only_pop, thr<t) %>% select(SNP) %>% pull()
  non_overlapping <- filter(only_pop, thr==t, !(SNP %in% th))
  filtered <- bind_rows(filtered, non_overlapping)
}
  filtered_tibble <- bind_rows(filtered_tibble, filtered) %>% distinct()
}
}
else{
  filtered_tibble <- tibble(location=NULL, SNP=NULL, thr=NULL, perc_shared_in_pop=NULL, number_external=NULL)
}
filtered_tibble
}
```

Function which takes as input the tibble returned from
“different_thresholds” and creates a barplot.

``` r
private_plot <- function(data, metadata, titlee){
  all_countries <- metadata %>% select(country) %>% distinct() %>% pull()
  all_thresholds <- seq(0.05, 0.5, +0.05)
  
  if (!is_empty(data)){
  plottable <- data %>% group_by(location, thr) %>% summarise(count = n()) %>% mutate(thr=1-thr)
  countries <- data %>% select(location) %>% distinct() %>% pull()
  thresholds <- data %>% select(thr) %>% distinct() %>% pull()
  missing_countries <- setdiff(all_countries, countries)
  missing_thresholds <- setdiff(all_thresholds, thresholds)
  for (t in missing_thresholds){
  missing <- tibble(location = missing_countries, thr = t, count = 0)
  plottable <- bind_rows(plottable, missing)}
  plottable$thr <- round((plottable$thr), 2)
  plottable$thr <- as.character(plottable$thr)
  }
  
  else{
    plottable <- tibble(location=NULL, thr=NULL, count=NULL)
    for (t in all_thresholds){
  missing <- tibble(location = all_countries, thr = t, count = 0)
  plottable <- bind_rows(plottable, missing)}
  plottable$thr <- round((plottable$thr), 2)
  plottable$thr <- as.character(plottable$thr)
  }
  
  ggplot(plottable, aes(x = thr, y = count, fill = location)) +
  geom_bar(stat = "identity") +
  labs(x = "Frequency threshold", y = "Private SNPs", fill = "Location") +
  ggtitle(titlee) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) #+ xlim(0,0.6)
}
```

``` r
cn_filtered <- HGDP %>% group_by(familyname) %>% summarise(mean_cn = mean(copynumber)) %>% filter(mean_cn > 5) %>% mutate(te = paste0(familyname, "_te")) %>% select(te) %>% pull()

folder = "/home/riccardo/Desktop/human-te/single-TEs/pcaable/"
files <- list.files(path = folder, full.names = FALSE)
meta_country <- HGDP_meta %>% select(-pop, -sex) %>% rename("location"="country")

for (file in files){
  te <- sub("^[^_]*_", "", file)
  if (te %in% cn_filtered){
    path = paste0(folder, file)
    test <- read_tsv(path) %>% separate(familyname_position, into = c("ID", "pop"), sep = "-") %>% select(-pop)
    joined <- inner_join(meta_country, test, by="ID")
    tib_test <- different_thresholds(joined, 33)
    plo <- private_plot(tib_test, HGDP_meta, te)
    print(plo)
  }
}
```

    ## Rows: 275 Columns: 330
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (329): 6kbHsap_te_8, 6kbHsap_te_29, 6kbHsap_te_44, 6kbHsap_te_54, 6kbHsa...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.
    ## Rows: 275 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (22): ALR__te_15, ALR__te_29, ALR__te_30, ALR__te_31, ALR__te_37, ALR__t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

    ## Rows: 275 Columns: 16
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (15): ALR1_te_16, ALR1_te_19, ALR1_te_20, ALR1_te_22, ALR1_te_41, ALR1_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

    ## Rows: 275 Columns: 11
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (10): ALR2_te_63, ALR2_te_67, ALR2_te_101, ALR2_te_102, ALR2_te_119, ALR...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

    ## Rows: 275 Columns: 16
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (15): ALRa__te_13, ALRa__te_14, ALRa__te_21, ALRa__te_32, ALRa__te_43, A...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->

    ## Rows: 275 Columns: 18
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (17): ALRa_te_15, ALRa_te_30, ALRa_te_48, ALRa_te_65, ALRa_te_83, ALRa_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-5.png)<!-- -->

    ## Rows: 275 Columns: 21
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (20): ALRb_te_8, ALRb_te_15, ALRb_te_32, ALRb_te_37, ALRb_te_38, ALRb_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-6.png)<!-- -->

    ## Rows: 275 Columns: 51
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (50): ALU_te_9, ALU_te_11, ALU_te_20, ALU_te_57, ALU_te_62, ALU_te_63, A...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-7.png)<!-- -->

    ## Rows: 275 Columns: 44
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (43): BSRf_te_13, BSRf_te_20, BSRf_te_30, BSRf_te_34, BSRf_te_57, BSRf_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-8.png)<!-- -->

    ## Rows: 275 Columns: 21
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (20): CER_te_10, CER_te_17, CER_te_25, CER_te_58, CER_te_65, CER_te_67, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-9.png)<!-- -->

    ## Rows: 275 Columns: 369
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (368): CHARLIE1_te_108, CHARLIE1_te_114, CHARLIE1_te_115, CHARLIE1_te_12...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-10.png)<!-- -->

    ## Rows: 275 Columns: 186
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (185): Charlie12_te_17, Charlie12_te_18, Charlie12_te_19, Charlie12_te_2...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-11.png)<!-- -->

    ## Rows: 275 Columns: 211
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (210): CHARLIE1A_te_44, CHARLIE1A_te_48, CHARLIE1A_te_55, CHARLIE1A_te_5...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-12.png)<!-- -->

    ## Rows: 275 Columns: 207
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (206): CHARLIE3_te_18, CHARLIE3_te_19, CHARLIE3_te_25, CHARLIE3_te_26, C...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-13.png)<!-- -->

    ## Rows: 275 Columns: 418
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (417): CHARLIE5_te_4, CHARLIE5_te_16, CHARLIE5_te_17, CHARLIE5_te_67, CH...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-14.png)<!-- -->

    ## Rows: 275 Columns: 30
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (29): CHESHIRE_A_te_8, CHESHIRE_A_te_9, CHESHIRE_A_te_44, CHESHIRE_A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-15.png)<!-- -->

    ## Rows: 275 Columns: 290
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (289): CHESHIRE_te_33, CHESHIRE_te_44, CHESHIRE_te_52, CHESHIRE_te_53, C...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-16.png)<!-- -->

    ## Rows: 275 Columns: 1564
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1563): ERV24_Prim_te_179, ERV24_Prim_te_180, ERV24_Prim_te_185, ERV24_P...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-17.png)<!-- -->

    ## Rows: 275 Columns: 674
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (673): ERVL_te_115, ERVL_te_118, ERVL_te_126, ERVL_te_142, ERVL_te_143, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-18.png)<!-- -->

    ## Rows: 275 Columns: 438
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (437): ERVL-B4_te_12, ERVL-B4_te_13, ERVL-B4_te_110, ERVL-B4_te_111, ERV...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-19.png)<!-- -->

    ## Rows: 275 Columns: 51
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (50): GOLEM_A_te_9, GOLEM_A_te_10, GOLEM_A_te_19, GOLEM_A_te_20, GOLEM_A...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-20.png)<!-- -->

    ## Rows: 275 Columns: 88
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (87): GOLEM_B_te_9, GOLEM_B_te_13, GOLEM_B_te_19, GOLEM_B_te_20, GOLEM_B...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-21.png)<!-- -->

    ## Rows: 275 Columns: 34
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (33): GOLEM_C_te_22, GOLEM_C_te_27, GOLEM_C_te_28, GOLEM_C_te_35, GOLEM_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-22.png)<!-- -->

    ## Rows: 275 Columns: 185
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (184): GOLEM_te_27, GOLEM_te_28, GOLEM_te_34, GOLEM_te_35, GOLEM_te_37, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-23.png)<!-- -->

    ## Rows: 275 Columns: 13
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (12): GSAT_te_55, GSAT_te_56, GSAT_te_59, GSAT_te_60, GSAT_te_105, GSAT_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-24.png)<!-- -->

    ## Rows: 275 Columns: 14
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (13): GSATII_te_19, GSATII_te_27, GSATII_te_42, GSATII_te_56, GSATII_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-25.png)<!-- -->

    ## Rows: 275 Columns: 24
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (23): GSATX_te_11, GSATX_te_22, GSATX_te_49, GSATX_te_55, GSATX_te_60, G...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-26.png)<!-- -->

    ## Rows: 275 Columns: 363
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (362): HAL1M8_te_107, HAL1M8_te_118, HAL1M8_te_126, HAL1M8_te_127, HAL1M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-27.png)<!-- -->

    ## Rows: 275 Columns: 551
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (550): HARLEQUIN_te_19, HARLEQUIN_te_25, HARLEQUIN_te_35, HARLEQUIN_te_3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-28.png)<!-- -->

    ## Rows: 275 Columns: 49
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (48): HARLEQUINLTR_te_61, HARLEQUINLTR_te_96, HARLEQUINLTR_te_108, HARLE...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-29.png)<!-- -->

    ## Rows: 275 Columns: 526
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (525): HERV-K14CI_te_13, HERV-K14CI_te_14, HERV-K14CI_te_15, HERV-K14CI_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-30.png)<!-- -->

    ## Rows: 275 Columns: 556
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (555): HERV-K14I_te_16, HERV-K14I_te_22, HERV-K14I_te_34, HERV-K14I_te_3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-31.png)<!-- -->

    ## Rows: 275 Columns: 38
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (37): HERV1_LTR_te_44, HERV1_LTR_te_47, HERV1_LTR_te_81, HERV1_LTR_te_84...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-32.png)<!-- -->

    ## Rows: 275 Columns: 1182
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1181): HERV15I_te_30, HERV15I_te_53, HERV15I_te_66, HERV15I_te_71, HERV...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-33.png)<!-- -->

    ## Rows: 275 Columns: 515
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (514): HERV17_te_7, HERV17_te_8, HERV17_te_13, HERV17_te_14, HERV17_te_3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-34.png)<!-- -->

    ## Rows: 275 Columns: 1098
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1097): HERV18_te_7, HERV18_te_17, HERV18_te_18, HERV18_te_30, HERV18_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-35.png)<!-- -->

    ## Rows: 275 Columns: 618
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (617): HERV19I_te_9, HERV19I_te_16, HERV19I_te_17, HERV19I_te_27, HERV19...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-36.png)<!-- -->

    ## Rows: 275 Columns: 1681
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1680): HERV3_te_4, HERV3_te_29, HERV3_te_30, HERV3_te_36, HERV3_te_39, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-37.png)<!-- -->

    ## Rows: 275 Columns: 1159
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1158): HERV30I_te_7, HERV30I_te_22, HERV30I_te_44, HERV30I_te_45, HERV3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-38.png)<!-- -->

    ## Rows: 275 Columns: 811
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (810): HERV35I_te_115, HERV35I_te_118, HERV35I_te_119, HERV35I_te_126, H...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-39.png)<!-- -->

    ## Rows: 275 Columns: 315
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (314): HERV38I_te_88, HERV38I_te_92, HERV38I_te_93, HERV38I_te_98, HERV3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-40.png)<!-- -->

    ## Rows: 275 Columns: 268
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (267): HERV39_te_487, HERV39_te_512, HERV39_te_514, HERV39_te_515, HERV3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-41.png)<!-- -->

    ## Rows: 275 Columns: 716
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (715): HERV4_I_te_9, HERV4_I_te_27, HERV4_I_te_39, HERV4_I_te_54, HERV4_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-42.png)<!-- -->

    ## Rows: 275 Columns: 1012
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1011): HERV46I_te_17, HERV46I_te_18, HERV46I_te_20, HERV46I_te_21, HERV...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-43.png)<!-- -->

    ## Rows: 275 Columns: 705
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (704): HERV49I_te_12, HERV49I_te_14, HERV49I_te_20, HERV49I_te_21, HERV4...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-44.png)<!-- -->

    ## Rows: 275 Columns: 993
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (992): HERV57I_te_27, HERV57I_te_30, HERV57I_te_33, HERV57I_te_34, HERV5...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-45.png)<!-- -->

    ## Rows: 275 Columns: 503
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (502): HERV9_te_8, HERV9_te_22, HERV9_te_26, HERV9_te_32, HERV9_te_39, H...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-46.png)<!-- -->

    ## Rows: 275 Columns: 626
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (625): HERVE_te_16, HERVE_te_18, HERVE_te_19, HERVE_te_26, HERVE_te_31, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-47.png)<!-- -->

    ## Rows: 275 Columns: 1001
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1000): HERVFH19I_te_16, HERVFH19I_te_17, HERVFH19I_te_19, HERVFH19I_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-48.png)<!-- -->

    ## Rows: 275 Columns: 1034
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1033): HERVFH21I_te_10, HERVFH21I_te_15, HERVFH21I_te_18, HERVFH21I_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-49.png)<!-- -->

    ## Rows: 275 Columns: 535
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (534): HERVH_te_10, HERVH_te_17, HERVH_te_61, HERVH_te_72, HERVH_te_90, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-50.png)<!-- -->

    ## Rows: 275 Columns: 607
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (606): HERVH48I_te_9, HERVH48I_te_10, HERVH48I_te_16, HERVH48I_te_17, HE...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-51.png)<!-- -->

    ## Rows: 275 Columns: 780
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (779): HERVIP10F_te_13, HERVIP10F_te_14, HERVIP10F_te_17, HERVIP10F_te_1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-52.png)<!-- -->

    ## Rows: 275 Columns: 369
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (368): HERVIP10FH_te_12, HERVIP10FH_te_16, HERVIP10FH_te_34, HERVIP10FH_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-53.png)<!-- -->

    ## Rows: 275 Columns: 572
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (571): HERVK_te_18, HERVK_te_19, HERVK_te_73, HERVK_te_90, HERVK_te_106,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-54.png)<!-- -->

    ## Rows: 275 Columns: 320
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (319): HERVK11DI_te_7, HERVK11DI_te_13, HERVK11DI_te_14, HERVK11DI_te_33...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-55.png)<!-- -->

    ## Rows: 275 Columns: 773
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (772): HERVK11I_te_6, HERVK11I_te_7, HERVK11I_te_13, HERVK11I_te_41, HER...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-56.png)<!-- -->

    ## Rows: 275 Columns: 686
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (685): HERVK22I_te_11, HERVK22I_te_12, HERVK22I_te_22, HERVK22I_te_23, H...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-57.png)<!-- -->

    ## Rows: 275 Columns: 732
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (731): HERVK3I_te_7, HERVK3I_te_8, HERVK3I_te_12, HERVK3I_te_13, HERVK3I...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-58.png)<!-- -->

    ## Rows: 275 Columns: 399
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (398): HERVK9I_te_23, HERVK9I_te_26, HERVK9I_te_33, HERVK9I_te_34, HERVK...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-59.png)<!-- -->

    ## Rows: 275 Columns: 750
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (749): HERVKC4_te_23, HERVKC4_te_24, HERVKC4_te_33, HERVKC4_te_39, HERVK...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-60.png)<!-- -->

    ## Rows: 275 Columns: 242
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (241): HERVL_te_57, HERVL_te_58, HERVL_te_105, HERVL_te_159, HERVL_te_16...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-61.png)<!-- -->

    ## Rows: 275 Columns: 578
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (577): HERVL66I_te_8, HERVL66I_te_9, HERVL66I_te_14, HERVL66I_te_20, HER...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-62.png)<!-- -->

    ## Rows: 275 Columns: 684
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (683): HERVP71A_I_te_13, HERVP71A_I_te_14, HERVP71A_I_te_17, HERVP71A_I_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-63.png)<!-- -->

    ## Rows: 275 Columns: 1389
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1388): HERVS71_te_13, HERVS71_te_16, HERVS71_te_21, HERVS71_te_22, HERV...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-64.png)<!-- -->

    ## Rows: 275 Columns: 63
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (62): HSATI_te_4, HSATI_te_17, HSATI_te_25, HSATI_te_33, HSATI_te_36, HS...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-65.png)<!-- -->

    ## Rows: 275 Columns: 24
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (23): HSATII_te_9, HSATII_te_12, HSATII_te_18, HSATII_te_29, HSATII_te_3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-66.png)<!-- -->

    ## Rows: 275 Columns: 150
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (149): HSMAR1_te_24, HSMAR1_te_25, HSMAR1_te_41, HSMAR1_te_49, HSMAR1_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-67.png)<!-- -->

    ## Rows: 275 Columns: 63
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (62): HSMAR2_te_32, HSMAR2_te_33, HSMAR2_te_66, HSMAR2_te_88, HSMAR2_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-68.png)<!-- -->

    ## Rows: 275 Columns: 380
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (379): HUERS-P1_te_14, HUERS-P1_te_17, HUERS-P1_te_18, HUERS-P1_te_43, H...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-69.png)<!-- -->

    ## Rows: 275 Columns: 305
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (304): HUERS-P2_te_55, HUERS-P2_te_59, HUERS-P2_te_70, HUERS-P2_te_71, H...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-70.png)<!-- -->

    ## Rows: 275 Columns: 1107
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1106): HUERS-P3_te_14, HUERS-P3_te_17, HUERS-P3_te_41, HUERS-P3_te_42, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-71.png)<!-- -->

    ## Rows: 275 Columns: 595
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (594): HUERS-P3B_te_50, HUERS-P3B_te_81, HUERS-P3B_te_86, HUERS-P3B_te_8...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-72.png)<!-- -->

    ## Rows: 275 Columns: 75
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (74): IN25_te_5, IN25_te_16, IN25_te_43, IN25_te_44, IN25_te_66, IN25_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-73.png)<!-- -->

    ## Rows: 275 Columns: 266
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (265): L1_te_6, L1_te_21, L1_te_22, L1_te_54, L1_te_55, L1_te_60, L1_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-74.png)<!-- -->

    ## Rows: 275 Columns: 75
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (74): L1HS_te_145, L1HS_te_151, L1HS_te_155, L1HS_te_199, L1HS_te_292, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-75.png)<!-- -->

    ## Rows: 275 Columns: 256
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (255): L1M1_5_te_21, L1M1_5_te_22, L1M1_5_te_35, L1M1_5_te_36, L1M1_5_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-76.png)<!-- -->

    ## Rows: 275 Columns: 106
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (105): L1M1B_5_te_9, L1M1B_5_te_10, L1M1B_5_te_23, L1M1B_5_te_24, L1M1B_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-77.png)<!-- -->

    ## Rows: 275 Columns: 388
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (387): L1M2_5_te_42, L1M2_5_te_43, L1M2_5_te_72, L1M2_5_te_73, L1M2_5_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-78.png)<!-- -->

    ## Rows: 275 Columns: 289
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (288): L1M2A_5_te_25, L1M2A_5_te_29, L1M2A_5_te_34, L1M2A_5_te_39, L1M2A...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-79.png)<!-- -->

    ## Rows: 275 Columns: 241
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (240): L1M2A1_5_te_20, L1M2A1_5_te_24, L1M2A1_5_te_34, L1M2A1_5_te_38, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-80.png)<!-- -->

    ## Rows: 275 Columns: 305
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (304): L1M2B_5_te_12, L1M2B_5_te_37, L1M2B_5_te_45, L1M2B_5_te_64, L1M2B...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-81.png)<!-- -->

    ## Rows: 275 Columns: 303
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (302): L1M2C_5_te_35, L1M2C_5_te_42, L1M2C_5_te_43, L1M2C_5_te_44, L1M2C...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-82.png)<!-- -->

    ## Rows: 275 Columns: 289
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (288): L1M3A_5_te_288, L1M3A_5_te_290, L1M3A_5_te_319, L1M3A_5_te_320, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-83.png)<!-- -->

    ## Rows: 275 Columns: 167
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (166): L1M3C_5_te_26, L1M3C_5_te_98, L1M3C_5_te_99, L1M3C_5_te_100, L1M3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-84.png)<!-- -->

    ## Rows: 275 Columns: 243
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (242): L1M3DE_5_te_39, L1M3DE_5_te_41, L1M3DE_5_te_48, L1M3DE_5_te_49, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-85.png)<!-- -->

    ## Rows: 275 Columns: 273
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (272): L1M4B_te_41, L1M4B_te_49, L1M4B_te_56, L1M4B_te_57, L1M4B_te_58, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-86.png)<!-- -->

    ## Rows: 275 Columns: 16
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (15): L1MA1_te_59, L1MA1_te_60, L1MA1_te_79, L1MA1_te_127, L1MA1_te_218,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-87.png)<!-- -->

    ## Rows: 275 Columns: 68
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (67): L1MA10_te_22, L1MA10_te_27, L1MA10_te_39, L1MA10_te_42, L1MA10_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-88.png)<!-- -->

    ## Rows: 275 Columns: 58
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (57): L1MA2_te_8, L1MA2_te_39, L1MA2_te_59, L1MA2_te_60, L1MA2_te_64, L1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-89.png)<!-- -->

    ## Rows: 275 Columns: 30
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (29): L1MA3_te_37, L1MA3_te_59, L1MA3_te_60, L1MA3_te_127, L1MA3_te_151,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-90.png)<!-- -->

    ## Rows: 275 Columns: 46
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (45): L1MA4_te_18, L1MA4_te_39, L1MA4_te_59, L1MA4_te_60, L1MA4_te_69, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-91.png)<!-- -->

    ## Rows: 275 Columns: 28
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (27): L1MA4A_te_41, L1MA4A_te_42, L1MA4A_te_162, L1MA4A_te_246, L1MA4A_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-92.png)<!-- -->

    ## Rows: 275 Columns: 39
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (38): L1MA5_te_39, L1MA5_te_59, L1MA5_te_60, L1MA5_te_64, L1MA5_te_125, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-93.png)<!-- -->

    ## Rows: 275 Columns: 27
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (26): L1MA5A_te_59, L1MA5A_te_64, L1MA5A_te_141, L1MA5A_te_250, L1MA5A_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-94.png)<!-- -->

    ## Rows: 275 Columns: 33
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (32): L1MA6_te_42, L1MA6_te_45, L1MA6_te_57, L1MA6_te_64, L1MA6_te_163, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-95.png)<!-- -->

    ## Rows: 275 Columns: 51
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (50): L1MA7_te_9, L1MA7_te_20, L1MA7_te_21, L1MA7_te_22, L1MA7_te_24, L1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-96.png)<!-- -->

    ## Rows: 275 Columns: 107
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (106): L1MA8_te_59, L1MA8_te_60, L1MA8_te_64, L1MA8_te_87, L1MA8_te_162,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-97.png)<!-- -->

    ## Rows: 275 Columns: 217
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (216): L1MA9_5_te_8, L1MA9_5_te_12, L1MA9_5_te_17, L1MA9_5_te_18, L1MA9_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-98.png)<!-- -->

    ## Rows: 275 Columns: 36
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (35): L1MA9_te_29, L1MA9_te_59, L1MA9_te_60, L1MA9_te_64, L1MA9_te_78, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-99.png)<!-- -->

    ## Rows: 275 Columns: 35
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (34): L1MB1_te_6, L1MB1_te_28, L1MB1_te_39, L1MB1_te_41, L1MB1_te_64, L1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-100.png)<!-- -->

    ## Rows: 275 Columns: 35
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (34): L1MB2_te_37, L1MB2_te_45, L1MB2_te_54, L1MB2_te_59, L1MB2_te_60, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-101.png)<!-- -->

    ## Rows: 275 Columns: 112
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (111): L1MB3_5_te_12, L1MB3_5_te_15, L1MB3_5_te_39, L1MB3_5_te_169, L1MB...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-102.png)<!-- -->

    ## Rows: 275 Columns: 45
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (44): L1MB3_te_19, L1MB3_te_22, L1MB3_te_25, L1MB3_te_28, L1MB3_te_34, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-103.png)<!-- -->

    ## Rows: 275 Columns: 135
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (134): L1MB4_5_te_12, L1MB4_5_te_27, L1MB4_5_te_35, L1MB4_5_te_45, L1MB4...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-104.png)<!-- -->

    ## Rows: 275 Columns: 27
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (26): L1MB4_te_54, L1MB4_te_64, L1MB4_te_180, L1MB4_te_196, L1MB4_te_206...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-105.png)<!-- -->

    ## Rows: 275 Columns: 40
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (39): L1MB5_te_64, L1MB5_te_65, L1MB5_te_78, L1MB5_te_81, L1MB5_te_91, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-106.png)<!-- -->

    ## Rows: 275 Columns: 166
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (165): L1MB6_5_te_12, L1MB6_5_te_17, L1MB6_5_te_21, L1MB6_5_te_22, L1MB6...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-107.png)<!-- -->

    ## Rows: 275 Columns: 33
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (32): L1MB7_te_78, L1MB7_te_81, L1MB7_te_87, L1MB7_te_141, L1MB7_te_146,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-108.png)<!-- -->

    ## Rows: 275 Columns: 35
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (34): L1MB8_te_77, L1MB8_te_78, L1MB8_te_81, L1MB8_te_82, L1MB8_te_84, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-109.png)<!-- -->

    ## Rows: 275 Columns: 32
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (31): L1MC1_te_59, L1MC1_te_60, L1MC1_te_64, L1MC1_te_65, L1MC1_te_249, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-110.png)<!-- -->

    ## Rows: 275 Columns: 31
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (30): L1MC2_te_50, L1MC2_te_65, L1MC2_te_141, L1MC2_te_192, L1MC2_te_216...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-111.png)<!-- -->

    ## Rows: 275 Columns: 175
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (174): L1MC3_te_45, L1MC3_te_87, L1MC3_te_93, L1MC3_te_129, L1MC3_te_141...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-112.png)<!-- -->

    ## Rows: 275 Columns: 139
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (138): L1MC4_te_84, L1MC4_te_163, L1MC4_te_165, L1MC4_te_172, L1MC4_te_1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-113.png)<!-- -->

    ## Rows: 275 Columns: 108
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (107): L1MCA_5_te_55, L1MCA_5_te_73, L1MCA_5_te_79, L1MCA_5_te_80, L1MCA...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-114.png)<!-- -->

    ## Rows: 275 Columns: 93
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (92): L1MCB_5_te_58, L1MCB_5_te_61, L1MCB_5_te_62, L1MCB_5_te_73, L1MCB_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-115.png)<!-- -->

    ## Rows: 275 Columns: 62
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (61): L1MD1_te_45, L1MD1_te_59, L1MD1_te_65, L1MD1_te_114, L1MD1_te_155,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-116.png)<!-- -->

    ## Rows: 275 Columns: 49
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (48): L1MD2_te_87, L1MD2_te_100, L1MD2_te_200, L1MD2_te_251, L1MD2_te_25...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-117.png)<!-- -->

    ## Rows: 275 Columns: 84
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (83): L1MD3_te_21, L1MD3_te_25, L1MD3_te_30, L1MD3_te_34, L1MD3_te_36, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-118.png)<!-- -->

    ## Rows: 275 Columns: 202
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (201): L1MDA_5_te_18, L1MDA_5_te_24, L1MDA_5_te_30, L1MDA_5_te_47, L1MDA...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-119.png)<!-- -->

    ## Rows: 275 Columns: 84
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (83): L1MDB_5_te_47, L1MDB_5_te_54, L1MDB_5_te_56, L1MDB_5_te_58, L1MDB_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-120.png)<!-- -->

    ## Rows: 275 Columns: 259
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (258): L1ME_ORF2_te_18, L1ME_ORF2_te_34, L1ME_ORF2_te_52, L1ME_ORF2_te_5...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-121.png)<!-- -->

    ## Rows: 275 Columns: 39
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (38): L1ME1_te_76, L1ME1_te_94, L1ME1_te_109, L1ME1_te_194, L1ME1_te_232...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-122.png)<!-- -->

    ## Rows: 275 Columns: 54
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (53): L1ME2_te_167, L1ME2_te_169, L1ME2_te_172, L1ME2_te_174, L1ME2_te_1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-123.png)<!-- -->

    ## Rows: 275 Columns: 72
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (71): L1ME3_te_203, L1ME3_te_206, L1ME3_te_207, L1ME3_te_209, L1ME3_te_2...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-124.png)<!-- -->

    ## Rows: 275 Columns: 83
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (82): L1ME3A_te_12, L1ME3A_te_20, L1ME3A_te_21, L1ME3A_te_25, L1ME3A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-125.png)<!-- -->

    ## Rows: 275 Columns: 18
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (17): L1ME5_te_255, L1ME5_te_259, L1ME5_te_262, L1ME5_te_264, L1ME5_te_2...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-126.png)<!-- -->

    ## Rows: 275 Columns: 157
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (156): L1MEC_5_te_42, L1MEC_5_te_72, L1MEC_5_te_73, L1MEC_5_te_87, L1MEC...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-127.png)<!-- -->

    ## Rows: 275 Columns: 214
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (213): L1MEf_5end_te_67, L1MEf_5end_te_72, L1MEf_5end_te_98, L1MEf_5end_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-128.png)<!-- -->

    ## Rows: 275 Columns: 178
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (177): L1MEg_5end_te_73, L1MEg_5end_te_74, L1MEg_5end_te_86, L1MEg_5end_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-129.png)<!-- -->

    ## Rows: 275 Columns: 433
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (432): L1P_MA2_te_6, L1P_MA2_te_7, L1P_MA2_te_10, L1P_MA2_te_20, L1P_MA2...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-130.png)<!-- -->

    ## Rows: 275 Columns: 300
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (299): L1P4a_5end_te_15, L1P4a_5end_te_16, L1P4a_5end_te_18, L1P4a_5end_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-131.png)<!-- -->

    ## Rows: 275 Columns: 300
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (299): L1P4c_5end_te_14, L1P4c_5end_te_15, L1P4c_5end_te_21, L1P4c_5end_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-132.png)<!-- -->

    ## Rows: 275 Columns: 343
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (342): L1P4d_5end_te_9, L1P4d_5end_te_12, L1P4d_5end_te_13, L1P4d_5end_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-133.png)<!-- -->

    ## Rows: 275 Columns: 147
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (146): L1P4e_5end_te_14, L1P4e_5end_te_15, L1P4e_5end_te_47, L1P4e_5end_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-134.png)<!-- -->

    ## Rows: 275 Columns: 36
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (35): L1PA10_te_68, L1PA10_te_112, L1PA10_te_116, L1PA10_te_117, L1PA10_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-135.png)<!-- -->

    ## Rows: 275 Columns: 33
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (32): L1PA11_te_12, L1PA11_te_200, L1PA11_te_201, L1PA11_te_250, L1PA11_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-136.png)<!-- -->

    ## Rows: 275 Columns: 345
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (344): L1PA12_5_te_19, L1PA12_5_te_20, L1PA12_5_te_40, L1PA12_5_te_41, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-137.png)<!-- -->

    ## Rows: 275 Columns: 35
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (34): L1PA12_te_57, L1PA12_te_71, L1PA12_te_78, L1PA12_te_91, L1PA12_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-138.png)<!-- -->

    ## Rows: 275 Columns: 263
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (262): L1PA13_5_te_10, L1PA13_5_te_11, L1PA13_5_te_18, L1PA13_5_te_19, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-139.png)<!-- -->

    ## Rows: 275 Columns: 29
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (28): L1PA13_te_96, L1PA13_te_117, L1PA13_te_250, L1PA13_te_251, L1PA13_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-140.png)<!-- -->

    ## Rows: 275 Columns: 210
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (209): L1PA14_5_te_29, L1PA14_5_te_30, L1PA14_5_te_47, L1PA14_5_te_48, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-141.png)<!-- -->

    ## Rows: 275 Columns: 38
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (37): L1PA14_te_60, L1PA14_te_96, L1PA14_te_117, L1PA14_te_146, L1PA14_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-142.png)<!-- -->

    ## Rows: 275 Columns: 53
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (52): L1PA15_te_69, L1PA15_te_96, L1PA15_te_115, L1PA15_te_117, L1PA15_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-143.png)<!-- -->

    ## Rows: 275 Columns: 338
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (337): L1PA16_5_te_8, L1PA16_5_te_18, L1PA16_5_te_19, L1PA16_5_te_33, L1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-144.png)<!-- -->

    ## Rows: 275 Columns: 53
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (52): L1PA16_te_34, L1PA16_te_54, L1PA16_te_69, L1PA16_te_78, L1PA16_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-145.png)<!-- -->

    ## Rows: 275 Columns: 227
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (226): L1PA17_5_te_27, L1PA17_5_te_29, L1PA17_5_te_30, L1PA17_5_te_32, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-146.png)<!-- -->

    ## Rows: 275 Columns: 17
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (16): L1PA2_te_21, L1PA2_te_25, L1PA2_te_28, L1PA2_te_74, L1PA2_te_75, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-147.png)<!-- -->

    ## Rows: 275 Columns: 16
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (15): L1PA3_te_28, L1PA3_te_50, L1PA3_te_61, L1PA3_te_63, L1PA3_te_75, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-148.png)<!-- -->

    ## Rows: 275 Columns: 22
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (21): L1PA4_te_66, L1PA4_te_74, L1PA4_te_215, L1PA4_te_322, L1PA4_te_331...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-149.png)<!-- -->

    ## Rows: 275 Columns: 33
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (32): L1PA6_te_49, L1PA6_te_129, L1PA6_te_200, L1PA6_te_201, L1PA6_te_21...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-150.png)<!-- -->

    ## Rows: 275 Columns: 179
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (178): L1PA7_5_te_19, L1PA7_5_te_20, L1PA7_5_te_35, L1PA7_5_te_36, L1PA7...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-151.png)<!-- -->

    ## Rows: 275 Columns: 38
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (37): L1PA7_te_96, L1PA7_te_98, L1PA7_te_116, L1PA7_te_117, L1PA7_te_200...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-152.png)<!-- -->

    ## Rows: 275 Columns: 62
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (61): L1PA8_te_21, L1PA8_te_34, L1PA8_te_43, L1PA8_te_48, L1PA8_te_57, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-153.png)<!-- -->

    ## Rows: 275 Columns: 37
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (36): L1PB1_te_54, L1PB1_te_133, L1PB1_te_187, L1PB1_te_188, L1PB1_te_19...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-154.png)<!-- -->

    ## Rows: 275 Columns: 33
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (32): L1PB2_te_18, L1PB2_te_54, L1PB2_te_116, L1PB2_te_117, L1PB2_te_223...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-155.png)<!-- -->

    ## Rows: 275 Columns: 435
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (434): L1PB2c_te_12, L1PB2c_te_13, L1PB2c_te_16, L1PB2c_te_17, L1PB2c_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-156.png)<!-- -->

    ## Rows: 275 Columns: 45
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (44): L1PB3_te_49, L1PB3_te_63, L1PB3_te_87, L1PB3_te_201, L1PB3_te_204,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-157.png)<!-- -->

    ## Rows: 275 Columns: 48
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (47): L1PB4_te_18, L1PB4_te_54, L1PB4_te_115, L1PB4_te_116, L1PB4_te_117...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-158.png)<!-- -->

    ## Rows: 275 Columns: 274
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (273): L1PBA_5_te_10, L1PBA_5_te_11, L1PBA_5_te_45, L1PBA_5_te_46, L1PBA...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-159.png)<!-- -->

    ## Rows: 275 Columns: 49
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (48): L1PBA1_5_te_38, L1PBA1_5_te_39, L1PBA1_5_te_54, L1PBA1_5_te_55, L1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-160.png)<!-- -->

    ## Rows: 275 Columns: 80
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (79): L1PBB_5_te_32, L1PBB_5_te_63, L1PBB_5_te_64, L1PBB_5_te_98, L1PBB_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-161.png)<!-- -->

    ## Rows: 275 Columns: 431
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (430): L1PREC1_te_6, L1PREC1_te_21, L1PREC1_te_29, L1PREC1_te_33, L1PREC...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-162.png)<!-- -->

    ## Rows: 275 Columns: 680
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (679): L1PREC2_te_17, L1PREC2_te_18, L1PREC2_te_38, L1PREC2_te_39, L1PRE...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-163.png)<!-- -->

    ## Rows: 275 Columns: 210
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (209): L2_te_219, L2_te_488, L2_te_490, L2_te_593, L2_te_599, L2_te_601,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-164.png)<!-- -->

    ## Rows: 275 Columns: 26
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (25): L2B_te_9, L2B_te_51, L2B_te_56, L2B_te_64, L2B_te_66, L2B_te_67, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-165.png)<!-- -->

    ## Rows: 275 Columns: 35
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (34): LOR1_te_7, LOR1_te_8, LOR1_te_33, LOR1_te_38, LOR1_te_39, LOR1_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-166.png)<!-- -->

    ## Rows: 275 Columns: 35
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (34): LOR1a_LTR_te_85, LOR1a_LTR_te_86, LOR1a_LTR_te_87, LOR1a_LTR_te_97...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-167.png)<!-- -->

    ## Rows: 275 Columns: 45
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (44): LOR1b_LTR_te_102, LOR1b_LTR_te_114, LOR1b_LTR_te_131, LOR1b_LTR_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-168.png)<!-- -->

    ## Rows: 275 Columns: 792
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (791): LOR1I_te_38, LOR1I_te_44, LOR1I_te_47, LOR1I_te_48, LOR1I_te_57, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-169.png)<!-- -->

    ## Rows: 275 Columns: 93
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (92): LSAU_te_21, LSAU_te_22, LSAU_te_32, LSAU_te_35, LSAU_te_74, LSAU_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-170.png)<!-- -->

    ## Rows: 275 Columns: 59
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (58): LTR06_te_25, LTR06_te_26, LTR06_te_57, LTR06_te_58, LTR06_te_60, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-171.png)<!-- -->

    ## Rows: 275 Columns: 73
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (72): LTR1_te_31, LTR1_te_35, LTR1_te_40, LTR1_te_41, LTR1_te_42, LTR1_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-172.png)<!-- -->

    ## Rows: 275 Columns: 61
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (60): LTR10A_te_12, LTR10A_te_24, LTR10A_te_31, LTR10A_te_39, LTR10A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-173.png)<!-- -->

    ## Rows: 275 Columns: 69
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (68): LTR10B_te_19, LTR10B_te_43, LTR10B_te_58, LTR10B_te_67, LTR10B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-174.png)<!-- -->

    ## Rows: 275 Columns: 67
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (66): LTR10B1_te_37, LTR10B1_te_39, LTR10B1_te_40, LTR10B1_te_43, LTR10B...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-175.png)<!-- -->

    ## Rows: 275 Columns: 65
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (64): LTR10B2_te_17, LTR10B2_te_20, LTR10B2_te_31, LTR10B2_te_49, LTR10B...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-176.png)<!-- -->

    ## Rows: 275 Columns: 22
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (21): LTR10C_te_46, LTR10C_te_130, LTR10C_te_131, LTR10C_te_254, LTR10C_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-177.png)<!-- -->

    ## Rows: 275 Columns: 38
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (37): LTR10D_te_32, LTR10D_te_33, LTR10D_te_44, LTR10D_te_72, LTR10D_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-178.png)<!-- -->

    ## Rows: 275 Columns: 60
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (59): LTR10E_te_88, LTR10E_te_97, LTR10E_te_159, LTR10E_te_162, LTR10E_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-179.png)<!-- -->

    ## Rows: 275 Columns: 40
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (39): LTR10F_te_26, LTR10F_te_37, LTR10F_te_91, LTR10F_te_107, LTR10F_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-180.png)<!-- -->

    ## Rows: 275 Columns: 55
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (54): LTR10G_te_14, LTR10G_te_18, LTR10G_te_20, LTR10G_te_40, LTR10G_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-181.png)<!-- -->

    ## Rows: 275 Columns: 84
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (83): LTR12_te_21, LTR12_te_32, LTR12_te_33, LTR12_te_34, LTR12_te_67, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-182.png)<!-- -->

    ## Rows: 275 Columns: 67
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (66): LTR12B_te_32, LTR12B_te_33, LTR12B_te_67, LTR12B_te_96, LTR12B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-183.png)<!-- -->

    ## Rows: 275 Columns: 94
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (93): LTR12C_te_25, LTR12C_te_42, LTR12C_te_52, LTR12C_te_105, LTR12C_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-184.png)<!-- -->

    ## Rows: 275 Columns: 125
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (124): LTR12D_te_34, LTR12D_te_40, LTR12D_te_42, LTR12D_te_43, LTR12D_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-185.png)<!-- -->

    ## Rows: 275 Columns: 116
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (115): LTR12E_te_13, LTR12E_te_14, LTR12E_te_20, LTR12E_te_25, LTR12E_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-186.png)<!-- -->

    ## Rows: 275 Columns: 54
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (53): LTR12F_te_32, LTR12F_te_92, LTR12F_te_93, LTR12F_te_105, LTR12F_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-187.png)<!-- -->

    ## Rows: 275 Columns: 85
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (84): LTR13_te_8, LTR13_te_23, LTR13_te_27, LTR13_te_44, LTR13_te_48, LT...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-188.png)<!-- -->

    ## Rows: 275 Columns: 67
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (66): LTR13A_te_9, LTR13A_te_19, LTR13A_te_93, LTR13A_te_143, LTR13A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-189.png)<!-- -->

    ## Rows: 275 Columns: 36
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (35): LTR14_te_30, LTR14_te_122, LTR14_te_132, LTR14_te_144, LTR14_te_18...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-190.png)<!-- -->

    ## Rows: 275 Columns: 25
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (24): LTR14A_te_23, LTR14A_te_99, LTR14A_te_130, LTR14A_te_131, LTR14A_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-191.png)<!-- -->

    ## Rows: 275 Columns: 52
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (51): LTR14B_te_25, LTR14B_te_82, LTR14B_te_83, LTR14B_te_106, LTR14B_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-192.png)<!-- -->

    ## Rows: 275 Columns: 41
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (40): LTR14C_te_30, LTR14C_te_65, LTR14C_te_75, LTR14C_te_76, LTR14C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-193.png)<!-- -->

    ## Rows: 275 Columns: 35
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (34): LTR15_te_40, LTR15_te_80, LTR15_te_126, LTR15_te_127, LTR15_te_185...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-194.png)<!-- -->

    ## Rows: 275 Columns: 38
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (37): LTR16C_te_96, LTR16C_te_97, LTR16C_te_100, LTR16C_te_106, LTR16C_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-195.png)<!-- -->

    ## Rows: 275 Columns: 91
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (90): LTR17_te_32, LTR17_te_33, LTR17_te_43, LTR17_te_68, LTR17_te_69, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-196.png)<!-- -->

    ## Rows: 275 Columns: 62
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (61): LTR18A_te_8, LTR18A_te_18, LTR18A_te_38, LTR18A_te_39, LTR18A_te_4...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-197.png)<!-- -->

    ## Rows: 275 Columns: 81
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (80): LTR18B_te_39, LTR18B_te_40, LTR18B_te_44, LTR18B_te_45, LTR18B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-198.png)<!-- -->

    ## Rows: 275 Columns: 31
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (30): LTR18C_te_56, LTR18C_te_57, LTR18C_te_71, LTR18C_te_82, LTR18C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-199.png)<!-- -->

    ## Rows: 275 Columns: 45
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (44): LTR19A_te_18, LTR19A_te_19, LTR19A_te_37, LTR19A_te_38, LTR19A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-200.png)<!-- -->

    ## Rows: 275 Columns: 27
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (26): LTR19B_te_19, LTR19B_te_26, LTR19B_te_27, LTR19B_te_37, LTR19B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-201.png)<!-- -->

    ## Rows: 275 Columns: 81
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (80): LTR19C_te_13, LTR19C_te_36, LTR19C_te_46, LTR19C_te_50, LTR19C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-202.png)<!-- -->

    ## Rows: 275 Columns: 121
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (120): LTR1A1_te_6, LTR1A1_te_7, LTR1A1_te_42, LTR1A1_te_43, LTR1A1_te_5...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-203.png)<!-- -->

    ## Rows: 275 Columns: 99
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (98): LTR1A2_te_6, LTR1A2_te_27, LTR1A2_te_29, LTR1A2_te_30, LTR1A2_te_4...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-204.png)<!-- -->

    ## Rows: 275 Columns: 101
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (100): LTR1B_te_43, LTR1B_te_50, LTR1B_te_51, LTR1B_te_76, LTR1B_te_78, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-205.png)<!-- -->

    ## Rows: 275 Columns: 117
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (116): LTR1B0_te_6, LTR1B0_te_7, LTR1B0_te_21, LTR1B0_te_42, LTR1B0_te_4...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-206.png)<!-- -->

    ## Rows: 275 Columns: 139
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (138): LTR1B1_te_58, LTR1B1_te_77, LTR1B1_te_79, LTR1B1_te_83, LTR1B1_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-207.png)<!-- -->

    ## Rows: 275 Columns: 116
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (115): LTR1C_te_23, LTR1C_te_24, LTR1C_te_31, LTR1C_te_32, LTR1C_te_52, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-208.png)<!-- -->

    ## Rows: 275 Columns: 101
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (100): LTR1C1_te_11, LTR1C1_te_12, LTR1C1_te_23, LTR1C1_te_24, LTR1C1_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-209.png)<!-- -->

    ## Rows: 275 Columns: 99
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (98): LTR1C3_te_11, LTR1C3_te_12, LTR1C3_te_23, LTR1C3_te_24, LTR1C3_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-210.png)<!-- -->

    ## Rows: 275 Columns: 139
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (138): LTR1D_te_42, LTR1D_te_43, LTR1D_te_47, LTR1D_te_51, LTR1D_te_54, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-211.png)<!-- -->

    ## Rows: 275 Columns: 179
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (178): LTR1D1_te_42, LTR1D1_te_43, LTR1D1_te_50, LTR1D1_te_51, LTR1D1_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-212.png)<!-- -->

    ## Rows: 275 Columns: 124
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (123): LTR1E_te_51, LTR1E_te_52, LTR1E_te_54, LTR1E_te_61, LTR1E_te_81, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-213.png)<!-- -->

    ## Rows: 275 Columns: 109
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (108): LTR1F_te_6, LTR1F_te_7, LTR1F_te_42, LTR1F_te_43, LTR1F_te_50, LT...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-214.png)<!-- -->

    ## Rows: 275 Columns: 105
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (104): LTR1F1_te_6, LTR1F1_te_7, LTR1F1_te_42, LTR1F1_te_43, LTR1F1_te_5...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-215.png)<!-- -->

    ## Rows: 275 Columns: 115
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (114): LTR1F2_te_6, LTR1F2_te_7, LTR1F2_te_29, LTR1F2_te_37, LTR1F2_te_4...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-216.png)<!-- -->

    ## Rows: 275 Columns: 50
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (49): LTR2_te_34, LTR2_te_40, LTR2_te_41, LTR2_te_44, LTR2_te_53, LTR2_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-217.png)<!-- -->

    ## Rows: 275 Columns: 36
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (35): LTR21A_te_24, LTR21A_te_61, LTR21A_te_67, LTR21A_te_68, LTR21A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-218.png)<!-- -->

    ## Rows: 275 Columns: 54
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (53): LTR21B_te_11, LTR21B_te_15, LTR21B_te_16, LTR21B_te_18, LTR21B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-219.png)<!-- -->

    ## Rows: 275 Columns: 47
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (46): LTR21C_te_21, LTR21C_te_33, LTR21C_te_34, LTR21C_te_41, LTR21C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-220.png)<!-- -->

    ## Rows: 275 Columns: 37
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (36): LTR22_te_64, LTR22_te_91, LTR22_te_136, LTR22_te_137, LTR22_te_184...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-221.png)<!-- -->

    ## Rows: 275 Columns: 46
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (45): LTR22A_te_11, LTR22A_te_12, LTR22A_te_28, LTR22A_te_46, LTR22A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-222.png)<!-- -->

    ## Rows: 275 Columns: 69
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (68): LTR22B_te_13, LTR22B_te_51, LTR22B_te_56, LTR22B_te_63, LTR22B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-223.png)<!-- -->

    ## Rows: 275 Columns: 71
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (70): LTR22B1_te_7, LTR22B1_te_41, LTR22B1_te_53, LTR22B1_te_57, LTR22B1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-224.png)<!-- -->

    ## Rows: 275 Columns: 44
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (43): LTR22B2_te_56, LTR22B2_te_59, LTR22B2_te_74, LTR22B2_te_78, LTR22B...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-225.png)<!-- -->

    ## Rows: 275 Columns: 74
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (73): LTR22C_te_7, LTR22C_te_9, LTR22C_te_12, LTR22C_te_28, LTR22C_te_34...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-226.png)<!-- -->

    ## Rows: 275 Columns: 68
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (67): LTR22C2_te_7, LTR22C2_te_9, LTR22C2_te_14, LTR22C2_te_18, LTR22C2_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-227.png)<!-- -->

    ## Rows: 275 Columns: 64
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (63): LTR22E_te_5, LTR22E_te_12, LTR22E_te_35, LTR22E_te_59, LTR22E_te_6...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-228.png)<!-- -->

    ## Rows: 275 Columns: 20
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (19): LTR23_te_97, LTR23_te_101, LTR23_te_103, LTR23_te_104, LTR23_te_10...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-229.png)<!-- -->

    ## Rows: 275 Columns: 47
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (46): LTR24_te_50, LTR24_te_77, LTR24_te_78, LTR24_te_89, LTR24_te_90, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-230.png)<!-- -->

    ## Rows: 275 Columns: 36
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (35): LTR24B_te_8, LTR24B_te_78, LTR24B_te_98, LTR24B_te_113, LTR24B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-231.png)<!-- -->

    ## Rows: 275 Columns: 41
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (40): LTR24C_te_25, LTR24C_te_40, LTR24C_te_45, LTR24C_te_61, LTR24C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-232.png)<!-- -->

    ## Rows: 275 Columns: 104
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (103): LTR25_te_13, LTR25_te_17, LTR25_te_33, LTR25_te_49, LTR25_te_51, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-233.png)<!-- -->

    ## Rows: 275 Columns: 1111
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1110): LTR25-int_te_5, LTR25-int_te_8, LTR25-int_te_27, LTR25-int_te_31...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-234.png)<!-- -->

    ## Rows: 275 Columns: 51
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (50): LTR26_te_7, LTR26_te_8, LTR26_te_45, LTR26_te_58, LTR26_te_59, LTR...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-235.png)<!-- -->

    ## Rows: 275 Columns: 60
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (59): LTR26B_te_73, LTR26B_te_80, LTR26B_te_107, LTR26B_te_110, LTR26B_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-236.png)<!-- -->

    ## Rows: 275 Columns: 49
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (48): LTR26C_te_7, LTR26C_te_8, LTR26C_te_15, LTR26C_te_34, LTR26C_te_10...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-237.png)<!-- -->

    ## Rows: 275 Columns: 46
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (45): LTR26D_te_38, LTR26D_te_39, LTR26D_te_89, LTR26D_te_118, LTR26D_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-238.png)<!-- -->

    ## Rows: 275 Columns: 70
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (69): LTR26E_te_12, LTR26E_te_13, LTR26E_te_15, LTR26E_te_20, LTR26E_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-239.png)<!-- -->

    ## Rows: 275 Columns: 67
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (66): LTR27_te_34, LTR27_te_38, LTR27_te_55, LTR27_te_61, LTR27_te_70, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-240.png)<!-- -->

    ## Rows: 275 Columns: 260
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (259): LTR2752_te_7, LTR2752_te_8, LTR2752_te_24, LTR2752_te_33, LTR2752...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-241.png)<!-- -->

    ## Rows: 275 Columns: 73
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (72): LTR27B_te_7, LTR27B_te_43, LTR27B_te_51, LTR27B_te_52, LTR27B_te_1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-242.png)<!-- -->

    ## Rows: 275 Columns: 80
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (79): LTR27C_te_26, LTR27C_te_37, LTR27C_te_43, LTR27C_te_44, LTR27C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-243.png)<!-- -->

    ## Rows: 275 Columns: 98
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (97): LTR27D_te_7, LTR27D_te_8, LTR27D_te_43, LTR27D_te_44, LTR27D_te_51...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-244.png)<!-- -->

    ## Rows: 275 Columns: 120
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (119): LTR27E_te_7, LTR27E_te_8, LTR27E_te_9, LTR27E_te_42, LTR27E_te_43...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-245.png)<!-- -->

    ## Rows: 275 Columns: 109
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (108): LTR28_te_7, LTR28_te_8, LTR28_te_20, LTR28_te_43, LTR28_te_44, LT...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-246.png)<!-- -->

    ## Rows: 275 Columns: 155
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (154): LTR28B_te_42, LTR28B_te_43, LTR28B_te_50, LTR28B_te_51, LTR28B_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-247.png)<!-- -->

    ## Rows: 275 Columns: 185
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (184): LTR28C_te_7, LTR28C_te_8, LTR28C_te_26, LTR28C_te_37, LTR28C_te_6...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-248.png)<!-- -->

    ## Rows: 275 Columns: 47
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (46): LTR29_te_18, LTR29_te_19, LTR29_te_36, LTR29_te_37, LTR29_te_38, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-249.png)<!-- -->

    ## Rows: 275 Columns: 25
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (24): LTR2B_te_91, LTR2B_te_95, LTR2B_te_112, LTR2B_te_124, LTR2B_te_156...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-250.png)<!-- -->

    ## Rows: 275 Columns: 54
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (53): LTR2C_te_16, LTR2C_te_45, LTR2C_te_53, LTR2C_te_67, LTR2C_te_72, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-251.png)<!-- -->

    ## Rows: 275 Columns: 37
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (36): LTR3_te_11, LTR3_te_26, LTR3_te_27, LTR3_te_63, LTR3_te_93, LTR3_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-252.png)<!-- -->

    ## Rows: 275 Columns: 64
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (63): LTR30_te_11, LTR30_te_12, LTR30_te_33, LTR30_te_89, LTR30_te_90, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-253.png)<!-- -->

    ## Rows: 275 Columns: 84
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (83): LTR31_te_10, LTR31_te_11, LTR31_te_32, LTR31_te_33, LTR31_te_36, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-254.png)<!-- -->

    ## Rows: 275 Columns: 50
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (49): LTR32_te_22, LTR32_te_24, LTR32_te_26, LTR32_te_35, LTR32_te_58, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-255.png)<!-- -->

    ## Rows: 275 Columns: 45
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (44): LTR34_te_61, LTR34_te_106, LTR34_te_110, LTR34_te_117, LTR34_te_13...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-256.png)<!-- -->

    ## Rows: 275 Columns: 72
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (71): LTR35_te_14, LTR35_te_32, LTR35_te_38, LTR35_te_40, LTR35_te_42, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-257.png)<!-- -->

    ## Rows: 275 Columns: 38
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (37): LTR35A_te_14, LTR35A_te_17, LTR35A_te_18, LTR35A_te_102, LTR35A_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-258.png)<!-- -->

    ## Rows: 275 Columns: 92
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (91): LTR35B_te_35, LTR35B_te_51, LTR35B_te_53, LTR35B_te_66, LTR35B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-259.png)<!-- -->

    ## Rows: 275 Columns: 65
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (64): LTR36_te_7, LTR36_te_13, LTR36_te_16, LTR36_te_17, LTR36_te_36, LT...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-260.png)<!-- -->

    ## Rows: 275 Columns: 49
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (48): LTR37A_te_15, LTR37A_te_23, LTR37A_te_30, LTR37A_te_45, LTR37A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-261.png)<!-- -->

    ## Rows: 275 Columns: 61
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (60): LTR37B_te_20, LTR37B_te_27, LTR37B_te_35, LTR37B_te_42, LTR37B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-262.png)<!-- -->

    ## Rows: 275 Columns: 63
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (62): LTR38_te_7, LTR38_te_15, LTR38_te_51, LTR38_te_56, LTR38_te_120, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-263.png)<!-- -->

    ## Rows: 275 Columns: 72
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (71): LTR38A1_te_7, LTR38A1_te_8, LTR38A1_te_62, LTR38A1_te_66, LTR38A1_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-264.png)<!-- -->

    ## Rows: 275 Columns: 47
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (46): LTR38B_te_7, LTR38B_te_48, LTR38B_te_60, LTR38B_te_103, LTR38B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-265.png)<!-- -->

    ## Rows: 275 Columns: 90
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (89): LTR38C_te_17, LTR38C_te_29, LTR38C_te_30, LTR38C_te_40, LTR38C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-266.png)<!-- -->

    ## Rows: 275 Columns: 73
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (72): LTR39_te_17, LTR39_te_29, LTR39_te_36, LTR39_te_37, LTR39_te_43, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-267.png)<!-- -->

    ## Rows: 275 Columns: 44
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (43): LTR3A_te_12, LTR3A_te_19, LTR3A_te_20, LTR3A_te_26, LTR3A_te_99, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-268.png)<!-- -->

    ## Rows: 275 Columns: 51
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (50): LTR3B_te_6, LTR3B_te_12, LTR3B_te_18, LTR3B_te_27, LTR3B_te_63, LT...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-269.png)<!-- -->

    ## Rows: 275 Columns: 76
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (75): LTR4_te_14, LTR4_te_20, LTR4_te_27, LTR4_te_31, LTR4_te_54, LTR4_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-270.png)<!-- -->

    ## Rows: 275 Columns: 53
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (52): LTR40A_te_30, LTR40A_te_31, LTR40A_te_40, LTR40A_te_43, LTR40A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-271.png)<!-- -->

    ## Rows: 275 Columns: 32
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (31): LTR40B_te_79, LTR40B_te_84, LTR40B_te_90, LTR40B_te_96, LTR40B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-272.png)<!-- -->

    ## Rows: 275 Columns: 47
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (46): LTR41_te_47, LTR41_te_48, LTR41_te_51, LTR41_te_52, LTR41_te_61, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-273.png)<!-- -->

    ## Rows: 275 Columns: 44
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (43): LTR42_te_24, LTR42_te_25, LTR42_te_26, LTR42_te_37, LTR42_te_47, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-274.png)<!-- -->

    ## Rows: 275 Columns: 46
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (45): LTR43_te_28, LTR43_te_36, LTR43_te_39, LTR43_te_40, LTR43_te_44, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-275.png)<!-- -->

    ## Rows: 275 Columns: 65
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (64): LTR43B_te_42, LTR43B_te_55, LTR43B_te_59, LTR43B_te_60, LTR43B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-276.png)<!-- -->

    ## Rows: 275 Columns: 54
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (53): LTR44_te_41, LTR44_te_48, LTR44_te_49, LTR44_te_52, LTR44_te_60, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-277.png)<!-- -->

    ## Rows: 275 Columns: 55
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (54): LTR45_te_7, LTR45_te_8, LTR45_te_9, LTR45_te_47, LTR45_te_54, LTR4...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-278.png)<!-- -->

    ## Rows: 275 Columns: 39
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (38): LTR45B_te_94, LTR45B_te_95, LTR45B_te_124, LTR45B_te_125, LTR45B_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-279.png)<!-- -->

    ## Rows: 275 Columns: 51
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (50): LTR45C_te_29, LTR45C_te_42, LTR45C_te_74, LTR45C_te_105, LTR45C_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-280.png)<!-- -->

    ## Rows: 275 Columns: 46
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (45): LTR46_te_31, LTR46_te_35, LTR46_te_38, LTR46_te_70, LTR46_te_73, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-281.png)<!-- -->

    ## Rows: 275 Columns: 35
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (34): LTR47A_te_38, LTR47A_te_39, LTR47A_te_48, LTR47A_te_64, LTR47A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-282.png)<!-- -->

    ## Rows: 275 Columns: 31
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (30): LTR47A2_te_48, LTR47A2_te_49, LTR47A2_te_54, LTR47A2_te_55, LTR47A...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-283.png)<!-- -->

    ## Rows: 275 Columns: 76
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (75): LTR47B_te_15, LTR47B_te_30, LTR47B_te_38, LTR47B_te_39, LTR47B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-284.png)<!-- -->

    ## Rows: 275 Columns: 45
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (44): LTR47B2_te_38, LTR47B2_te_39, LTR47B2_te_69, LTR47B2_te_70, LTR47B...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-285.png)<!-- -->

    ## Rows: 275 Columns: 61
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (60): LTR47B3_te_38, LTR47B3_te_39, LTR47B3_te_69, LTR47B3_te_70, LTR47B...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-286.png)<!-- -->

    ## Rows: 275 Columns: 48
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (47): LTR47B4_te_38, LTR47B4_te_39, LTR47B4_te_46, LTR47B4_te_62, LTR47B...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-287.png)<!-- -->

    ## Rows: 275 Columns: 68
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (67): LTR48_te_50, LTR48_te_59, LTR48_te_67, LTR48_te_80, LTR48_te_87, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-288.png)<!-- -->

    ## Rows: 275 Columns: 58
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (57): LTR48B_te_14, LTR48B_te_15, LTR48B_te_16, LTR48B_te_52, LTR48B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-289.png)<!-- -->

    ## Rows: 275 Columns: 29
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (28): LTR49_te_11, LTR49_te_21, LTR49_te_35, LTR49_te_69, LTR49_te_73, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-290.png)<!-- -->

    ## Rows: 275 Columns: 43
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (42): LTR5_Hs_te_13, LTR5_Hs_te_79, LTR5_Hs_te_107, LTR5_Hs_te_111, LTR5...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-291.png)<!-- -->

    ## Rows: 275 Columns: 84
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (83): LTR5_te_21, LTR5_te_27, LTR5_te_41, LTR5_te_43, LTR5_te_62, LTR5_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-292.png)<!-- -->

    ## Rows: 275 Columns: 52
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (51): LTR51_te_25, LTR51_te_26, LTR51_te_45, LTR51_te_46, LTR51_te_47, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-293.png)<!-- -->

    ## Rows: 275 Columns: 38
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (37): LTR52_te_64, LTR52_te_72, LTR52_te_77, LTR52_te_112, LTR52_te_130,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-294.png)<!-- -->

    ## Rows: 275 Columns: 58
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (57): LTR53_te_27, LTR53_te_57, LTR53_te_62, LTR53_te_63, LTR53_te_65, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-295.png)<!-- -->

    ## Rows: 275 Columns: 68
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (67): LTR53B_te_11, LTR53B_te_12, LTR53B_te_14, LTR53B_te_18, LTR53B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-296.png)<!-- -->

    ## Rows: 275 Columns: 43
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (42): LTR54_te_13, LTR54_te_30, LTR54_te_39, LTR54_te_97, LTR54_te_106, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-297.png)<!-- -->

    ## Rows: 275 Columns: 30
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (29): LTR54B_te_13, LTR54B_te_78, LTR54B_te_107, LTR54B_te_108, LTR54B_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-298.png)<!-- -->

    ## Rows: 275 Columns: 48
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (47): LTR55_te_16, LTR55_te_20, LTR55_te_33, LTR55_te_64, LTR55_te_74, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-299.png)<!-- -->

    ## Rows: 275 Columns: 21
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (20): LTR56_te_8, LTR56_te_16, LTR56_te_36, LTR56_te_37, LTR56_te_42, LT...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-300.png)<!-- -->

    ## Rows: 275 Columns: 40
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (39): LTR57_te_12, LTR57_te_21, LTR57_te_47, LTR57_te_48, LTR57_te_53, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-301.png)<!-- -->

    ## Rows: 275 Columns: 102
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (101): LTR58_te_8, LTR58_te_10, LTR58_te_14, LTR58_te_26, LTR58_te_36, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-302.png)<!-- -->

    ## Rows: 275 Columns: 54
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (53): LTR59_te_15, LTR59_te_27, LTR59_te_49, LTR59_te_58, LTR59_te_59, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-303.png)<!-- -->

    ## Rows: 275 Columns: 63
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (62): LTR5A_te_8, LTR5A_te_31, LTR5A_te_49, LTR5A_te_52, LTR5A_te_107, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-304.png)<!-- -->

    ## Rows: 275 Columns: 86
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (85): LTR5B_te_50, LTR5B_te_66, LTR5B_te_77, LTR5B_te_85, LTR5B_te_96, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-305.png)<!-- -->

    ## Rows: 275 Columns: 76
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (75): LTR60_te_34, LTR60_te_40, LTR60_te_57, LTR60_te_64, LTR60_te_89, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-306.png)<!-- -->

    ## Rows: 275 Columns: 73
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (72): LTR60B_te_39, LTR60B_te_40, LTR60B_te_50, LTR60B_te_83, LTR60B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-307.png)<!-- -->

    ## Rows: 275 Columns: 56
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (55): LTR61_te_39, LTR61_te_64, LTR61_te_89, LTR61_te_104, LTR61_te_110,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-308.png)<!-- -->

    ## Rows: 275 Columns: 76
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (75): LTR62_te_24, LTR62_te_26, LTR62_te_63, LTR62_te_70, LTR62_te_71, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-309.png)<!-- -->

    ## Rows: 275 Columns: 45
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (44): LTR64_te_14, LTR64_te_34, LTR64_te_54, LTR64_te_85, LTR64_te_88, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-310.png)<!-- -->

    ## Rows: 275 Columns: 106
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (105): LTR65_te_8, LTR65_te_22, LTR65_te_25, LTR65_te_26, LTR65_te_29, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-311.png)<!-- -->

    ## Rows: 275 Columns: 48
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (47): LTR66_te_43, LTR66_te_71, LTR66_te_72, LTR66_te_100, LTR66_te_110,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-312.png)<!-- -->

    ## Rows: 275 Columns: 77
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (76): LTR69_te_7, LTR69_te_25, LTR69_te_26, LTR69_te_32, LTR69_te_65, LT...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-313.png)<!-- -->

    ## Rows: 275 Columns: 72
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (71): LTR6A_te_10, LTR6A_te_24, LTR6A_te_25, LTR6A_te_53, LTR6A_te_56, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-314.png)<!-- -->

    ## Rows: 275 Columns: 32
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (31): LTR6B_te_16, LTR6B_te_52, LTR6B_te_222, LTR6B_te_226, LTR6B_te_233...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-315.png)<!-- -->

    ## Rows: 275 Columns: 88
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (87): LTR70_te_6, LTR70_te_14, LTR70_te_46, LTR70_te_50, LTR70_te_65, LT...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-316.png)<!-- -->

    ## Rows: 275 Columns: 50
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (49): LTR71A_te_51, LTR71A_te_69, LTR71A_te_70, LTR71A_te_75, LTR71A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-317.png)<!-- -->

    ## Rows: 275 Columns: 47
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (46): LTR71B_te_49, LTR71B_te_53, LTR71B_te_54, LTR71B_te_112, LTR71B_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-318.png)<!-- -->

    ## Rows: 275 Columns: 47
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (46): LTR72_te_33, LTR72_te_85, LTR72_te_100, LTR72_te_103, LTR72_te_111...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-319.png)<!-- -->

    ## Rows: 275 Columns: 45
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (44): LTR72B_te_22, LTR72B_te_23, LTR72B_te_31, LTR72B_te_37, LTR72B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-320.png)<!-- -->

    ## Rows: 275 Columns: 49
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (48): LTR73_te_52, LTR73_te_83, LTR73_te_84, LTR73_te_89, LTR73_te_90, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-321.png)<!-- -->

    ## Rows: 275 Columns: 81
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (80): LTR75_1_te_13, LTR75_1_te_14, LTR75_1_te_33, LTR75_1_te_41, LTR75_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-322.png)<!-- -->

    ## Rows: 275 Columns: 57
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (56): LTR76_te_22, LTR76_te_48, LTR76_te_58, LTR76_te_80, LTR76_te_133, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-323.png)<!-- -->

    ## Rows: 275 Columns: 63
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (62): LTR77_te_27, LTR77_te_44, LTR77_te_60, LTR77_te_61, LTR77_te_62, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-324.png)<!-- -->

    ## Rows: 275 Columns: 60
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (59): LTR77B_te_52, LTR77B_te_53, LTR77B_te_57, LTR77B_te_100, LTR77B_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-325.png)<!-- -->

    ## Rows: 275 Columns: 28
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (27): LTR7A_te_32, LTR7A_te_52, LTR7A_te_57, LTR7A_te_101, LTR7A_te_114,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-326.png)<!-- -->

    ## Rows: 275 Columns: 16
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (15): LTR7B_te_52, LTR7B_te_55, LTR7B_te_85, LTR7B_te_86, LTR7B_te_142, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-327.png)<!-- -->

    ## Rows: 275 Columns: 49
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (48): LTR7C_te_17, LTR7C_te_18, LTR7C_te_32, LTR7C_te_57, LTR7C_te_58, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-328.png)<!-- -->

    ## Rows: 275 Columns: 19
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (18): LTR7Y_te_20, LTR7Y_te_25, LTR7Y_te_58, LTR7Y_te_59, LTR7Y_te_104, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-329.png)<!-- -->

    ## Rows: 275 Columns: 46
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (45): LTR8_te_7, LTR8_te_8, LTR8_te_47, LTR8_te_50, LTR8_te_55, LTR8_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-330.png)<!-- -->

    ## Rows: 275 Columns: 60
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (59): LTR8A_te_61, LTR8A_te_72, LTR8A_te_91, LTR8A_te_104, LTR8A_te_105,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-331.png)<!-- -->

    ## Rows: 275 Columns: 69
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (68): LTR8B_te_7, LTR8B_te_8, LTR8B_te_38, LTR8B_te_43, LTR8B_te_88, LTR...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-332.png)<!-- -->

    ## Rows: 275 Columns: 88
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (87): LTR9_te_6, LTR9_te_24, LTR9_te_40, LTR9_te_50, LTR9_te_57, LTR9_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-333.png)<!-- -->

    ## Rows: 275 Columns: 88
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (87): LTR9A1_te_12, LTR9A1_te_45, LTR9A1_te_57, LTR9A1_te_58, LTR9A1_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-334.png)<!-- -->

    ## Rows: 275 Columns: 65
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (64): LTR9B_te_60, LTR9B_te_68, LTR9B_te_69, LTR9B_te_95, LTR9B_te_96, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-335.png)<!-- -->

    ## Rows: 275 Columns: 68
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (67): LTR9C_te_24, LTR9C_te_59, LTR9C_te_60, LTR9C_te_128, LTR9C_te_135,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-336.png)<!-- -->

    ## Rows: 275 Columns: 75
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (74): LTR9D_te_6, LTR9D_te_12, LTR9D_te_35, LTR9D_te_49, LTR9D_te_50, LT...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-337.png)<!-- -->

    ## Rows: 275 Columns: 1158
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1157): MER101_I_te_456, MER101_I_te_457, MER101_I_te_461, MER101_I_te_4...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-338.png)<!-- -->

    ## Rows: 275 Columns: 36
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (35): MER101_te_30, MER101_te_48, MER101_te_49, MER101_te_87, MER101_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-339.png)<!-- -->

    ## Rows: 275 Columns: 78
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (77): MER101B_te_18, MER101B_te_21, MER101B_te_22, MER101B_te_42, MER101...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-340.png)<!-- -->

    ## Rows: 275 Columns: 9
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): familyname_position
    ## dbl (8): MER103_te_23, MER103_te_25, MER103_te_32, MER103_te_34, MER103_te_3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-341.png)<!-- -->

    ## Rows: 275 Columns: 16
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (15): MER105_te_37, MER105_te_39, MER105_te_42, MER105_te_67, MER105_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-342.png)<!-- -->

    ## Rows: 275 Columns: 24
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (23): MER106_te_36, MER106_te_49, MER106_te_57, MER106_te_60, MER106_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-343.png)<!-- -->

    ## Rows: 275 Columns: 5
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): familyname_position
    ## dbl (4): MER107_te_91, MER107_te_92, MER107_te_181, MER107_te_182
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-344.png)<!-- -->

    ## Rows: 275 Columns: 48
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (47): MER119_te_132, MER119_te_139, MER119_te_169, MER119_te_170, MER119...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-345.png)<!-- -->

    ## Rows: 275 Columns: 94
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (93): MER11A_te_6, MER11A_te_7, MER11A_te_23, MER11A_te_24, MER11A_te_27...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-346.png)<!-- -->

    ## Rows: 275 Columns: 112
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (111): MER11B_te_6, MER11B_te_7, MER11B_te_24, MER11B_te_27, MER11B_te_2...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-347.png)<!-- -->

    ## Rows: 275 Columns: 61
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (60): MER11C_te_28, MER11C_te_37, MER11C_te_50, MER11C_te_59, MER11C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-348.png)<!-- -->

    ## Rows: 275 Columns: 51
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (50): MER11D_te_30, MER11D_te_31, MER11D_te_120, MER11D_te_131, MER11D_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-349.png)<!-- -->

    ## Rows: 275 Columns: 7
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): familyname_position
    ## dbl (6): MER121_te_99, MER121_te_108, MER121_te_137, MER121_te_207, MER121_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-350.png)<!-- -->

    ## Rows: 275 Columns: 61
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (60): MER122_te_5, MER122_te_14, MER122_te_19, MER122_te_35, MER122_te_3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-351.png)<!-- -->

    ## Rows: 275 Columns: 80
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (79): MER1A_te_19, MER1A_te_25, MER1A_te_26, MER1A_te_27, MER1A_te_30, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-352.png)<!-- -->

    ## Rows: 275 Columns: 52
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (51): MER1B_te_18, MER1B_te_19, MER1B_te_23, MER1B_te_24, MER1B_te_25, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-353.png)<!-- -->

    ## Rows: 275 Columns: 36
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (35): MER2_te_12, MER2_te_13, MER2_te_19, MER2_te_20, MER2_te_42, MER2_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-354.png)<!-- -->

    ## Rows: 275 Columns: 6
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): familyname_position
    ## dbl (5): MER20_te_15, MER20_te_16, MER20_te_92, MER20_te_179, MER20_te_180
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-355.png)<!-- -->

    ## Rows: 275 Columns: 104
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (103): MER21_te_33, MER21_te_37, MER21_te_43, MER21_te_50, MER21_te_52, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-356.png)<!-- -->

    ## Rows: 275 Columns: 71
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (70): MER21A_te_37, MER21A_te_97, MER21A_te_98, MER21A_te_111, MER21A_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-357.png)<!-- -->

    ## Rows: 275 Columns: 67
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (66): MER21B_te_31, MER21B_te_32, MER21B_te_33, MER21B_te_37, MER21B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-358.png)<!-- -->

    ## Rows: 275 Columns: 70
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (69): MER21C_BT_te_30, MER21C_BT_te_39, MER21C_BT_te_65, MER21C_BT_te_78...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-359.png)<!-- -->

    ## Rows: 275 Columns: 59
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (58): MER21C_te_35, MER21C_te_40, MER21C_te_94, MER21C_te_101, MER21C_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-360.png)<!-- -->

    ## Rows: 275 Columns: 502
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (501): MER21I_te_17, MER21I_te_18, MER21I_te_19, MER21I_te_20, MER21I_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-361.png)<!-- -->

    ## Rows: 275 Columns: 100
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (99): MER22_te_44, MER22_te_72, MER22_te_92, MER22_te_108, MER22_te_112,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-362.png)<!-- -->

    ## Rows: 275 Columns: 41
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (40): MER28_te_34, MER28_te_35, MER28_te_39, MER28_te_50, MER28_te_51, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-363.png)<!-- -->

    ## Rows: 275 Columns: 27
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (26): MER2B_te_19, MER2B_te_45, MER2B_te_46, MER2B_te_47, MER2B_te_60, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-364.png)<!-- -->

    ## Rows: 275 Columns: 11
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (10): MER3_te_4, MER3_te_28, MER3_te_51, MER3_te_52, MER3_te_63, MER3_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-365.png)<!-- -->

    ## Rows: 275 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (22): MER30_te_37, MER30_te_88, MER30_te_89, MER30_te_104, MER30_te_117,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-366.png)<!-- -->

    ## Rows: 275 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (22): MER30B_te_11, MER30B_te_20, MER30B_te_21, MER30B_te_25, MER30B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-367.png)<!-- -->

    ## Rows: 275 Columns: 50
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (49): MER31_te_15, MER31_te_32, MER31_te_33, MER31_te_37, MER31_te_38, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-368.png)<!-- -->

    ## Rows: 275 Columns: 26
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (25): MER31A_te_14, MER31A_te_21, MER31A_te_40, MER31A_te_133, MER31A_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-369.png)<!-- -->

    ## Rows: 275 Columns: 43
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (42): MER31B_te_11, MER31B_te_42, MER31B_te_55, MER31B_te_70, MER31B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-370.png)<!-- -->

    ## Rows: 275 Columns: 16
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (15): MER33_te_4, MER33_te_16, MER33_te_17, MER33_te_78, MER33_te_97, ME...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-371.png)<!-- -->

    ## Rows: 275 Columns: 63
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (62): MER34_te_52, MER34_te_79, MER34_te_89, MER34_te_99, MER34_te_103, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-372.png)<!-- -->

    ## Rows: 275 Columns: 56
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (55): MER34A_te_72, MER34A_te_88, MER34A_te_95, MER34A_te_109, MER34A_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-373.png)<!-- -->

    ## Rows: 275 Columns: 47
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (46): MER34A1_te_142, MER34A1_te_143, MER34A1_te_154, MER34A1_te_180, ME...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-374.png)<!-- -->

    ## Rows: 275 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (22): MER34B_te_42, MER34B_te_52, MER34B_te_66, MER34B_te_165, MER34B_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-375.png)<!-- -->

    ## Rows: 275 Columns: 21
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (20): MER34C_te_80, MER34C_te_132, MER34C_te_133, MER34C_te_186, MER34C_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-376.png)<!-- -->

    ## Rows: 275 Columns: 32
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (31): MER34C2_te_70, MER34C2_te_72, MER34C2_te_107, MER34C2_te_146, MER3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-377.png)<!-- -->

    ## Rows: 275 Columns: 39
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (38): MER34D_te_7, MER34D_te_10, MER34D_te_17, MER34D_te_20, MER34D_te_2...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-378.png)<!-- -->

    ## Rows: 275 Columns: 36
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (35): MER39_te_18, MER39_te_19, MER39_te_64, MER39_te_139, MER39_te_140,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-379.png)<!-- -->

    ## Rows: 275 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (22): MER39B_te_43, MER39B_te_114, MER39B_te_115, MER39B_te_154, MER39B_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-380.png)<!-- -->

    ## Rows: 275 Columns: 48
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (47): MER41A_te_9, MER41A_te_11, MER41A_te_36, MER41A_te_56, MER41A_te_6...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-381.png)<!-- -->

    ## Rows: 275 Columns: 49
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (48): MER41B_te_10, MER41B_te_11, MER41B_te_56, MER41B_te_87, MER41B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-382.png)<!-- -->

    ## Rows: 275 Columns: 54
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (53): MER41C_te_11, MER41C_te_14, MER41C_te_15, MER41C_te_24, MER41C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-383.png)<!-- -->

    ## Rows: 275 Columns: 66
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (65): MER41D_te_11, MER41D_te_32, MER41D_te_63, MER41D_te_94, MER41D_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-384.png)<!-- -->

    ## Rows: 275 Columns: 61
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (60): MER41E_te_12, MER41E_te_13, MER41E_te_26, MER41E_te_27, MER41E_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-385.png)<!-- -->

    ## Rows: 275 Columns: 56
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (55): MER41F_te_27, MER41F_te_33, MER41F_te_38, MER41F_te_40, MER41F_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-386.png)<!-- -->

    ## Rows: 275 Columns: 90
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (89): MER41G_te_23, MER41G_te_26, MER41G_te_27, MER41G_te_64, MER41G_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-387.png)<!-- -->

    ## Rows: 275 Columns: 404
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (403): MER41I_te_28, MER41I_te_42, MER41I_te_43, MER41I_te_44, MER41I_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-388.png)<!-- -->

    ## Rows: 275 Columns: 22
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (21): MER44A_te_19, MER44A_te_20, MER44A_te_21, MER44A_te_22, MER44A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-389.png)<!-- -->

    ## Rows: 275 Columns: 49
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (48): MER44B_te_120, MER44B_te_121, MER44B_te_122, MER44B_te_123, MER44B...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-390.png)<!-- -->

    ## Rows: 275 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (22): MER45_te_30, MER45_te_35, MER45_te_36, MER45_te_45, MER45_te_46, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-391.png)<!-- -->

    ## Rows: 275 Columns: 87
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (86): MER45B_te_264, MER45B_te_266, MER45B_te_272, MER45B_te_274, MER45B...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-392.png)<!-- -->

    ## Rows: 275 Columns: 111
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (110): MER45R_te_513, MER45R_te_515, MER45R_te_516, MER45R_te_522, MER45...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-393.png)<!-- -->

    ## Rows: 275 Columns: 63
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (62): MER47B_te_35, MER47B_te_36, MER47B_te_47, MER47B_te_48, MER47B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-394.png)<!-- -->

    ## Rows: 275 Columns: 49
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (48): MER48_te_16, MER48_te_31, MER48_te_33, MER48_te_35, MER48_te_36, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-395.png)<!-- -->

    ## Rows: 275 Columns: 72
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (71): MER49_te_21, MER49_te_26, MER49_te_31, MER49_te_32, MER49_te_34, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-396.png)<!-- -->

    ## Rows: 275 Columns: 54
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (53): MER4A_te_21, MER4A_te_22, MER4A_te_30, MER4A_te_33, MER4A_te_70, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-397.png)<!-- -->

    ## Rows: 275 Columns: 33
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (32): MER4A1_LTR_te_21, MER4A1_LTR_te_22, MER4A1_LTR_te_30, MER4A1_LTR_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-398.png)<!-- -->

    ## Rows: 275 Columns: 35
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (34): MER4A1_te_30, MER4A1_te_107, MER4A1_te_127, MER4A1_te_154, MER4A1_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-399.png)<!-- -->

    ## Rows: 275 Columns: 60
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (59): MER4B_te_36, MER4B_te_37, MER4B_te_57, MER4B_te_58, MER4B_te_120, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-400.png)<!-- -->

    ## Rows: 275 Columns: 571
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (570): MER4BI_te_6, MER4BI_te_10, MER4BI_te_11, MER4BI_te_20, MER4BI_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-401.png)<!-- -->

    ## Rows: 275 Columns: 24
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (23): MER4C_te_22, MER4C_te_23, MER4C_te_39, MER4C_te_49, MER4C_te_84, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-402.png)<!-- -->

    ## Rows: 275 Columns: 70
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (69): MER4CL34_te_15, MER4CL34_te_21, MER4CL34_te_22, MER4CL34_te_76, ME...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-403.png)<!-- -->

    ## Rows: 275 Columns: 65
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (64): MER4D_LTR_te_35, MER4D_LTR_te_36, MER4D_LTR_te_37, MER4D_LTR_te_58...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-404.png)<!-- -->

    ## Rows: 275 Columns: 42
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (41): MER4D_te_29, MER4D_te_30, MER4D_te_35, MER4D_te_36, MER4D_te_48, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-405.png)<!-- -->

    ## Rows: 275 Columns: 31
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (30): MER4D1_te_37, MER4D1_te_57, MER4D1_te_58, MER4D1_te_78, MER4D1_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-406.png)<!-- -->

    ## Rows: 275 Columns: 35
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (34): MER4E_te_43, MER4E_te_55, MER4E_te_105, MER4E_te_111, MER4E_te_117...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-407.png)<!-- -->

    ## Rows: 275 Columns: 35
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (34): MER4E1_te_76, MER4E1_te_118, MER4E1_te_119, MER4E1_te_135, MER4E1_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-408.png)<!-- -->

    ## Rows: 275 Columns: 418
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (417): MER4I_te_31, MER4I_te_60, MER4I_te_61, MER4I_te_77, MER4I_te_99, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-409.png)<!-- -->

    ## Rows: 275 Columns: 75
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (74): MER50_te_6, MER50_te_20, MER50_te_24, MER50_te_58, MER50_te_75, ME...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-410.png)<!-- -->

    ## Rows: 275 Columns: 80
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (79): MER50B_te_24, MER50B_te_25, MER50B_te_43, MER50B_te_58, MER50B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-411.png)<!-- -->

    ## Rows: 275 Columns: 152
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (151): MER50C_te_6, MER50C_te_9, MER50C_te_15, MER50C_te_16, MER50C_te_2...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-412.png)<!-- -->

    ## Rows: 275 Columns: 881
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (880): MER50I_te_12, MER50I_te_13, MER50I_te_14, MER50I_te_21, MER50I_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-413.png)<!-- -->

    ## Rows: 275 Columns: 60
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (59): MER51A_te_42, MER51A_te_43, MER51A_te_44, MER51A_te_49, MER51A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-414.png)<!-- -->

    ## Rows: 275 Columns: 31
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (30): MER51B_te_42, MER51B_te_43, MER51B_te_61, MER51B_te_87, MER51B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-415.png)<!-- -->

    ## Rows: 275 Columns: 57
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (56): MER51C_te_6, MER51C_te_10, MER51C_te_33, MER51C_te_47, MER51C_te_7...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-416.png)<!-- -->

    ## Rows: 275 Columns: 89
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (88): MER51D_te_13, MER51D_te_16, MER51D_te_17, MER51D_te_25, MER51D_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-417.png)<!-- -->

    ## Rows: 275 Columns: 38
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (37): MER51E_te_45, MER51E_te_87, MER51E_te_95, MER51E_te_108, MER51E_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-418.png)<!-- -->

    ## Rows: 275 Columns: 665
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (664): MER51I_te_9, MER51I_te_10, MER51I_te_51, MER51I_te_62, MER51I_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-419.png)<!-- -->

    ## Rows: 275 Columns: 303
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (302): MER52A_te_19, MER52A_te_27, MER52A_te_32, MER52A_te_37, MER52A_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-420.png)<!-- -->

    ## Rows: 275 Columns: 798
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (797): MER52AI_te_5, MER52AI_te_17, MER52AI_te_18, MER52AI_te_27, MER52A...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-421.png)<!-- -->

    ## Rows: 275 Columns: 281
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (280): MER52B_te_10, MER52B_te_11, MER52B_te_14, MER52B_te_19, MER52B_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-422.png)<!-- -->

    ## Rows: 275 Columns: 150
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (149): MER52C_te_25, MER52C_te_28, MER52C_te_45, MER52C_te_46, MER52C_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-423.png)<!-- -->

    ## Rows: 275 Columns: 276
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (275): MER52D_te_31, MER52D_te_32, MER52D_te_37, MER52D_te_40, MER52D_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-424.png)<!-- -->

    ## Rows: 275 Columns: 147
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (146): MER54_te_32, MER54_te_33, MER54_te_53, MER54_te_56, MER54_te_83, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-425.png)<!-- -->

    ## Rows: 275 Columns: 102
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (101): MER54A_te_12, MER54A_te_15, MER54A_te_29, MER54A_te_32, MER54A_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-426.png)<!-- -->

    ## Rows: 275 Columns: 78
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (77): MER54B_te_30, MER54B_te_32, MER54B_te_33, MER54B_te_52, MER54B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-427.png)<!-- -->

    ## Rows: 275 Columns: 491
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (490): MER57A_I_te_6, MER57A_I_te_7, MER57A_I_te_22, MER57A_I_te_23, MER...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-428.png)<!-- -->

    ## Rows: 275 Columns: 38
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (37): MER57A1_te_10, MER57A1_te_37, MER57A1_te_38, MER57A1_te_65, MER57A...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-429.png)<!-- -->

    ## Rows: 275 Columns: 46
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (45): MER57B1_te_49, MER57B1_te_50, MER57B1_te_63, MER57B1_te_73, MER57B...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-430.png)<!-- -->

    ## Rows: 275 Columns: 47
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (46): MER57B2_te_49, MER57B2_te_50, MER57B2_te_66, MER57B2_te_67, MER57B...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-431.png)<!-- -->

    ## Rows: 275 Columns: 32
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (31): MER57C1_te_19, MER57C1_te_24, MER57C1_te_32, MER57C1_te_33, MER57C...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-432.png)<!-- -->

    ## Rows: 275 Columns: 22
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (21): MER57C2_te_24, MER57C2_te_36, MER57C2_te_37, MER57C2_te_38, MER57C...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-433.png)<!-- -->

    ## Rows: 275 Columns: 33
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (32): MER57D_te_45, MER57D_te_46, MER57D_te_73, MER57D_te_89, MER57D_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-434.png)<!-- -->

    ## Rows: 275 Columns: 54
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (53): MER57E1_te_24, MER57E1_te_31, MER57E1_te_37, MER57E1_te_38, MER57E...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-435.png)<!-- -->

    ## Rows: 275 Columns: 54
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (53): MER57E3_te_15, MER57E3_te_22, MER57E3_te_36, MER57E3_te_37, MER57E...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-436.png)<!-- -->

    ## Rows: 275 Columns: 35
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (34): MER57F_te_55, MER57F_te_80, MER57F_te_107, MER57F_te_108, MER57F_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-437.png)<!-- -->

    ## Rows: 275 Columns: 578
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (577): MER57I_te_6, MER57I_te_15, MER57I_te_23, MER57I_te_27, MER57I_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-438.png)<!-- -->

    ## Rows: 275 Columns: 46
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (45): MER58B_te_23, MER58B_te_24, MER58B_te_27, MER58B_te_36, MER58B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-439.png)<!-- -->

    ## Rows: 275 Columns: 32
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (31): MER58D_te_103, MER58D_te_115, MER58D_te_121, MER58D_te_131, MER58D...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-440.png)<!-- -->

    ## Rows: 275 Columns: 8
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): familyname_position
    ## dbl (7): MER5A_te_25, MER5A_te_26, MER5A_te_70, MER5A_te_76, MER5A_te_77, ME...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-441.png)<!-- -->

    ## Rows: 275 Columns: 5
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): familyname_position
    ## dbl (4): MER5A1_te_22, MER5A1_te_23, MER5A1_te_74, MER5A1_te_143
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-442.png)<!-- -->

    ## Rows: 275 Columns: 16
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (15): MER5B_te_5, MER5B_te_15, MER5B_te_16, MER5B_te_25, MER5B_te_30, ME...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-443.png)<!-- -->

    ## Rows: 275 Columns: 18
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (17): MER5C_te_68, MER5C_te_78, MER5C_te_164, MER5C_te_168, MER5C_te_181...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-444.png)<!-- -->

    ## Rows: 275 Columns: 65
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (64): MER6_te_203, MER6_te_206, MER6_te_208, MER6_te_212, MER6_te_218, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-445.png)<!-- -->

    ## Rows: 275 Columns: 46
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (45): MER61A_te_19, MER61A_te_28, MER61A_te_29, MER61A_te_50, MER61A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-446.png)<!-- -->

    ## Rows: 275 Columns: 56
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (55): MER61B_te_16, MER61B_te_32, MER61B_te_33, MER61B_te_35, MER61B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-447.png)<!-- -->

    ## Rows: 275 Columns: 45
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (44): MER61C_te_26, MER61C_te_35, MER61C_te_36, MER61C_te_46, MER61C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-448.png)<!-- -->

    ## Rows: 275 Columns: 58
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (57): MER61D_te_51, MER61D_te_52, MER61D_te_62, MER61D_te_63, MER61D_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-449.png)<!-- -->

    ## Rows: 275 Columns: 57
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (56): MER61E_te_24, MER61E_te_25, MER61E_te_53, MER61E_te_54, MER61E_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-450.png)<!-- -->

    ## Rows: 275 Columns: 64
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (63): MER61F_te_37, MER61F_te_42, MER61F_te_50, MER61F_te_56, MER61F_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-451.png)<!-- -->

    ## Rows: 275 Columns: 422
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (421): MER61I_te_16, MER61I_te_26, MER61I_te_29, MER61I_te_31, MER61I_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-452.png)<!-- -->

    ## Rows: 275 Columns: 36
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (35): MER65A_te_7, MER65A_te_8, MER65A_te_25, MER65A_te_34, MER65A_te_35...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-453.png)<!-- -->

    ## Rows: 275 Columns: 54
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (53): MER65C_te_14, MER65C_te_18, MER65C_te_19, MER65C_te_53, MER65C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-454.png)<!-- -->

    ## Rows: 275 Columns: 59
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (58): MER65D_te_44, MER65D_te_48, MER65D_te_54, MER65D_te_56, MER65D_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-455.png)<!-- -->

    ## Rows: 275 Columns: 908
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (907): MER66_I_te_100, MER66_I_te_110, MER66_I_te_178, MER66_I_te_210, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-456.png)<!-- -->

    ## Rows: 275 Columns: 65
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (64): MER66A_te_5, MER66A_te_19, MER66A_te_20, MER66A_te_21, MER66A_te_2...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-457.png)<!-- -->

    ## Rows: 275 Columns: 61
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (60): MER66B_te_20, MER66B_te_21, MER66B_te_26, MER66B_te_27, MER66B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-458.png)<!-- -->

    ## Rows: 275 Columns: 43
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (42): MER66C_te_14, MER66C_te_15, MER66C_te_85, MER66C_te_90, MER66C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-459.png)<!-- -->

    ## Rows: 275 Columns: 52
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (51): MER66D_te_31, MER66D_te_32, MER66D_te_36, MER66D_te_37, MER66D_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-460.png)<!-- -->

    ## Rows: 275 Columns: 52
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (51): MER67A_te_51, MER67A_te_57, MER67A_te_58, MER67A_te_105, MER67A_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-461.png)<!-- -->

    ## Rows: 275 Columns: 63
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (62): MER67B_te_19, MER67B_te_20, MER67B_te_22, MER67B_te_24, MER67B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-462.png)<!-- -->

    ## Rows: 275 Columns: 35
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (34): MER67C_te_9, MER67C_te_25, MER67C_te_68, MER67C_te_73, MER67C_te_1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-463.png)<!-- -->

    ## Rows: 275 Columns: 54
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (53): MER67D_te_19, MER67D_te_20, MER67D_te_51, MER67D_te_52, MER67D_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-464.png)<!-- -->

    ## Rows: 275 Columns: 42
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (41): MER68A_te_60, MER68A_te_100, MER68A_te_101, MER68A_te_123, MER68A_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-465.png)<!-- -->

    ## Rows: 275 Columns: 44
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (43): MER68B_te_43, MER68B_te_44, MER68B_te_47, MER68B_te_51, MER68B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-466.png)<!-- -->

    ## Rows: 275 Columns: 48
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (47): MER6A_te_164, MER6A_te_170, MER6A_te_172, MER6A_te_173, MER6A_te_1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-467.png)<!-- -->

    ## Rows: 275 Columns: 45
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (44): MER6B_te_21, MER6B_te_23, MER6B_te_25, MER6B_te_26, MER6B_te_34, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-468.png)<!-- -->

    ## Rows: 275 Columns: 37
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (36): MER6C_te_9, MER6C_te_14, MER6C_te_19, MER6C_te_20, MER6C_te_23, ME...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-469.png)<!-- -->

    ## Rows: 275 Columns: 69
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (68): MER70A_te_48, MER70A_te_54, MER70A_te_55, MER70A_te_56, MER70A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-470.png)<!-- -->

    ## Rows: 275 Columns: 74
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (73): MER70C_te_11, MER70C_te_12, MER70C_te_20, MER70C_te_28, MER70C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-471.png)<!-- -->

    ## Rows: 275 Columns: 65
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (64): MER72_te_20, MER72_te_44, MER72_te_49, MER72_te_83, MER72_te_86, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-472.png)<!-- -->

    ## Rows: 275 Columns: 77
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (76): MER72B_te_20, MER72B_te_28, MER72B_te_33, MER72B_te_35, MER72B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-473.png)<!-- -->

    ## Rows: 275 Columns: 63
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (62): MER73_te_33, MER73_te_34, MER73_te_63, MER73_te_77, MER73_te_78, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-474.png)<!-- -->

    ## Rows: 275 Columns: 79
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (78): MER74_te_71, MER74_te_77, MER74_te_78, MER74_te_79, MER74_te_92, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-475.png)<!-- -->

    ## Rows: 275 Columns: 53
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (52): MER74A_te_9, MER74A_te_15, MER74A_te_33, MER74A_te_71, MER74A_te_8...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-476.png)<!-- -->

    ## Rows: 275 Columns: 71
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (70): MER74C_te_51, MER74C_te_52, MER74C_te_77, MER74C_te_83, MER74C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-477.png)<!-- -->

    ## Rows: 275 Columns: 61
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (60): MER75_te_10, MER75_te_11, MER75_te_19, MER75_te_20, MER75_te_29, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-478.png)<!-- -->

    ## Rows: 275 Columns: 11
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (10): MER75B_te_10, MER75B_te_11, MER75B_te_19, MER75B_te_20, MER75B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-479.png)<!-- -->

    ## Rows: 275 Columns: 52
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (51): MER76_te_20, MER76_te_31, MER76_te_50, MER76_te_55, MER76_te_92, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-480.png)<!-- -->

    ## Rows: 275 Columns: 73
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (72): MER77_te_34, MER77_te_37, MER77_te_64, MER77_te_65, MER77_te_77, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-481.png)<!-- -->

    ## Rows: 275 Columns: 75
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (74): MER82_te_46, MER82_te_47, MER82_te_48, MER82_te_49, MER82_te_51, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-482.png)<!-- -->

    ## Rows: 275 Columns: 45
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (44): MER83_te_37, MER83_te_38, MER83_te_42, MER83_te_45, MER83_te_57, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-483.png)<!-- -->

    ## Rows: 275 Columns: 35
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (34): MER83B_te_8, MER83B_te_29, MER83B_te_38, MER83B_te_47, MER83B_te_5...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-484.png)<!-- -->

    ## Rows: 275 Columns: 630
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (629): MER83BI_te_72, MER83BI_te_73, MER83BI_te_74, MER83BI_te_90, MER83...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-485.png)<!-- -->

    ## Rows: 275 Columns: 39
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (38): MER83C_te_12, MER83C_te_33, MER83C_te_34, MER83C_te_44, MER83C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-486.png)<!-- -->

    ## Rows: 275 Columns: 65
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (64): MER84_te_16, MER84_te_34, MER84_te_38, MER84_te_44, MER84_te_45, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-487.png)<!-- -->

    ## Rows: 275 Columns: 661
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (660): MER84I_te_9, MER84I_te_16, MER84I_te_25, MER84I_te_48, MER84I_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-488.png)<!-- -->

    ## Rows: 275 Columns: 41
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (40): MER87B_te_76, MER87B_te_77, MER87B_te_110, MER87B_te_139, MER87B_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-489.png)<!-- -->

    ## Rows: 275 Columns: 67
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (66): MER88_te_24, MER88_te_25, MER88_te_67, MER88_te_73, MER88_te_74, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-490.png)<!-- -->

    ## Rows: 275 Columns: 33
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (32): MER89_te_21, MER89_te_22, MER89_te_44, MER89_te_46, MER89_te_81, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-491.png)<!-- -->

    ## Rows: 275 Columns: 49
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (48): MER9_te_5, MER9_te_79, MER9_te_88, MER9_te_90, MER9_te_97, MER9_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-492.png)<!-- -->

    ## Rows: 275 Columns: 67
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (66): MER90a_LTR_te_17, MER90a_LTR_te_36, MER90a_LTR_te_46, MER90a_LTR_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-493.png)<!-- -->

    ## Rows: 275 Columns: 66
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (65): MER92B_te_14, MER92B_te_32, MER92B_te_62, MER92B_te_65, MER92B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-494.png)<!-- -->

    ## Rows: 275 Columns: 54
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (53): MER95_te_41, MER95_te_43, MER95_te_54, MER95_te_55, MER95_te_58, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-495.png)<!-- -->

    ## Rows: 275 Columns: 21
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (20): MER96B_te_128, MER96B_te_129, MER96B_te_130, MER96B_te_134, MER96B...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-496.png)<!-- -->

    ## Rows: 275 Columns: 43
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (42): MER9B_te_20, MER9B_te_87, MER9B_te_96, MER9B_te_98, MER9B_te_100, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-497.png)<!-- -->

    ## Rows: 275 Columns: 12
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (11): MIR_te_17, MIR_te_21, MIR_te_27, MIR_te_55, MIR_te_56, MIR_te_62, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-498.png)<!-- -->

    ## Rows: 275 Columns: 9
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): familyname_position
    ## dbl (8): MIR3_te_70, MIR3_te_71, MIR3_te_159, MIR3_te_186, MIR3_te_187, MIR3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-499.png)<!-- -->

    ## Rows: 275 Columns: 18
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (17): MIRb_te_48, MIRb_te_49, MIRb_te_53, MIRb_te_62, MIRb_te_63, MIRb_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-500.png)<!-- -->

    ## Rows: 275 Columns: 140
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (139): MLT-int_te_1, MLT-int_te_12, MLT-int_te_13, MLT-int_te_57, MLT-in...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-501.png)<!-- -->

    ## Rows: 275 Columns: 125
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (124): MLT1_I_te_12, MLT1_I_te_13, MLT1_I_te_68, MLT1_I_te_84, MLT1_I_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-502.png)<!-- -->

    ## Rows: 275 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (22): MLT1A0_te_47, MLT1A0_te_48, MLT1A0_te_54, MLT1A0_te_73, MLT1A0_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-503.png)<!-- -->

    ## Rows: 275 Columns: 25
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (24): MLT1A1_te_20, MLT1A1_te_55, MLT1A1_te_56, MLT1A1_te_63, MLT1A1_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-504.png)<!-- -->

    ## Rows: 275 Columns: 13
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (12): MLT1B_te_91, MLT1B_te_213, MLT1B_te_229, MLT1B_te_230, MLT1B_te_23...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-505.png)<!-- -->

    ## Rows: 275 Columns: 21
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (20): MLT1C_te_92, MLT1C_te_93, MLT1C_te_197, MLT1C_te_198, MLT1C_te_212...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-506.png)<!-- -->

    ## Rows: 275 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (22): MLT1C1_te_26, MLT1C1_te_27, MLT1C1_te_69, MLT1C1_te_127, MLT1C1_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-507.png)<!-- -->

    ## Rows: 275 Columns: 17
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (16): MLT1C2_te_39, MLT1C2_te_53, MLT1C2_te_54, MLT1C2_te_220, MLT1C2_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-508.png)<!-- -->

    ## Rows: 275 Columns: 31
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (30): MLT1D_te_11, MLT1D_te_38, MLT1D_te_39, MLT1D_te_50, MLT1D_te_51, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-509.png)<!-- -->

    ## Rows: 275 Columns: 48
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (47): MLT1E_te_19, MLT1E_te_29, MLT1E_te_38, MLT1E_te_39, MLT1E_te_41, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-510.png)<!-- -->

    ## Rows: 275 Columns: 68
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (67): MLT1E1_te_39, MLT1E1_te_56, MLT1E1_te_57, MLT1E1_te_68, MLT1E1_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-511.png)<!-- -->

    ## Rows: 275 Columns: 49
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (48): MLT1E1A_te_64, MLT1E1A_te_76, MLT1E1A_te_81, MLT1E1A_te_84, MLT1E1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-512.png)<!-- -->

    ## Rows: 275 Columns: 47
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (46): MLT1E2_te_37, MLT1E2_te_67, MLT1E2_te_68, MLT1E2_te_86, MLT1E2_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-513.png)<!-- -->

    ## Rows: 275 Columns: 174
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (173): MLT1F_I_te_73, MLT1F_I_te_79, MLT1F_I_te_102, MLT1F_I_te_104, MLT...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-514.png)<!-- -->

    ## Rows: 275 Columns: 45
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (44): MLT1F_te_27, MLT1F_te_31, MLT1F_te_40, MLT1F_te_50, MLT1F_te_80, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-515.png)<!-- -->

    ## Rows: 275 Columns: 37
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (36): MLT1F1_te_44, MLT1F1_te_49, MLT1F1_te_52, MLT1F1_te_70, MLT1F1_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-516.png)<!-- -->

    ## Rows: 275 Columns: 31
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (30): MLT1F2_te_14, MLT1F2_te_175, MLT1F2_te_180, MLT1F2_te_187, MLT1F2_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-517.png)<!-- -->

    ## Rows: 275 Columns: 75
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (74): MLT1G1_te_10, MLT1G1_te_30, MLT1G1_te_31, MLT1G1_te_48, MLT1G1_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-518.png)<!-- -->

    ## Rows: 275 Columns: 94
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (93): MLT1G2_te_90, MLT1G2_te_100, MLT1G2_te_102, MLT1G2_te_104, MLT1G2_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-519.png)<!-- -->

    ## Rows: 275 Columns: 27
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (26): MLT1G3_te_192, MLT1G3_te_207, MLT1G3_te_238, MLT1G3_te_239, MLT1G3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-520.png)<!-- -->

    ## Rows: 275 Columns: 52
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (51): MLT1H_te_28, MLT1H_te_30, MLT1H_te_33, MLT1H_te_56, MLT1H_te_98, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-521.png)<!-- -->

    ## Rows: 275 Columns: 55
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (54): MLT1H1_te_29, MLT1H1_te_34, MLT1H1_te_35, MLT1H1_te_38, MLT1H1_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-522.png)<!-- -->

    ## Rows: 275 Columns: 37
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (36): MLT1H2_te_27, MLT1H2_te_34, MLT1H2_te_40, MLT1H2_te_47, MLT1H2_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-523.png)<!-- -->

    ## Rows: 275 Columns: 9
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): familyname_position
    ## dbl (8): MLT1I_te_143, MLT1I_te_144, MLT1I_te_156, MLT1I_te_159, MLT1I_te_16...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-524.png)<!-- -->

    ## Rows: 275 Columns: 19
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (18): MLT1J2_te_50, MLT1J2_te_101, MLT1J2_te_135, MLT1J2_te_148, MLT1J2_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-525.png)<!-- -->

    ## Rows: 275 Columns: 28
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (27): MLT2A1_te_88, MLT2A1_te_128, MLT2A1_te_152, MLT2A1_te_170, MLT2A1_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-526.png)<!-- -->

    ## Rows: 275 Columns: 41
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (40): MLT2A2_te_15, MLT2A2_te_56, MLT2A2_te_112, MLT2A2_te_119, MLT2A2_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-527.png)<!-- -->

    ## Rows: 275 Columns: 42
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (41): MLT2B2_te_60, MLT2B2_te_97, MLT2B2_te_98, MLT2B2_te_101, MLT2B2_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-528.png)<!-- -->

    ## Rows: 275 Columns: 58
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (57): MLT2B3_te_11, MLT2B3_te_17, MLT2B3_te_101, MLT2B3_te_123, MLT2B3_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-529.png)<!-- -->

    ## Rows: 275 Columns: 33
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (32): MLT2B4_te_39, MLT2B4_te_40, MLT2B4_te_60, MLT2B4_te_79, MLT2B4_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-530.png)<!-- -->

    ## Rows: 275 Columns: 41
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (40): MLT2C2_te_7, MLT2C2_te_37, MLT2C2_te_43, MLT2C2_te_44, MLT2C2_te_4...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-531.png)<!-- -->

    ## Rows: 275 Columns: 36
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (35): MLT2D_te_18, MLT2D_te_20, MLT2D_te_25, MLT2D_te_33, MLT2D_te_37, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-532.png)<!-- -->

    ## Rows: 275 Columns: 17
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (16): MSR1_te_54, MSR1_te_65, MSR1_te_67, MSR1_te_71, MSR1_te_77, MSR1_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-533.png)<!-- -->

    ## Rows: 275 Columns: 157
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (156): MST_I_te_23, MST_I_te_57, MST_I_te_58, MST_I_te_59, MST_I_te_76, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-534.png)<!-- -->

    ## Rows: 275 Columns: 31
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (30): MSTA_te_14, MSTA_te_41, MSTA_te_57, MSTA_te_67, MSTA_te_139, MSTA_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-535.png)<!-- -->

    ## Rows: 275 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (22): MSTA1_te_58, MSTA1_te_143, MSTA1_te_145, MSTA1_te_146, MSTA1_te_14...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-536.png)<!-- -->

    ## Rows: 275 Columns: 25
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (24): MSTA2_te_63, MSTA2_te_69, MSTA2_te_70, MSTA2_te_145, MSTA2_te_155,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-537.png)<!-- -->

    ## Rows: 275 Columns: 44
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (43): MSTB_te_48, MSTB_te_56, MSTB_te_90, MSTB_te_100, MSTB_te_101, MSTB...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-538.png)<!-- -->

    ## Rows: 275 Columns: 24
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (23): MSTB1_te_32, MSTB1_te_56, MSTB1_te_100, MSTB1_te_101, MSTB1_te_133...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-539.png)<!-- -->

    ## Rows: 275 Columns: 24
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (23): MSTC_te_16, MSTC_te_45, MSTC_te_98, MSTC_te_150, MSTC_te_182, MSTC...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-540.png)<!-- -->

    ## Rows: 275 Columns: 39
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (38): MSTD_te_69, MSTD_te_78, MSTD_te_83, MSTD_te_85, MSTD_te_133, MSTD_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-541.png)<!-- -->

    ## Rows: 275 Columns: 19
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (18): ORSL_te_13, ORSL_te_14, ORSL_te_22, ORSL_te_31, ORSL_te_50, ORSL_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-542.png)<!-- -->

    ## Rows: 275 Columns: 32
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (31): PABL_A_te_24, PABL_A_te_25, PABL_A_te_51, PABL_A_te_88, PABL_A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-543.png)<!-- -->

    ## Rows: 275 Columns: 313
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (312): PABL_AI_te_11, PABL_AI_te_13, PABL_AI_te_17, PABL_AI_te_21, PABL_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-544.png)<!-- -->

    ## Rows: 275 Columns: 43
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (42): PABL_B_te_47, PABL_B_te_98, PABL_B_te_109, PABL_B_te_158, PABL_B_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-545.png)<!-- -->

    ## Rows: 275 Columns: 990
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (989): PRIMA4_I_te_12, PRIMA4_I_te_13, PRIMA4_I_te_14, PRIMA4_I_te_31, P...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-546.png)<!-- -->

    ## Rows: 275 Columns: 34
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (33): PRIMA4_LTR_te_19, PRIMA4_LTR_te_39, PRIMA4_LTR_te_75, PRIMA4_LTR_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-547.png)<!-- -->

    ## Rows: 275 Columns: 751
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (750): PRIMA41_te_10, PRIMA41_te_11, PRIMA41_te_17, PRIMA41_te_18, PRIMA...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-548.png)<!-- -->

    ## Rows: 275 Columns: 58
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (57): PrimLTR79_te_14, PrimLTR79_te_24, PrimLTR79_te_33, PrimLTR79_te_58...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-549.png)<!-- -->

    ## Rows: 275 Columns: 87
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (86): PTR5_te_14, PTR5_te_43, PTR5_te_48, PTR5_te_51, PTR5_te_52, PTR5_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-550.png)<!-- -->

    ## Rows: 275 Columns: 220
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (219): RICKSHA_0_te_34, RICKSHA_0_te_35, RICKSHA_0_te_38, RICKSHA_0_te_3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-551.png)<!-- -->

    ## Rows: 275 Columns: 179
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (178): Ricksha_a_te_17, Ricksha_a_te_19, Ricksha_a_te_20, Ricksha_a_te_2...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-552.png)<!-- -->

    ## Rows: 275 Columns: 224
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (223): RICKSHA_te_18, RICKSHA_te_19, RICKSHA_te_20, RICKSHA_te_25, RICKS...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-553.png)<!-- -->

    ## Rows: 275 Columns: 39
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (38): SATR1_te_21, SATR1_te_30, SATR1_te_35, SATR1_te_77, SATR1_te_106, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-554.png)<!-- -->

    ## Rows: 275 Columns: 43
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (42): SATR2_te_25, SATR2_te_35, SATR2_te_41, SATR2_te_43, SATR2_te_44, S...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-555.png)<!-- -->

    ## Rows: 275 Columns: 21
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (20): SN5_te_41, SN5_te_44, SN5_te_51, SN5_te_140, SN5_te_142, SN5_te_18...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-556.png)<!-- -->

    ## Rows: 275 Columns: 42
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (41): SVA_A_te_92, SVA_A_te_155, SVA_A_te_166, SVA_A_te_251, SVA_A_te_26...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-557.png)<!-- -->

    ## Rows: 275 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (22): SVA2_te_2, SVA2_te_13, SVA2_te_32, SVA2_te_36, SVA2_te_40, SVA2_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-558.png)<!-- -->

    ## Rows: 275 Columns: 157
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (156): TAR1_te_27, TAR1_te_44, TAR1_te_49, TAR1_te_53, TAR1_te_71, TAR1_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-559.png)<!-- -->

    ## Rows: 275 Columns: 123
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (122): THE1_I_te_26, THE1_I_te_27, THE1_I_te_28, THE1_I_te_34, THE1_I_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-560.png)<!-- -->

    ## Rows: 275 Columns: 27
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (26): THE1A_te_64, THE1A_te_65, THE1A_te_110, THE1A_te_111, THE1A_te_120...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-561.png)<!-- -->

    ## Rows: 275 Columns: 33
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (32): THE1B_te_64, THE1B_te_65, THE1B_te_69, THE1B_te_70, THE1B_te_82, T...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-562.png)<!-- -->

    ## Rows: 275 Columns: 33
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (32): THE1C_te_64, THE1C_te_65, THE1C_te_69, THE1C_te_70, THE1C_te_82, T...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-563.png)<!-- -->

    ## Rows: 275 Columns: 27
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (26): THE1D_te_40, THE1D_te_65, THE1D_te_66, THE1D_te_71, THE1D_te_78, T...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-564.png)<!-- -->

    ## Rows: 275 Columns: 14
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (13): THER1_te_75, THER1_te_84, THER1_te_90, THER1_te_96, THER1_te_98, T...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-565.png)<!-- -->

    ## Rows: 275 Columns: 127
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (126): TIGGER1_te_12, TIGGER1_te_13, TIGGER1_te_22, TIGGER1_te_23, TIGGE...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-566.png)<!-- -->

    ## Rows: 275 Columns: 165
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (164): TIGGER2_te_380, TIGGER2_te_393, TIGGER2_te_399, TIGGER2_te_411, T...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-567.png)<!-- -->

    ## Rows: 275 Columns: 54
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (53): Tigger2b_Pri_te_405, Tigger2b_Pri_te_407, Tigger2b_Pri_te_411, Tig...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-568.png)<!-- -->

    ## Rows: 275 Columns: 82
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (81): Tigger3b_te_22, Tigger3b_te_23, Tigger3b_te_27, Tigger3b_te_28, Ti...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-569.png)<!-- -->

    ## Rows: 275 Columns: 43
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (42): Tigger3d_te_23, Tigger3d_te_27, Tigger3d_te_34, Tigger3d_te_35, Ti...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-570.png)<!-- -->

    ## Rows: 275 Columns: 29
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (28): Tigger4a_te_21, Tigger4a_te_22, Tigger4a_te_29, Tigger4a_te_30, Ti...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-571.png)<!-- -->

    ## Rows: 275 Columns: 31
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (30): TIGGER5_A_te_125, TIGGER5_A_te_130, TIGGER5_A_te_132, TIGGER5_A_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-572.png)<!-- -->

    ## Rows: 275 Columns: 36
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (35): TIGGER5A_te_5, TIGGER5A_te_6, TIGGER5A_te_12, TIGGER5A_te_13, TIGG...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-573.png)<!-- -->

    ## Rows: 275 Columns: 185
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (184): TIGGER6A_te_27, TIGGER6A_te_28, TIGGER6A_te_52, TIGGER6A_te_53, T...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-574.png)<!-- -->

    ## Rows: 275 Columns: 269
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (268): TIGGER6B_te_50, TIGGER6B_te_53, TIGGER6B_te_54, TIGGER6B_te_57, T...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-575.png)<!-- -->

    ## Rows: 275 Columns: 221
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (220): TIGGER7_te_291, TIGGER7_te_292, TIGGER7_te_295, TIGGER7_te_297, T...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-576.png)<!-- -->

    ## Rows: 275 Columns: 107
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (106): Tigger9b_te_12, Tigger9b_te_16, Tigger9b_te_19, Tigger9b_te_20, T...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-577.png)<!-- -->

    ## Rows: 275 Columns: 25
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (24): ZOMBI_A_te_21, ZOMBI_A_te_22, ZOMBI_A_te_29, ZOMBI_A_te_30, ZOMBI_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-578.png)<!-- -->

    ## Rows: 275 Columns: 16
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (15): ZOMBI_B_te_29, ZOMBI_B_te_30, ZOMBI_B_te_67, ZOMBI_B_te_116, ZOMBI...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-579.png)<!-- -->

    ## Rows: 275 Columns: 162
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (161): ZOMBI_te_29, ZOMBI_te_30, ZOMBI_te_67, ZOMBI_te_68, ZOMBI_te_78, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-580.png)<!-- -->

    ## Rows: 553 Columns: 398
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (397): 6kbHsap_te_7, 6kbHsap_te_8, 6kbHsap_te_29, 6kbHsap_te_44, 6kbHsap...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-581.png)<!-- -->

    ## Rows: 553 Columns: 25
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (24): ALR__te_15, ALR__te_29, ALR__te_30, ALR__te_31, ALR__te_37, ALR__t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-582.png)<!-- -->

    ## Rows: 553 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): familyname_position
    ## dbl (2): ALR_te_138, ALR_te_157
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-583.png)<!-- -->

    ## Rows: 553 Columns: 18
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (17): ALR1_te_16, ALR1_te_19, ALR1_te_20, ALR1_te_22, ALR1_te_41, ALR1_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-584.png)<!-- -->

    ## Rows: 553 Columns: 11
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (10): ALR2_te_63, ALR2_te_67, ALR2_te_101, ALR2_te_102, ALR2_te_119, ALR...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-585.png)<!-- -->

    ## Rows: 553 Columns: 17
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (16): ALRa__te_13, ALRa__te_14, ALRa__te_21, ALRa__te_32, ALRa__te_43, A...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-586.png)<!-- -->

    ## Rows: 553 Columns: 18
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (17): ALRa_te_15, ALRa_te_30, ALRa_te_48, ALRa_te_65, ALRa_te_83, ALRa_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-587.png)<!-- -->

    ## Rows: 553 Columns: 26
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (25): ALRb_te_6, ALRb_te_8, ALRb_te_15, ALRb_te_16, ALRb_te_19, ALRb_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-588.png)<!-- -->

    ## Rows: 553 Columns: 53
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (52): ALU_te_9, ALU_te_11, ALU_te_20, ALU_te_48, ALU_te_54, ALU_te_57, A...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-589.png)<!-- -->

    ## Rows: 553 Columns: 46
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (45): BSRf_te_13, BSRf_te_20, BSRf_te_30, BSRf_te_34, BSRf_te_57, BSRf_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-590.png)<!-- -->

    ## Rows: 553 Columns: 21
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (20): CER_te_10, CER_te_17, CER_te_25, CER_te_58, CER_te_65, CER_te_67, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-591.png)<!-- -->

    ## Rows: 553 Columns: 415
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (414): CHARLIE1_te_108, CHARLIE1_te_114, CHARLIE1_te_115, CHARLIE1_te_12...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-592.png)<!-- -->

    ## Rows: 553 Columns: 194
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (193): Charlie12_te_17, Charlie12_te_18, Charlie12_te_19, Charlie12_te_2...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-593.png)<!-- -->

    ## Rows: 553 Columns: 218
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (217): CHARLIE1A_te_44, CHARLIE1A_te_48, CHARLIE1A_te_55, CHARLIE1A_te_5...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-594.png)<!-- -->

    ## Rows: 553 Columns: 211
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (210): CHARLIE3_te_18, CHARLIE3_te_19, CHARLIE3_te_25, CHARLIE3_te_26, C...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-595.png)<!-- -->

    ## Rows: 553 Columns: 468
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (467): CHARLIE5_te_4, CHARLIE5_te_16, CHARLIE5_te_17, CHARLIE5_te_67, CH...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-596.png)<!-- -->

    ## Rows: 553 Columns: 33
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (32): CHESHIRE_A_te_8, CHESHIRE_A_te_9, CHESHIRE_A_te_44, CHESHIRE_A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-597.png)<!-- -->

    ## Rows: 553 Columns: 329
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (328): CHESHIRE_te_33, CHESHIRE_te_40, CHESHIRE_te_44, CHESHIRE_te_52, C...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-598.png)<!-- -->

    ## Rows: 553 Columns: 1735
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1734): ERV24_Prim_te_179, ERV24_Prim_te_180, ERV24_Prim_te_185, ERV24_P...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-599.png)<!-- -->

    ## Rows: 553 Columns: 750
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (749): ERVL_te_118, ERVL_te_126, ERVL_te_131, ERVL_te_142, ERVL_te_143, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-600.png)<!-- -->

    ## Rows: 553 Columns: 468
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (467): ERVL-B4_te_12, ERVL-B4_te_13, ERVL-B4_te_44, ERVL-B4_te_110, ERVL...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-601.png)<!-- -->

    ## Rows: 553 Columns: 54
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (53): GOLEM_A_te_9, GOLEM_A_te_10, GOLEM_A_te_19, GOLEM_A_te_20, GOLEM_A...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-602.png)<!-- -->

    ## Rows: 553 Columns: 90
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (89): GOLEM_B_te_9, GOLEM_B_te_13, GOLEM_B_te_19, GOLEM_B_te_20, GOLEM_B...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-603.png)<!-- -->

    ## Rows: 553 Columns: 39
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (38): GOLEM_C_te_22, GOLEM_C_te_27, GOLEM_C_te_28, GOLEM_C_te_35, GOLEM_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-604.png)<!-- -->

    ## Rows: 553 Columns: 207
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (206): GOLEM_te_28, GOLEM_te_34, GOLEM_te_35, GOLEM_te_37, GOLEM_te_38, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-605.png)<!-- -->

    ## Rows: 553 Columns: 14
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (13): GSAT_te_55, GSAT_te_56, GSAT_te_59, GSAT_te_60, GSAT_te_105, GSAT_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-606.png)<!-- -->

    ## Rows: 553 Columns: 14
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (13): GSATII_te_19, GSATII_te_27, GSATII_te_42, GSATII_te_56, GSATII_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-607.png)<!-- -->

    ## Rows: 553 Columns: 26
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (25): GSATX_te_11, GSATX_te_14, GSATX_te_17, GSATX_te_22, GSATX_te_49, G...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-608.png)<!-- -->

    ## Rows: 553 Columns: 410
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (409): HAL1M8_te_107, HAL1M8_te_118, HAL1M8_te_126, HAL1M8_te_127, HAL1M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-609.png)<!-- -->

    ## Rows: 553 Columns: 569
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (568): HARLEQUIN_te_19, HARLEQUIN_te_25, HARLEQUIN_te_35, HARLEQUIN_te_3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-610.png)<!-- -->

    ## Rows: 553 Columns: 49
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (48): HARLEQUINLTR_te_61, HARLEQUINLTR_te_96, HARLEQUINLTR_te_107, HARLE...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-611.png)<!-- -->

    ## Rows: 553 Columns: 526
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (525): HERV-K14CI_te_13, HERV-K14CI_te_14, HERV-K14CI_te_15, HERV-K14CI_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-612.png)<!-- -->

    ## Rows: 553 Columns: 579
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (578): HERV-K14I_te_16, HERV-K14I_te_22, HERV-K14I_te_34, HERV-K14I_te_3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-613.png)<!-- -->

    ## Rows: 553 Columns: 40
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (39): HERV1_LTR_te_44, HERV1_LTR_te_47, HERV1_LTR_te_81, HERV1_LTR_te_84...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-614.png)<!-- -->

    ## Rows: 553 Columns: 1307
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1306): HERV15I_te_6, HERV15I_te_12, HERV15I_te_15, HERV15I_te_17, HERV1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-615.png)<!-- -->

    ## Rows: 553 Columns: 525
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (524): HERV17_te_7, HERV17_te_8, HERV17_te_13, HERV17_te_30, HERV17_te_3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-616.png)<!-- -->

    ## Rows: 553 Columns: 1133
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1132): HERV18_te_7, HERV18_te_17, HERV18_te_18, HERV18_te_30, HERV18_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-617.png)<!-- -->

    ## Rows: 553 Columns: 562
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (561): HERV19I_te_9, HERV19I_te_16, HERV19I_te_17, HERV19I_te_27, HERV19...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-618.png)<!-- -->

    ## Rows: 553 Columns: 1792
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1791): HERV3_te_4, HERV3_te_29, HERV3_te_30, HERV3_te_36, HERV3_te_39, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-619.png)<!-- -->

    ## Rows: 553 Columns: 1202
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1201): HERV30I_te_7, HERV30I_te_22, HERV30I_te_44, HERV30I_te_45, HERV3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-620.png)<!-- -->

    ## Rows: 553 Columns: 883
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (882): HERV35I_te_115, HERV35I_te_118, HERV35I_te_119, HERV35I_te_126, H...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-621.png)<!-- -->

    ## Rows: 553 Columns: 354
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (353): HERV38I_te_88, HERV38I_te_92, HERV38I_te_93, HERV38I_te_98, HERV3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-622.png)<!-- -->

    ## Rows: 553 Columns: 298
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (297): HERV39_te_487, HERV39_te_512, HERV39_te_514, HERV39_te_515, HERV3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-623.png)<!-- -->

    ## Rows: 553 Columns: 765
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (764): HERV4_I_te_9, HERV4_I_te_27, HERV4_I_te_39, HERV4_I_te_54, HERV4_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-624.png)<!-- -->

    ## Rows: 553 Columns: 1210
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1209): HERV46I_te_16, HERV46I_te_17, HERV46I_te_18, HERV46I_te_20, HERV...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-625.png)<!-- -->

    ## Rows: 553 Columns: 756
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (755): HERV49I_te_12, HERV49I_te_14, HERV49I_te_20, HERV49I_te_21, HERV4...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-626.png)<!-- -->

    ## Rows: 553 Columns: 1089
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1088): HERV57I_te_27, HERV57I_te_30, HERV57I_te_33, HERV57I_te_34, HERV...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-627.png)<!-- -->

    ## Rows: 553 Columns: 513
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (512): HERV9_te_8, HERV9_te_13, HERV9_te_22, HERV9_te_26, HERV9_te_32, H...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-628.png)<!-- -->

    ## Rows: 553 Columns: 662
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (661): HERVE_te_16, HERVE_te_18, HERVE_te_19, HERVE_te_26, HERVE_te_31, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-629.png)<!-- -->

    ## Rows: 553 Columns: 1124
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1123): HERVFH19I_te_16, HERVFH19I_te_17, HERVFH19I_te_19, HERVFH19I_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-630.png)<!-- -->

    ## Rows: 553 Columns: 1137
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1136): HERVFH21I_te_10, HERVFH21I_te_15, HERVFH21I_te_18, HERVFH21I_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-631.png)<!-- -->

    ## Rows: 553 Columns: 552
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (551): HERVH_te_10, HERVH_te_17, HERVH_te_21, HERVH_te_61, HERVH_te_72, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-632.png)<!-- -->

    ## Rows: 553 Columns: 620
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (619): HERVH48I_te_9, HERVH48I_te_10, HERVH48I_te_16, HERVH48I_te_17, HE...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-633.png)<!-- -->

    ## Rows: 553 Columns: 824
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (823): HERVIP10F_te_13, HERVIP10F_te_14, HERVIP10F_te_17, HERVIP10F_te_1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-634.png)<!-- -->

    ## Rows: 553 Columns: 387
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (386): HERVIP10FH_te_12, HERVIP10FH_te_16, HERVIP10FH_te_34, HERVIP10FH_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-635.png)<!-- -->

    ## Rows: 553 Columns: 587
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (586): HERVK_te_18, HERVK_te_19, HERVK_te_73, HERVK_te_90, HERVK_te_106,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-636.png)<!-- -->

    ## Rows: 553 Columns: 584
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (583): HERVK11DI_te_7, HERVK11DI_te_13, HERVK11DI_te_14, HERVK11DI_te_33...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-637.png)<!-- -->

    ## Rows: 553 Columns: 799
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (798): HERVK11I_te_6, HERVK11I_te_7, HERVK11I_te_13, HERVK11I_te_41, HER...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-638.png)<!-- -->

    ## Rows: 553 Columns: 709
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (708): HERVK22I_te_11, HERVK22I_te_12, HERVK22I_te_22, HERVK22I_te_23, H...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-639.png)<!-- -->

    ## Rows: 553 Columns: 758
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (757): HERVK3I_te_7, HERVK3I_te_8, HERVK3I_te_12, HERVK3I_te_13, HERVK3I...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-640.png)<!-- -->

    ## Rows: 553 Columns: 405
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (404): HERVK9I_te_23, HERVK9I_te_26, HERVK9I_te_33, HERVK9I_te_34, HERVK...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-641.png)<!-- -->

    ## Rows: 553 Columns: 797
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (796): HERVKC4_te_23, HERVKC4_te_24, HERVKC4_te_33, HERVKC4_te_39, HERVK...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-642.png)<!-- -->

    ## Rows: 553 Columns: 243
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (242): HERVL_te_57, HERVL_te_58, HERVL_te_105, HERVL_te_159, HERVL_te_16...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-643.png)<!-- -->

    ## Rows: 553 Columns: 596
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (595): HERVL66I_te_8, HERVL66I_te_9, HERVL66I_te_13, HERVL66I_te_14, HER...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-644.png)<!-- -->

    ## Rows: 553 Columns: 724
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (723): HERVP71A_I_te_12, HERVP71A_I_te_13, HERVP71A_I_te_14, HERVP71A_I_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-645.png)<!-- -->

    ## Rows: 553 Columns: 1497
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1496): HERVS71_te_13, HERVS71_te_16, HERVS71_te_21, HERVS71_te_27, HERV...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-646.png)<!-- -->

    ## Rows: 553 Columns: 22
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (21): HSATI_te_4, HSATI_te_14, HSATI_te_50, HSATI_te_57, HSATI_te_79, HS...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-647.png)<!-- -->

    ## Rows: 553 Columns: 26
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (25): HSATII_te_9, HSATII_te_12, HSATII_te_18, HSATII_te_29, HSATII_te_3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-648.png)<!-- -->

    ## Rows: 553 Columns: 150
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (149): HSMAR1_te_24, HSMAR1_te_25, HSMAR1_te_41, HSMAR1_te_49, HSMAR1_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-649.png)<!-- -->

    ## Rows: 553 Columns: 69
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (68): HSMAR2_te_32, HSMAR2_te_33, HSMAR2_te_66, HSMAR2_te_88, HSMAR2_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-650.png)<!-- -->

    ## Rows: 553 Columns: 390
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (389): HUERS-P1_te_14, HUERS-P1_te_17, HUERS-P1_te_18, HUERS-P1_te_43, H...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-651.png)<!-- -->

    ## Rows: 553 Columns: 323
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (322): HUERS-P2_te_54, HUERS-P2_te_55, HUERS-P2_te_59, HUERS-P2_te_70, H...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-652.png)<!-- -->

    ## Rows: 553 Columns: 1154
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1153): HUERS-P3_te_14, HUERS-P3_te_17, HUERS-P3_te_41, HUERS-P3_te_42, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-653.png)<!-- -->

    ## Rows: 553 Columns: 624
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (623): HUERS-P3B_te_50, HUERS-P3B_te_81, HUERS-P3B_te_86, HUERS-P3B_te_8...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-654.png)<!-- -->

    ## Rows: 553 Columns: 85
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (84): IN25_te_5, IN25_te_11, IN25_te_16, IN25_te_43, IN25_te_44, IN25_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-655.png)<!-- -->

    ## Rows: 553 Columns: 273
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (272): L1_te_6, L1_te_21, L1_te_22, L1_te_54, L1_te_55, L1_te_60, L1_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-656.png)<!-- -->

    ## Rows: 553 Columns: 80
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (79): L1HS_te_145, L1HS_te_151, L1HS_te_155, L1HS_te_199, L1HS_te_292, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-657.png)<!-- -->

    ## Rows: 553 Columns: 260
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (259): L1M1_5_te_13, L1M1_5_te_21, L1M1_5_te_22, L1M1_5_te_35, L1M1_5_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-658.png)<!-- -->

    ## Rows: 553 Columns: 112
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (111): L1M1B_5_te_9, L1M1B_5_te_10, L1M1B_5_te_23, L1M1B_5_te_24, L1M1B_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-659.png)<!-- -->

    ## Rows: 553 Columns: 406
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (405): L1M2_5_te_42, L1M2_5_te_43, L1M2_5_te_72, L1M2_5_te_73, L1M2_5_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-660.png)<!-- -->

    ## Rows: 553 Columns: 306
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (305): L1M2A_5_te_25, L1M2A_5_te_29, L1M2A_5_te_34, L1M2A_5_te_39, L1M2A...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-661.png)<!-- -->

    ## Rows: 553 Columns: 267
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (266): L1M2A1_5_te_20, L1M2A1_5_te_21, L1M2A1_5_te_24, L1M2A1_5_te_34, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-662.png)<!-- -->

    ## Rows: 553 Columns: 322
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (321): L1M2B_5_te_9, L1M2B_5_te_11, L1M2B_5_te_12, L1M2B_5_te_37, L1M2B_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-663.png)<!-- -->

    ## Rows: 553 Columns: 331
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (330): L1M2C_5_te_35, L1M2C_5_te_42, L1M2C_5_te_43, L1M2C_5_te_44, L1M2C...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-664.png)<!-- -->

    ## Rows: 553 Columns: 302
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (301): L1M3A_5_te_288, L1M3A_5_te_290, L1M3A_5_te_297, L1M3A_5_te_319, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-665.png)<!-- -->

    ## Rows: 553 Columns: 179
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (178): L1M3C_5_te_26, L1M3C_5_te_32, L1M3C_5_te_98, L1M3C_5_te_99, L1M3C...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-666.png)<!-- -->

    ## Rows: 553 Columns: 264
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (263): L1M3DE_5_te_39, L1M3DE_5_te_41, L1M3DE_5_te_48, L1M3DE_5_te_49, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-667.png)<!-- -->

    ## Rows: 553 Columns: 306
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (305): L1M4B_te_41, L1M4B_te_49, L1M4B_te_56, L1M4B_te_57, L1M4B_te_58, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-668.png)<!-- -->

    ## Rows: 553 Columns: 16
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (15): L1MA1_te_59, L1MA1_te_60, L1MA1_te_79, L1MA1_te_127, L1MA1_te_218,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-669.png)<!-- -->

    ## Rows: 553 Columns: 77
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (76): L1MA10_te_16, L1MA10_te_27, L1MA10_te_36, L1MA10_te_39, L1MA10_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-670.png)<!-- -->

    ## Rows: 553 Columns: 79
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (78): L1MA2_te_6, L1MA2_te_21, L1MA2_te_22, L1MA2_te_24, L1MA2_te_39, L1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-671.png)<!-- -->

    ## Rows: 553 Columns: 33
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (32): L1MA3_te_27, L1MA3_te_37, L1MA3_te_59, L1MA3_te_60, L1MA3_te_127, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-672.png)<!-- -->

    ## Rows: 553 Columns: 47
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (46): L1MA4_te_18, L1MA4_te_21, L1MA4_te_39, L1MA4_te_59, L1MA4_te_60, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-673.png)<!-- -->

    ## Rows: 553 Columns: 29
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (28): L1MA4A_te_41, L1MA4A_te_42, L1MA4A_te_96, L1MA4A_te_162, L1MA4A_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-674.png)<!-- -->

    ## Rows: 553 Columns: 39
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (38): L1MA5_te_39, L1MA5_te_59, L1MA5_te_60, L1MA5_te_64, L1MA5_te_125, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-675.png)<!-- -->

    ## Rows: 553 Columns: 27
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (26): L1MA5A_te_59, L1MA5A_te_64, L1MA5A_te_141, L1MA5A_te_250, L1MA5A_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-676.png)<!-- -->

    ## Rows: 553 Columns: 33
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (32): L1MA6_te_39, L1MA6_te_42, L1MA6_te_45, L1MA6_te_57, L1MA6_te_64, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-677.png)<!-- -->

    ## Rows: 553 Columns: 53
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (52): L1MA7_te_9, L1MA7_te_20, L1MA7_te_21, L1MA7_te_22, L1MA7_te_24, L1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-678.png)<!-- -->

    ## Rows: 553 Columns: 115
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (114): L1MA8_te_59, L1MA8_te_60, L1MA8_te_64, L1MA8_te_87, L1MA8_te_162,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-679.png)<!-- -->

    ## Rows: 553 Columns: 235
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (234): L1MA9_5_te_8, L1MA9_5_te_12, L1MA9_5_te_17, L1MA9_5_te_18, L1MA9_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-680.png)<!-- -->

    ## Rows: 553 Columns: 40
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (39): L1MA9_te_29, L1MA9_te_59, L1MA9_te_60, L1MA9_te_64, L1MA9_te_78, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-681.png)<!-- -->

    ## Rows: 553 Columns: 36
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (35): L1MB1_te_6, L1MB1_te_8, L1MB1_te_28, L1MB1_te_39, L1MB1_te_41, L1M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-682.png)<!-- -->

    ## Rows: 553 Columns: 38
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (37): L1MB2_te_25, L1MB2_te_37, L1MB2_te_45, L1MB2_te_54, L1MB2_te_59, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-683.png)<!-- -->

    ## Rows: 553 Columns: 117
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (116): L1MB3_5_te_12, L1MB3_5_te_15, L1MB3_5_te_39, L1MB3_5_te_135, L1MB...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-684.png)<!-- -->

    ## Rows: 553 Columns: 48
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (47): L1MB3_te_19, L1MB3_te_22, L1MB3_te_25, L1MB3_te_28, L1MB3_te_34, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-685.png)<!-- -->

    ## Rows: 553 Columns: 151
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (150): L1MB4_5_te_12, L1MB4_5_te_27, L1MB4_5_te_35, L1MB4_5_te_36, L1MB4...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-686.png)<!-- -->

    ## Rows: 553 Columns: 27
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (26): L1MB4_te_54, L1MB4_te_64, L1MB4_te_180, L1MB4_te_196, L1MB4_te_206...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-687.png)<!-- -->

    ## Rows: 553 Columns: 41
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (40): L1MB5_te_64, L1MB5_te_65, L1MB5_te_78, L1MB5_te_81, L1MB5_te_91, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-688.png)<!-- -->

    ## Rows: 553 Columns: 175
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (174): L1MB6_5_te_12, L1MB6_5_te_13, L1MB6_5_te_17, L1MB6_5_te_21, L1MB6...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-689.png)<!-- -->

    ## Rows: 553 Columns: 37
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (36): L1MB7_te_78, L1MB7_te_81, L1MB7_te_84, L1MB7_te_87, L1MB7_te_123, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-690.png)<!-- -->

    ## Rows: 553 Columns: 38
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (37): L1MB8_te_77, L1MB8_te_78, L1MB8_te_81, L1MB8_te_82, L1MB8_te_84, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-691.png)<!-- -->

    ## Rows: 553 Columns: 32
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (31): L1MC1_te_59, L1MC1_te_60, L1MC1_te_64, L1MC1_te_65, L1MC1_te_249, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-692.png)<!-- -->

    ## Rows: 553 Columns: 32
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (31): L1MC2_te_50, L1MC2_te_65, L1MC2_te_141, L1MC2_te_192, L1MC2_te_216...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-693.png)<!-- -->

    ## Rows: 553 Columns: 192
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (191): L1MC3_te_45, L1MC3_te_50, L1MC3_te_87, L1MC3_te_93, L1MC3_te_129,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-694.png)<!-- -->

    ## Rows: 553 Columns: 162
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (161): L1MC4_te_84, L1MC4_te_163, L1MC4_te_165, L1MC4_te_172, L1MC4_te_1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-695.png)<!-- -->

    ## Rows: 553 Columns: 112
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (111): L1MCA_5_te_55, L1MCA_5_te_73, L1MCA_5_te_79, L1MCA_5_te_80, L1MCA...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-696.png)<!-- -->

    ## Rows: 553 Columns: 103
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (102): L1MCB_5_te_58, L1MCB_5_te_61, L1MCB_5_te_62, L1MCB_5_te_74, L1MCB...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-697.png)<!-- -->

    ## Rows: 553 Columns: 68
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (67): L1MD1_te_45, L1MD1_te_59, L1MD1_te_65, L1MD1_te_114, L1MD1_te_155,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-698.png)<!-- -->

    ## Rows: 553 Columns: 51
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (50): L1MD2_te_87, L1MD2_te_100, L1MD2_te_200, L1MD2_te_251, L1MD2_te_25...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-699.png)<!-- -->

    ## Rows: 553 Columns: 94
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (93): L1MD3_te_21, L1MD3_te_22, L1MD3_te_25, L1MD3_te_30, L1MD3_te_34, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-700.png)<!-- -->

    ## Rows: 553 Columns: 208
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (207): L1MDA_5_te_24, L1MDA_5_te_30, L1MDA_5_te_47, L1MDA_5_te_48, L1MDA...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-701.png)<!-- -->

    ## Rows: 553 Columns: 100
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (99): L1MDB_5_te_47, L1MDB_5_te_54, L1MDB_5_te_56, L1MDB_5_te_58, L1MDB_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-702.png)<!-- -->

    ## Rows: 553 Columns: 273
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (272): L1ME_ORF2_te_18, L1ME_ORF2_te_20, L1ME_ORF2_te_34, L1ME_ORF2_te_5...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-703.png)<!-- -->

    ## Rows: 553 Columns: 43
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (42): L1ME1_te_76, L1ME1_te_94, L1ME1_te_109, L1ME1_te_194, L1ME1_te_232...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-704.png)<!-- -->

    ## Rows: 553 Columns: 61
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (60): L1ME2_te_167, L1ME2_te_168, L1ME2_te_169, L1ME2_te_172, L1ME2_te_1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-705.png)<!-- -->

    ## Rows: 553 Columns: 84
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (83): L1ME3_te_203, L1ME3_te_206, L1ME3_te_207, L1ME3_te_209, L1ME3_te_2...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-706.png)<!-- -->

    ## Rows: 553 Columns: 101
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (100): L1ME3A_te_20, L1ME3A_te_25, L1ME3A_te_26, L1ME3A_te_28, L1ME3A_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-707.png)<!-- -->

    ## Rows: 553 Columns: 21
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (20): L1ME5_te_255, L1ME5_te_259, L1ME5_te_262, L1ME5_te_264, L1ME5_te_2...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-708.png)<!-- -->

    ## Rows: 553 Columns: 174
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (173): L1MEC_5_te_42, L1MEC_5_te_72, L1MEC_5_te_73, L1MEC_5_te_81, L1MEC...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-709.png)<!-- -->

    ## Rows: 553 Columns: 231
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (230): L1MEf_5end_te_67, L1MEf_5end_te_72, L1MEf_5end_te_76, L1MEf_5end_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-710.png)<!-- -->

    ## Rows: 553 Columns: 201
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (200): L1MEg_5end_te_73, L1MEg_5end_te_74, L1MEg_5end_te_86, L1MEg_5end_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-711.png)<!-- -->

    ## Rows: 553 Columns: 451
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (450): L1P_MA2_te_6, L1P_MA2_te_7, L1P_MA2_te_10, L1P_MA2_te_20, L1P_MA2...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-712.png)<!-- -->

    ## Rows: 553 Columns: 316
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (315): L1P4a_5end_te_15, L1P4a_5end_te_16, L1P4a_5end_te_18, L1P4a_5end_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-713.png)<!-- -->

    ## Rows: 553 Columns: 318
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (317): L1P4c_5end_te_14, L1P4c_5end_te_15, L1P4c_5end_te_21, L1P4c_5end_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-714.png)<!-- -->

    ## Rows: 553 Columns: 370
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (369): L1P4d_5end_te_9, L1P4d_5end_te_12, L1P4d_5end_te_13, L1P4d_5end_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-715.png)<!-- -->

    ## Rows: 553 Columns: 157
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (156): L1P4e_5end_te_14, L1P4e_5end_te_15, L1P4e_5end_te_47, L1P4e_5end_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-716.png)<!-- -->

    ## Rows: 553 Columns: 37
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (36): L1PA10_te_68, L1PA10_te_112, L1PA10_te_116, L1PA10_te_117, L1PA10_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-717.png)<!-- -->

    ## Rows: 553 Columns: 34
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (33): L1PA11_te_12, L1PA11_te_111, L1PA11_te_200, L1PA11_te_201, L1PA11_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-718.png)<!-- -->

    ## Rows: 553 Columns: 355
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (354): L1PA12_5_te_19, L1PA12_5_te_20, L1PA12_5_te_40, L1PA12_5_te_41, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-719.png)<!-- -->

    ## Rows: 553 Columns: 38
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (37): L1PA12_te_9, L1PA12_te_53, L1PA12_te_57, L1PA12_te_71, L1PA12_te_7...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-720.png)<!-- -->

    ## Rows: 553 Columns: 270
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (269): L1PA13_5_te_10, L1PA13_5_te_11, L1PA13_5_te_18, L1PA13_5_te_19, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-721.png)<!-- -->

    ## Rows: 553 Columns: 29
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (28): L1PA13_te_96, L1PA13_te_117, L1PA13_te_250, L1PA13_te_251, L1PA13_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-722.png)<!-- -->

    ## Rows: 553 Columns: 217
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (216): L1PA14_5_te_29, L1PA14_5_te_30, L1PA14_5_te_47, L1PA14_5_te_48, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-723.png)<!-- -->

    ## Rows: 553 Columns: 39
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (38): L1PA14_te_60, L1PA14_te_73, L1PA14_te_96, L1PA14_te_117, L1PA14_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-724.png)<!-- -->

    ## Rows: 553 Columns: 53
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (52): L1PA15_te_69, L1PA15_te_96, L1PA15_te_115, L1PA15_te_117, L1PA15_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-725.png)<!-- -->

    ## Rows: 553 Columns: 354
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (353): L1PA16_5_te_8, L1PA16_5_te_18, L1PA16_5_te_19, L1PA16_5_te_33, L1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-726.png)<!-- -->

    ## Rows: 553 Columns: 53
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (52): L1PA16_te_34, L1PA16_te_54, L1PA16_te_69, L1PA16_te_78, L1PA16_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-727.png)<!-- -->

    ## Rows: 553 Columns: 248
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (247): L1PA17_5_te_26, L1PA17_5_te_27, L1PA17_5_te_29, L1PA17_5_te_30, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-728.png)<!-- -->

    ## Rows: 553 Columns: 19
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (18): L1PA2_te_15, L1PA2_te_16, L1PA2_te_18, L1PA2_te_21, L1PA2_te_25, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-729.png)<!-- -->

    ## Rows: 553 Columns: 18
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (17): L1PA3_te_27, L1PA3_te_28, L1PA3_te_39, L1PA3_te_50, L1PA3_te_61, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-730.png)<!-- -->

    ## Rows: 553 Columns: 22
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (21): L1PA4_te_66, L1PA4_te_74, L1PA4_te_215, L1PA4_te_322, L1PA4_te_331...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-731.png)<!-- -->

    ## Rows: 553 Columns: 33
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (32): L1PA6_te_49, L1PA6_te_129, L1PA6_te_200, L1PA6_te_201, L1PA6_te_21...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-732.png)<!-- -->

    ## Rows: 553 Columns: 179
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (178): L1PA7_5_te_19, L1PA7_5_te_20, L1PA7_5_te_35, L1PA7_5_te_36, L1PA7...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-733.png)<!-- -->

    ## Rows: 553 Columns: 40
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (39): L1PA7_te_84, L1PA7_te_96, L1PA7_te_98, L1PA7_te_116, L1PA7_te_117,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-734.png)<!-- -->

    ## Rows: 553 Columns: 66
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (65): L1PA8_te_21, L1PA8_te_34, L1PA8_te_43, L1PA8_te_48, L1PA8_te_57, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-735.png)<!-- -->

    ## Rows: 553 Columns: 37
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (36): L1PB1_te_54, L1PB1_te_133, L1PB1_te_187, L1PB1_te_188, L1PB1_te_19...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-736.png)<!-- -->

    ## Rows: 553 Columns: 34
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (33): L1PB2_te_18, L1PB2_te_54, L1PB2_te_116, L1PB2_te_117, L1PB2_te_223...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-737.png)<!-- -->

    ## Rows: 553 Columns: 439
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (438): L1PB2c_te_12, L1PB2c_te_13, L1PB2c_te_16, L1PB2c_te_17, L1PB2c_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-738.png)<!-- -->

    ## Rows: 553 Columns: 45
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (44): L1PB3_te_49, L1PB3_te_63, L1PB3_te_87, L1PB3_te_201, L1PB3_te_204,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-739.png)<!-- -->

    ## Rows: 553 Columns: 51
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (50): L1PB4_te_18, L1PB4_te_54, L1PB4_te_115, L1PB4_te_116, L1PB4_te_117...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-740.png)<!-- -->

    ## Rows: 553 Columns: 278
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (277): L1PBA_5_te_10, L1PBA_5_te_11, L1PBA_5_te_45, L1PBA_5_te_46, L1PBA...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-741.png)<!-- -->

    ## Rows: 553 Columns: 50
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (49): L1PBA1_5_te_38, L1PBA1_5_te_39, L1PBA1_5_te_54, L1PBA1_5_te_55, L1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-742.png)<!-- -->

    ## Rows: 553 Columns: 84
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (83): L1PBB_5_te_1, L1PBB_5_te_32, L1PBB_5_te_63, L1PBB_5_te_64, L1PBB_5...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-743.png)<!-- -->

    ## Rows: 553 Columns: 434
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (433): L1PREC1_te_6, L1PREC1_te_21, L1PREC1_te_29, L1PREC1_te_33, L1PREC...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-744.png)<!-- -->

    ## Rows: 553 Columns: 704
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (703): L1PREC2_te_17, L1PREC2_te_18, L1PREC2_te_38, L1PREC2_te_39, L1PRE...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-745.png)<!-- -->

    ## Rows: 553 Columns: 226
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (225): L2_te_219, L2_te_373, L2_te_392, L2_te_488, L2_te_490, L2_te_593,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-746.png)<!-- -->

    ## Rows: 553 Columns: 28
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (27): L2B_te_9, L2B_te_40, L2B_te_51, L2B_te_56, L2B_te_61, L2B_te_64, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-747.png)<!-- -->

    ## Rows: 553 Columns: 37
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (36): LOR1_te_7, LOR1_te_8, LOR1_te_9, LOR1_te_33, LOR1_te_38, LOR1_te_3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-748.png)<!-- -->

    ## Rows: 553 Columns: 37
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (36): LOR1a_LTR_te_85, LOR1a_LTR_te_86, LOR1a_LTR_te_87, LOR1a_LTR_te_97...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-749.png)<!-- -->

    ## Rows: 553 Columns: 47
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (46): LOR1b_LTR_te_102, LOR1b_LTR_te_114, LOR1b_LTR_te_115, LOR1b_LTR_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-750.png)<!-- -->

    ## Rows: 553 Columns: 860
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (859): LOR1I_te_38, LOR1I_te_44, LOR1I_te_47, LOR1I_te_48, LOR1I_te_57, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-751.png)<!-- -->

    ## Rows: 553 Columns: 115
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (114): LSAU_te_21, LSAU_te_22, LSAU_te_32, LSAU_te_35, LSAU_te_74, LSAU_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-752.png)<!-- -->

    ## Rows: 553 Columns: 66
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (65): LTR06_te_25, LTR06_te_26, LTR06_te_57, LTR06_te_58, LTR06_te_60, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-753.png)<!-- -->

    ## Rows: 553 Columns: 81
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (80): LTR1_te_28, LTR1_te_35, LTR1_te_37, LTR1_te_40, LTR1_te_41, LTR1_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-754.png)<!-- -->

    ## Rows: 553 Columns: 64
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (63): LTR10A_te_12, LTR10A_te_24, LTR10A_te_31, LTR10A_te_39, LTR10A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-755.png)<!-- -->

    ## Rows: 553 Columns: 76
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (75): LTR10B_te_19, LTR10B_te_43, LTR10B_te_58, LTR10B_te_67, LTR10B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-756.png)<!-- -->

    ## Rows: 553 Columns: 70
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (69): LTR10B1_te_37, LTR10B1_te_39, LTR10B1_te_40, LTR10B1_te_43, LTR10B...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-757.png)<!-- -->

    ## Rows: 553 Columns: 64
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (63): LTR10B2_te_20, LTR10B2_te_31, LTR10B2_te_49, LTR10B2_te_50, LTR10B...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-758.png)<!-- -->

    ## Rows: 553 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (22): LTR10C_te_46, LTR10C_te_130, LTR10C_te_131, LTR10C_te_254, LTR10C_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-759.png)<!-- -->

    ## Rows: 553 Columns: 40
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (39): LTR10D_te_32, LTR10D_te_33, LTR10D_te_44, LTR10D_te_72, LTR10D_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-760.png)<!-- -->

    ## Rows: 553 Columns: 60
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (59): LTR10E_te_88, LTR10E_te_97, LTR10E_te_162, LTR10E_te_167, LTR10E_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-761.png)<!-- -->

    ## Rows: 553 Columns: 41
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (40): LTR10F_te_26, LTR10F_te_37, LTR10F_te_91, LTR10F_te_107, LTR10F_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-762.png)<!-- -->

    ## Rows: 553 Columns: 56
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (55): LTR10G_te_14, LTR10G_te_18, LTR10G_te_20, LTR10G_te_40, LTR10G_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-763.png)<!-- -->

    ## Rows: 553 Columns: 84
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (83): LTR12_te_21, LTR12_te_32, LTR12_te_33, LTR12_te_34, LTR12_te_67, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-764.png)<!-- -->

    ## Rows: 553 Columns: 74
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (73): LTR12B_te_32, LTR12B_te_33, LTR12B_te_67, LTR12B_te_96, LTR12B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-765.png)<!-- -->

    ## Rows: 553 Columns: 96
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (95): LTR12C_te_25, LTR12C_te_42, LTR12C_te_52, LTR12C_te_105, LTR12C_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-766.png)<!-- -->

    ## Rows: 553 Columns: 134
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (133): LTR12D_te_34, LTR12D_te_39, LTR12D_te_40, LTR12D_te_42, LTR12D_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-767.png)<!-- -->

    ## Rows: 553 Columns: 122
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (121): LTR12E_te_13, LTR12E_te_14, LTR12E_te_20, LTR12E_te_25, LTR12E_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-768.png)<!-- -->

    ## Rows: 553 Columns: 55
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (54): LTR12F_te_32, LTR12F_te_92, LTR12F_te_93, LTR12F_te_105, LTR12F_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-769.png)<!-- -->

    ## Rows: 553 Columns: 91
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (90): LTR13_te_8, LTR13_te_23, LTR13_te_27, LTR13_te_44, LTR13_te_48, LT...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-770.png)<!-- -->

    ## Rows: 553 Columns: 69
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (68): LTR13A_te_9, LTR13A_te_19, LTR13A_te_93, LTR13A_te_143, LTR13A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-771.png)<!-- -->

    ## Rows: 553 Columns: 39
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (38): LTR14_te_12, LTR14_te_30, LTR14_te_122, LTR14_te_132, LTR14_te_144...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-772.png)<!-- -->

    ## Rows: 553 Columns: 26
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (25): LTR14A_te_23, LTR14A_te_99, LTR14A_te_130, LTR14A_te_131, LTR14A_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-773.png)<!-- -->

    ## Rows: 553 Columns: 53
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (52): LTR14B_te_25, LTR14B_te_82, LTR14B_te_83, LTR14B_te_106, LTR14B_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-774.png)<!-- -->

    ## Rows: 553 Columns: 44
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (43): LTR14C_te_23, LTR14C_te_30, LTR14C_te_65, LTR14C_te_75, LTR14C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-775.png)<!-- -->

    ## Rows: 553 Columns: 33
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (32): LTR15_te_40, LTR15_te_80, LTR15_te_127, LTR15_te_185, LTR15_te_221...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-776.png)<!-- -->

    ## Rows: 553 Columns: 40
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (39): LTR16C_te_96, LTR16C_te_97, LTR16C_te_100, LTR16C_te_106, LTR16C_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-777.png)<!-- -->

    ## Rows: 553 Columns: 93
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (92): LTR17_te_32, LTR17_te_33, LTR17_te_43, LTR17_te_68, LTR17_te_69, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-778.png)<!-- -->

    ## Rows: 553 Columns: 63
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (62): LTR18A_te_8, LTR18A_te_18, LTR18A_te_38, LTR18A_te_39, LTR18A_te_4...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-779.png)<!-- -->

    ## Rows: 553 Columns: 82
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (81): LTR18B_te_39, LTR18B_te_40, LTR18B_te_44, LTR18B_te_45, LTR18B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-780.png)<!-- -->

    ## Rows: 553 Columns: 32
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (31): LTR18C_te_56, LTR18C_te_57, LTR18C_te_71, LTR18C_te_82, LTR18C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-781.png)<!-- -->

    ## Rows: 553 Columns: 44
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (43): LTR19A_te_18, LTR19A_te_19, LTR19A_te_37, LTR19A_te_38, LTR19A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-782.png)<!-- -->

    ## Rows: 553 Columns: 28
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (27): LTR19B_te_19, LTR19B_te_26, LTR19B_te_27, LTR19B_te_37, LTR19B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-783.png)<!-- -->

    ## Rows: 553 Columns: 87
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (86): LTR19C_te_13, LTR19C_te_36, LTR19C_te_46, LTR19C_te_50, LTR19C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-784.png)<!-- -->

    ## Rows: 553 Columns: 121
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (120): LTR1A1_te_6, LTR1A1_te_7, LTR1A1_te_42, LTR1A1_te_43, LTR1A1_te_5...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-785.png)<!-- -->

    ## Rows: 553 Columns: 102
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (101): LTR1A2_te_6, LTR1A2_te_27, LTR1A2_te_29, LTR1A2_te_30, LTR1A2_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-786.png)<!-- -->

    ## Rows: 553 Columns: 108
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (107): LTR1B_te_43, LTR1B_te_50, LTR1B_te_51, LTR1B_te_76, LTR1B_te_78, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-787.png)<!-- -->

    ## Rows: 553 Columns: 118
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (117): LTR1B0_te_6, LTR1B0_te_7, LTR1B0_te_21, LTR1B0_te_42, LTR1B0_te_4...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-788.png)<!-- -->

    ## Rows: 553 Columns: 145
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (144): LTR1B1_te_20, LTR1B1_te_43, LTR1B1_te_58, LTR1B1_te_77, LTR1B1_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-789.png)<!-- -->

    ## Rows: 553 Columns: 120
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (119): LTR1C_te_23, LTR1C_te_24, LTR1C_te_31, LTR1C_te_32, LTR1C_te_52, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-790.png)<!-- -->

    ## Rows: 553 Columns: 104
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (103): LTR1C1_te_11, LTR1C1_te_12, LTR1C1_te_23, LTR1C1_te_24, LTR1C1_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-791.png)<!-- -->

    ## Rows: 553 Columns: 103
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (102): LTR1C3_te_11, LTR1C3_te_12, LTR1C3_te_23, LTR1C3_te_24, LTR1C3_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-792.png)<!-- -->

    ## Rows: 553 Columns: 143
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (142): LTR1D_te_42, LTR1D_te_43, LTR1D_te_47, LTR1D_te_51, LTR1D_te_54, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-793.png)<!-- -->

    ## Rows: 553 Columns: 179
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (178): LTR1D1_te_42, LTR1D1_te_43, LTR1D1_te_50, LTR1D1_te_51, LTR1D1_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-794.png)<!-- -->

    ## Rows: 553 Columns: 127
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (126): LTR1E_te_51, LTR1E_te_52, LTR1E_te_54, LTR1E_te_61, LTR1E_te_81, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-795.png)<!-- -->

    ## Rows: 553 Columns: 114
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (113): LTR1F_te_6, LTR1F_te_7, LTR1F_te_42, LTR1F_te_43, LTR1F_te_50, LT...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-796.png)<!-- -->

    ## Rows: 553 Columns: 107
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (106): LTR1F1_te_6, LTR1F1_te_7, LTR1F1_te_42, LTR1F1_te_43, LTR1F1_te_5...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-797.png)<!-- -->

    ## Rows: 553 Columns: 117
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (116): LTR1F2_te_6, LTR1F2_te_7, LTR1F2_te_29, LTR1F2_te_37, LTR1F2_te_4...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-798.png)<!-- -->

    ## Rows: 553 Columns: 55
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (54): LTR2_te_34, LTR2_te_40, LTR2_te_41, LTR2_te_44, LTR2_te_53, LTR2_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-799.png)<!-- -->

    ## Rows: 553 Columns: 45
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (44): LTR21A_te_24, LTR21A_te_43, LTR21A_te_49, LTR21A_te_55, LTR21A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-800.png)<!-- -->

    ## Rows: 553 Columns: 60
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (59): LTR21B_te_11, LTR21B_te_15, LTR21B_te_16, LTR21B_te_18, LTR21B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-801.png)<!-- -->

    ## Rows: 553 Columns: 51
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (50): LTR21C_te_21, LTR21C_te_33, LTR21C_te_34, LTR21C_te_41, LTR21C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-802.png)<!-- -->

    ## Rows: 553 Columns: 38
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (37): LTR22_te_64, LTR22_te_91, LTR22_te_92, LTR22_te_136, LTR22_te_137,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-803.png)<!-- -->

    ## Rows: 553 Columns: 47
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (46): LTR22A_te_11, LTR22A_te_12, LTR22A_te_28, LTR22A_te_46, LTR22A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-804.png)<!-- -->

    ## Rows: 553 Columns: 71
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (70): LTR22B_te_13, LTR22B_te_51, LTR22B_te_56, LTR22B_te_63, LTR22B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-805.png)<!-- -->

    ## Rows: 553 Columns: 73
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (72): LTR22B1_te_7, LTR22B1_te_41, LTR22B1_te_53, LTR22B1_te_57, LTR22B1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-806.png)<!-- -->

    ## Rows: 553 Columns: 43
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (42): LTR22B2_te_56, LTR22B2_te_59, LTR22B2_te_74, LTR22B2_te_78, LTR22B...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-807.png)<!-- -->

    ## Rows: 553 Columns: 93
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (92): LTR22C_te_7, LTR22C_te_9, LTR22C_te_12, LTR22C_te_28, LTR22C_te_34...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-808.png)<!-- -->

    ## Rows: 553 Columns: 72
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (71): LTR22C2_te_7, LTR22C2_te_9, LTR22C2_te_14, LTR22C2_te_18, LTR22C2_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-809.png)<!-- -->

    ## Rows: 553 Columns: 64
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (63): LTR22E_te_5, LTR22E_te_12, LTR22E_te_35, LTR22E_te_59, LTR22E_te_6...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-810.png)<!-- -->

    ## Rows: 553 Columns: 21
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (20): LTR23_te_97, LTR23_te_101, LTR23_te_104, LTR23_te_106, LTR23_te_11...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-811.png)<!-- -->

    ## Rows: 553 Columns: 49
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (48): LTR24_te_4, LTR24_te_50, LTR24_te_77, LTR24_te_78, LTR24_te_89, LT...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-812.png)<!-- -->

    ## Rows: 553 Columns: 37
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (36): LTR24B_te_8, LTR24B_te_78, LTR24B_te_98, LTR24B_te_113, LTR24B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-813.png)<!-- -->

    ## Rows: 553 Columns: 43
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (42): LTR24C_te_25, LTR24C_te_40, LTR24C_te_45, LTR24C_te_61, LTR24C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-814.png)<!-- -->

    ## Rows: 553 Columns: 105
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (104): LTR25_te_13, LTR25_te_15, LTR25_te_17, LTR25_te_33, LTR25_te_49, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-815.png)<!-- -->

    ## Rows: 553 Columns: 1186
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1185): LTR25-int_te_5, LTR25-int_te_8, LTR25-int_te_27, LTR25-int_te_31...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-816.png)<!-- -->

    ## Rows: 553 Columns: 53
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (52): LTR26_te_7, LTR26_te_8, LTR26_te_45, LTR26_te_58, LTR26_te_59, LTR...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-817.png)<!-- -->

    ## Rows: 553 Columns: 63
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (62): LTR26B_te_73, LTR26B_te_80, LTR26B_te_107, LTR26B_te_110, LTR26B_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-818.png)<!-- -->

    ## Rows: 553 Columns: 53
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (52): LTR26C_te_7, LTR26C_te_8, LTR26C_te_15, LTR26C_te_34, LTR26C_te_80...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-819.png)<!-- -->

    ## Rows: 553 Columns: 46
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (45): LTR26D_te_38, LTR26D_te_39, LTR26D_te_89, LTR26D_te_118, LTR26D_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-820.png)<!-- -->

    ## Rows: 553 Columns: 75
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (74): LTR26E_te_12, LTR26E_te_13, LTR26E_te_15, LTR26E_te_20, LTR26E_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-821.png)<!-- -->

    ## Rows: 553 Columns: 69
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (68): LTR27_te_34, LTR27_te_38, LTR27_te_52, LTR27_te_55, LTR27_te_61, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-822.png)<!-- -->

    ## Rows: 553 Columns: 267
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (266): LTR2752_te_7, LTR2752_te_8, LTR2752_te_24, LTR2752_te_33, LTR2752...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-823.png)<!-- -->

    ## Rows: 553 Columns: 74
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (73): LTR27B_te_7, LTR27B_te_43, LTR27B_te_51, LTR27B_te_52, LTR27B_te_1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-824.png)<!-- -->

    ## Rows: 553 Columns: 86
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (85): LTR27C_te_26, LTR27C_te_37, LTR27C_te_43, LTR27C_te_44, LTR27C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-825.png)<!-- -->

    ## Rows: 553 Columns: 100
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (99): LTR27D_te_7, LTR27D_te_8, LTR27D_te_43, LTR27D_te_44, LTR27D_te_51...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-826.png)<!-- -->

    ## Rows: 553 Columns: 126
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (125): LTR27E_te_7, LTR27E_te_8, LTR27E_te_9, LTR27E_te_16, LTR27E_te_42...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-827.png)<!-- -->

    ## Rows: 553 Columns: 112
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (111): LTR28_te_7, LTR28_te_8, LTR28_te_20, LTR28_te_43, LTR28_te_44, LT...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-828.png)<!-- -->

    ## Rows: 553 Columns: 163
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (162): LTR28B_te_42, LTR28B_te_43, LTR28B_te_50, LTR28B_te_51, LTR28B_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-829.png)<!-- -->

    ## Rows: 553 Columns: 190
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (189): LTR28C_te_7, LTR28C_te_8, LTR28C_te_26, LTR28C_te_37, LTR28C_te_6...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-830.png)<!-- -->

    ## Rows: 553 Columns: 50
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (49): LTR29_te_18, LTR29_te_19, LTR29_te_36, LTR29_te_37, LTR29_te_38, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-831.png)<!-- -->

    ## Rows: 553 Columns: 25
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (24): LTR2B_te_91, LTR2B_te_95, LTR2B_te_112, LTR2B_te_124, LTR2B_te_156...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-832.png)<!-- -->

    ## Rows: 553 Columns: 55
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (54): LTR2C_te_16, LTR2C_te_45, LTR2C_te_53, LTR2C_te_67, LTR2C_te_72, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-833.png)<!-- -->

    ## Rows: 553 Columns: 37
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (36): LTR3_te_11, LTR3_te_26, LTR3_te_27, LTR3_te_63, LTR3_te_93, LTR3_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-834.png)<!-- -->

    ## Rows: 553 Columns: 66
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (65): LTR30_te_11, LTR30_te_12, LTR30_te_32, LTR30_te_33, LTR30_te_89, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-835.png)<!-- -->

    ## Rows: 553 Columns: 87
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (86): LTR31_te_10, LTR31_te_11, LTR31_te_16, LTR31_te_32, LTR31_te_33, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-836.png)<!-- -->

    ## Rows: 553 Columns: 51
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (50): LTR32_te_22, LTR32_te_24, LTR32_te_26, LTR32_te_35, LTR32_te_58, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-837.png)<!-- -->

    ## Rows: 553 Columns: 44
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (43): LTR34_te_61, LTR34_te_106, LTR34_te_110, LTR34_te_117, LTR34_te_13...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-838.png)<!-- -->

    ## Rows: 553 Columns: 76
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (75): LTR35_te_14, LTR35_te_32, LTR35_te_38, LTR35_te_40, LTR35_te_42, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-839.png)<!-- -->

    ## Rows: 553 Columns: 38
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (37): LTR35A_te_14, LTR35A_te_17, LTR35A_te_18, LTR35A_te_102, LTR35A_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-840.png)<!-- -->

    ## Rows: 553 Columns: 102
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (101): LTR35B_te_15, LTR35B_te_24, LTR35B_te_26, LTR35B_te_35, LTR35B_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-841.png)<!-- -->

    ## Rows: 553 Columns: 66
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (65): LTR36_te_7, LTR36_te_13, LTR36_te_16, LTR36_te_17, LTR36_te_36, LT...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-842.png)<!-- -->

    ## Rows: 553 Columns: 51
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (50): LTR37A_te_15, LTR37A_te_23, LTR37A_te_30, LTR37A_te_44, LTR37A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-843.png)<!-- -->

    ## Rows: 553 Columns: 71
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (70): LTR37B_te_20, LTR37B_te_27, LTR37B_te_35, LTR37B_te_42, LTR37B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-844.png)<!-- -->

    ## Rows: 553 Columns: 70
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (69): LTR38_te_7, LTR38_te_15, LTR38_te_30, LTR38_te_51, LTR38_te_56, LT...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-845.png)<!-- -->

    ## Rows: 553 Columns: 77
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (76): LTR38A1_te_7, LTR38A1_te_8, LTR38A1_te_23, LTR38A1_te_31, LTR38A1_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-846.png)<!-- -->

    ## Rows: 553 Columns: 50
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (49): LTR38B_te_7, LTR38B_te_40, LTR38B_te_48, LTR38B_te_60, LTR38B_te_1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-847.png)<!-- -->

    ## Rows: 553 Columns: 97
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (96): LTR38C_te_13, LTR38C_te_17, LTR38C_te_29, LTR38C_te_30, LTR38C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-848.png)<!-- -->

    ## Rows: 553 Columns: 75
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (74): LTR39_te_14, LTR39_te_17, LTR39_te_29, LTR39_te_36, LTR39_te_37, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-849.png)<!-- -->

    ## Rows: 553 Columns: 45
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (44): LTR3A_te_12, LTR3A_te_19, LTR3A_te_20, LTR3A_te_99, LTR3A_te_110, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-850.png)<!-- -->

    ## Rows: 553 Columns: 52
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (51): LTR3B_te_6, LTR3B_te_12, LTR3B_te_18, LTR3B_te_27, LTR3B_te_63, LT...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-851.png)<!-- -->

    ## Rows: 553 Columns: 77
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (76): LTR4_te_14, LTR4_te_20, LTR4_te_27, LTR4_te_29, LTR4_te_31, LTR4_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-852.png)<!-- -->

    ## Rows: 553 Columns: 57
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (56): LTR40A_te_30, LTR40A_te_31, LTR40A_te_40, LTR40A_te_42, LTR40A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-853.png)<!-- -->

    ## Rows: 553 Columns: 39
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (38): LTR40B_te_79, LTR40B_te_84, LTR40B_te_90, LTR40B_te_96, LTR40B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-854.png)<!-- -->

    ## Rows: 553 Columns: 52
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (51): LTR41_te_47, LTR41_te_48, LTR41_te_51, LTR41_te_52, LTR41_te_61, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-855.png)<!-- -->

    ## Rows: 553 Columns: 48
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (47): LTR42_te_24, LTR42_te_25, LTR42_te_26, LTR42_te_37, LTR42_te_47, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-856.png)<!-- -->

    ## Rows: 553 Columns: 48
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (47): LTR43_te_28, LTR43_te_36, LTR43_te_39, LTR43_te_40, LTR43_te_44, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-857.png)<!-- -->

    ## Rows: 553 Columns: 69
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (68): LTR43B_te_42, LTR43B_te_55, LTR43B_te_59, LTR43B_te_60, LTR43B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-858.png)<!-- -->

    ## Rows: 553 Columns: 59
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (58): LTR44_te_41, LTR44_te_44, LTR44_te_48, LTR44_te_49, LTR44_te_52, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-859.png)<!-- -->

    ## Rows: 553 Columns: 57
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (56): LTR45_te_7, LTR45_te_8, LTR45_te_9, LTR45_te_47, LTR45_te_54, LTR4...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-860.png)<!-- -->

    ## Rows: 553 Columns: 40
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (39): LTR45B_te_94, LTR45B_te_95, LTR45B_te_124, LTR45B_te_125, LTR45B_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-861.png)<!-- -->

    ## Rows: 553 Columns: 52
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (51): LTR45C_te_29, LTR45C_te_42, LTR45C_te_74, LTR45C_te_105, LTR45C_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-862.png)<!-- -->

    ## Rows: 553 Columns: 48
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (47): LTR46_te_31, LTR46_te_35, LTR46_te_38, LTR46_te_70, LTR46_te_73, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-863.png)<!-- -->

    ## Rows: 553 Columns: 38
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (37): LTR47A_te_38, LTR47A_te_39, LTR47A_te_48, LTR47A_te_64, LTR47A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-864.png)<!-- -->

    ## Rows: 553 Columns: 32
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (31): LTR47A2_te_48, LTR47A2_te_49, LTR47A2_te_54, LTR47A2_te_55, LTR47A...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-865.png)<!-- -->

    ## Rows: 553 Columns: 73
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (72): LTR47B_te_15, LTR47B_te_30, LTR47B_te_38, LTR47B_te_39, LTR47B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-866.png)<!-- -->

    ## Rows: 553 Columns: 47
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (46): LTR47B2_te_38, LTR47B2_te_39, LTR47B2_te_69, LTR47B2_te_70, LTR47B...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-867.png)<!-- -->

    ## Rows: 553 Columns: 66
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (65): LTR47B3_te_38, LTR47B3_te_39, LTR47B3_te_69, LTR47B3_te_70, LTR47B...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-868.png)<!-- -->

    ## Rows: 553 Columns: 49
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (48): LTR47B4_te_38, LTR47B4_te_39, LTR47B4_te_46, LTR47B4_te_62, LTR47B...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-869.png)<!-- -->

    ## Rows: 553 Columns: 66
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (65): LTR48_te_50, LTR48_te_59, LTR48_te_67, LTR48_te_80, LTR48_te_87, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-870.png)<!-- -->

    ## Rows: 553 Columns: 63
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (62): LTR48B_te_14, LTR48B_te_15, LTR48B_te_16, LTR48B_te_52, LTR48B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-871.png)<!-- -->

    ## Rows: 553 Columns: 32
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (31): LTR49_te_11, LTR49_te_21, LTR49_te_35, LTR49_te_61, LTR49_te_69, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-872.png)<!-- -->

    ## Rows: 553 Columns: 46
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (45): LTR5_Hs_te_13, LTR5_Hs_te_79, LTR5_Hs_te_107, LTR5_Hs_te_111, LTR5...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-873.png)<!-- -->

    ## Rows: 553 Columns: 87
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (86): LTR5_te_21, LTR5_te_27, LTR5_te_41, LTR5_te_43, LTR5_te_62, LTR5_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-874.png)<!-- -->

    ## Rows: 553 Columns: 52
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (51): LTR51_te_25, LTR51_te_26, LTR51_te_45, LTR51_te_46, LTR51_te_47, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-875.png)<!-- -->

    ## Rows: 553 Columns: 36
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (35): LTR52_te_64, LTR52_te_72, LTR52_te_77, LTR52_te_112, LTR52_te_130,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-876.png)<!-- -->

    ## Rows: 553 Columns: 63
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (62): LTR53_te_27, LTR53_te_57, LTR53_te_62, LTR53_te_63, LTR53_te_65, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-877.png)<!-- -->

    ## Rows: 553 Columns: 77
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (76): LTR53B_te_11, LTR53B_te_12, LTR53B_te_14, LTR53B_te_16, LTR53B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-878.png)<!-- -->

    ## Rows: 553 Columns: 45
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (44): LTR54_te_13, LTR54_te_30, LTR54_te_39, LTR54_te_97, LTR54_te_106, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-879.png)<!-- -->

    ## Rows: 553 Columns: 36
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (35): LTR54B_te_13, LTR54B_te_78, LTR54B_te_107, LTR54B_te_108, LTR54B_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-880.png)<!-- -->

    ## Rows: 553 Columns: 51
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (50): LTR55_te_16, LTR55_te_20, LTR55_te_33, LTR55_te_64, LTR55_te_68, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-881.png)<!-- -->

    ## Rows: 553 Columns: 21
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (20): LTR56_te_8, LTR56_te_16, LTR56_te_36, LTR56_te_37, LTR56_te_42, LT...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-882.png)<!-- -->

    ## Rows: 553 Columns: 39
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (38): LTR57_te_21, LTR57_te_47, LTR57_te_48, LTR57_te_53, LTR57_te_68, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-883.png)<!-- -->

    ## Rows: 553 Columns: 117
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (116): LTR58_te_8, LTR58_te_10, LTR58_te_14, LTR58_te_20, LTR58_te_26, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-884.png)<!-- -->

    ## Rows: 553 Columns: 58
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (57): LTR59_te_15, LTR59_te_27, LTR59_te_49, LTR59_te_58, LTR59_te_59, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-885.png)<!-- -->

    ## Rows: 553 Columns: 68
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (67): LTR5A_te_8, LTR5A_te_31, LTR5A_te_49, LTR5A_te_52, LTR5A_te_107, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-886.png)<!-- -->

    ## Rows: 553 Columns: 88
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (87): LTR5B_te_50, LTR5B_te_66, LTR5B_te_77, LTR5B_te_85, LTR5B_te_96, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-887.png)<!-- -->

    ## Rows: 553 Columns: 92
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (91): LTR60_te_18, LTR60_te_34, LTR60_te_40, LTR60_te_57, LTR60_te_64, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-888.png)<!-- -->

    ## Rows: 553 Columns: 81
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (80): LTR60B_te_39, LTR60B_te_40, LTR60B_te_50, LTR60B_te_83, LTR60B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-889.png)<!-- -->

    ## Rows: 553 Columns: 56
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (55): LTR61_te_39, LTR61_te_64, LTR61_te_89, LTR61_te_104, LTR61_te_110,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-890.png)<!-- -->

    ## Rows: 553 Columns: 77
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (76): LTR62_te_24, LTR62_te_26, LTR62_te_63, LTR62_te_70, LTR62_te_71, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-891.png)<!-- -->

    ## Rows: 553 Columns: 45
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (44): LTR64_te_14, LTR64_te_34, LTR64_te_54, LTR64_te_85, LTR64_te_88, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-892.png)<!-- -->

    ## Rows: 553 Columns: 124
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (123): LTR65_te_8, LTR65_te_22, LTR65_te_25, LTR65_te_26, LTR65_te_29, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-893.png)<!-- -->

    ## Rows: 553 Columns: 49
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (48): LTR66_te_43, LTR66_te_71, LTR66_te_72, LTR66_te_100, LTR66_te_110,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-894.png)<!-- -->

    ## Rows: 553 Columns: 82
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (81): LTR69_te_7, LTR69_te_25, LTR69_te_26, LTR69_te_32, LTR69_te_65, LT...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-895.png)<!-- -->

    ## Rows: 553 Columns: 74
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (73): LTR6A_te_10, LTR6A_te_24, LTR6A_te_25, LTR6A_te_53, LTR6A_te_56, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-896.png)<!-- -->

    ## Rows: 553 Columns: 33
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (32): LTR6B_te_16, LTR6B_te_52, LTR6B_te_222, LTR6B_te_226, LTR6B_te_233...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-897.png)<!-- -->

    ## Rows: 553 Columns: 95
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (94): LTR70_te_6, LTR70_te_14, LTR70_te_46, LTR70_te_50, LTR70_te_65, LT...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-898.png)<!-- -->

    ## Rows: 553 Columns: 50
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (49): LTR71A_te_51, LTR71A_te_69, LTR71A_te_70, LTR71A_te_75, LTR71A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-899.png)<!-- -->

    ## Rows: 553 Columns: 48
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (47): LTR71B_te_49, LTR71B_te_53, LTR71B_te_54, LTR71B_te_112, LTR71B_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-900.png)<!-- -->

    ## Rows: 553 Columns: 49
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (48): LTR72_te_33, LTR72_te_85, LTR72_te_100, LTR72_te_103, LTR72_te_111...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-901.png)<!-- -->

    ## Rows: 553 Columns: 46
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (45): LTR72B_te_4, LTR72B_te_22, LTR72B_te_23, LTR72B_te_31, LTR72B_te_3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-902.png)<!-- -->

    ## Rows: 553 Columns: 54
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (53): LTR73_te_52, LTR73_te_83, LTR73_te_84, LTR73_te_89, LTR73_te_90, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-903.png)<!-- -->

    ## Rows: 553 Columns: 88
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (87): LTR75_1_te_13, LTR75_1_te_14, LTR75_1_te_33, LTR75_1_te_41, LTR75_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-904.png)<!-- -->

    ## Rows: 553 Columns: 58
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (57): LTR76_te_22, LTR76_te_48, LTR76_te_58, LTR76_te_80, LTR76_te_133, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-905.png)<!-- -->

    ## Rows: 553 Columns: 64
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (63): LTR77_te_27, LTR77_te_44, LTR77_te_60, LTR77_te_61, LTR77_te_62, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-906.png)<!-- -->

    ## Rows: 553 Columns: 66
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (65): LTR77B_te_52, LTR77B_te_53, LTR77B_te_57, LTR77B_te_89, LTR77B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-907.png)<!-- -->

    ## Rows: 553 Columns: 28
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (27): LTR7A_te_32, LTR7A_te_52, LTR7A_te_57, LTR7A_te_101, LTR7A_te_114,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-908.png)<!-- -->

    ## Rows: 553 Columns: 17
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (16): LTR7B_te_52, LTR7B_te_55, LTR7B_te_85, LTR7B_te_86, LTR7B_te_142, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-909.png)<!-- -->

    ## Rows: 553 Columns: 50
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (49): LTR7C_te_17, LTR7C_te_18, LTR7C_te_32, LTR7C_te_57, LTR7C_te_58, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-910.png)<!-- -->

    ## Rows: 553 Columns: 20
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (19): LTR7Y_te_20, LTR7Y_te_25, LTR7Y_te_58, LTR7Y_te_59, LTR7Y_te_104, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-911.png)<!-- -->

    ## Rows: 553 Columns: 46
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (45): LTR8_te_7, LTR8_te_8, LTR8_te_47, LTR8_te_50, LTR8_te_55, LTR8_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-912.png)<!-- -->

    ## Rows: 553 Columns: 60
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (59): LTR8A_te_61, LTR8A_te_72, LTR8A_te_91, LTR8A_te_104, LTR8A_te_105,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-913.png)<!-- -->

    ## Rows: 553 Columns: 70
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (69): LTR8B_te_7, LTR8B_te_8, LTR8B_te_38, LTR8B_te_43, LTR8B_te_88, LTR...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-914.png)<!-- -->

    ## Rows: 553 Columns: 87
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (86): LTR9_te_6, LTR9_te_24, LTR9_te_40, LTR9_te_50, LTR9_te_57, LTR9_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-915.png)<!-- -->

    ## Rows: 553 Columns: 90
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (89): LTR9A1_te_12, LTR9A1_te_45, LTR9A1_te_57, LTR9A1_te_58, LTR9A1_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-916.png)<!-- -->

    ## Rows: 553 Columns: 65
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (64): LTR9B_te_60, LTR9B_te_68, LTR9B_te_69, LTR9B_te_95, LTR9B_te_96, L...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-917.png)<!-- -->

    ## Rows: 553 Columns: 69
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (68): LTR9C_te_24, LTR9C_te_59, LTR9C_te_60, LTR9C_te_128, LTR9C_te_135,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-918.png)<!-- -->

    ## Rows: 553 Columns: 76
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (75): LTR9D_te_6, LTR9D_te_12, LTR9D_te_35, LTR9D_te_49, LTR9D_te_50, LT...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-919.png)<!-- -->

    ## Rows: 553 Columns: 1260
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1259): MER101_I_te_456, MER101_I_te_457, MER101_I_te_461, MER101_I_te_4...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-920.png)<!-- -->

    ## Rows: 553 Columns: 36
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (35): MER101_te_30, MER101_te_48, MER101_te_49, MER101_te_87, MER101_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-921.png)<!-- -->

    ## Rows: 553 Columns: 81
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (80): MER101B_te_18, MER101B_te_21, MER101B_te_22, MER101B_te_42, MER101...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-922.png)<!-- -->

    ## Rows: 553 Columns: 8
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): familyname_position
    ## dbl (7): MER103_te_23, MER103_te_25, MER103_te_32, MER103_te_34, MER103_te_3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-923.png)<!-- -->

    ## Rows: 553 Columns: 20
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (19): MER105_te_17, MER105_te_25, MER105_te_37, MER105_te_39, MER105_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-924.png)<!-- -->

    ## Rows: 553 Columns: 28
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (27): MER106_te_35, MER106_te_36, MER106_te_49, MER106_te_57, MER106_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-925.png)<!-- -->

    ## Rows: 553 Columns: 5
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): familyname_position
    ## dbl (4): MER107_te_91, MER107_te_92, MER107_te_181, MER107_te_182
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-926.png)<!-- -->

    ## Rows: 553 Columns: 52
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (51): MER119_te_132, MER119_te_139, MER119_te_170, MER119_te_171, MER119...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-927.png)<!-- -->

    ## Rows: 553 Columns: 95
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (94): MER11A_te_6, MER11A_te_7, MER11A_te_23, MER11A_te_24, MER11A_te_27...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-928.png)<!-- -->

    ## Rows: 553 Columns: 112
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (111): MER11B_te_6, MER11B_te_7, MER11B_te_24, MER11B_te_27, MER11B_te_2...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-929.png)<!-- -->

    ## Rows: 553 Columns: 64
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (63): MER11C_te_27, MER11C_te_28, MER11C_te_37, MER11C_te_50, MER11C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-930.png)<!-- -->

    ## Rows: 553 Columns: 53
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (52): MER11D_te_30, MER11D_te_31, MER11D_te_110, MER11D_te_120, MER11D_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-931.png)<!-- -->

    ## Rows: 553 Columns: 10
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): familyname_position
    ## dbl (9): MER121_te_81, MER121_te_87, MER121_te_99, MER121_te_108, MER121_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-932.png)<!-- -->

    ## Rows: 553 Columns: 70
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (69): MER122_te_5, MER122_te_14, MER122_te_19, MER122_te_35, MER122_te_3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-933.png)<!-- -->

    ## Rows: 553 Columns: 81
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (80): MER1A_te_19, MER1A_te_25, MER1A_te_26, MER1A_te_27, MER1A_te_30, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-934.png)<!-- -->

    ## Rows: 553 Columns: 53
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (52): MER1B_te_18, MER1B_te_19, MER1B_te_23, MER1B_te_24, MER1B_te_25, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-935.png)<!-- -->

    ## Rows: 553 Columns: 37
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (36): MER2_te_12, MER2_te_13, MER2_te_19, MER2_te_20, MER2_te_42, MER2_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-936.png)<!-- -->

    ## Rows: 553 Columns: 6
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): familyname_position
    ## dbl (5): MER20_te_15, MER20_te_16, MER20_te_92, MER20_te_179, MER20_te_180
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-937.png)<!-- -->

    ## Rows: 553 Columns: 113
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (112): MER21_te_33, MER21_te_37, MER21_te_43, MER21_te_50, MER21_te_52, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-938.png)<!-- -->

    ## Rows: 553 Columns: 73
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (72): MER21A_te_37, MER21A_te_97, MER21A_te_98, MER21A_te_111, MER21A_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-939.png)<!-- -->

    ## Rows: 553 Columns: 68
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (67): MER21B_te_31, MER21B_te_32, MER21B_te_33, MER21B_te_37, MER21B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-940.png)<!-- -->

    ## Rows: 553 Columns: 74
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (73): MER21C_BT_te_15, MER21C_BT_te_30, MER21C_BT_te_39, MER21C_BT_te_65...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-941.png)<!-- -->

    ## Rows: 553 Columns: 64
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (63): MER21C_te_22, MER21C_te_35, MER21C_te_40, MER21C_te_94, MER21C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-942.png)<!-- -->

    ## Rows: 553 Columns: 551
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (550): MER21I_te_17, MER21I_te_18, MER21I_te_19, MER21I_te_20, MER21I_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-943.png)<!-- -->

    ## Rows: 553 Columns: 106
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (105): MER22_te_44, MER22_te_72, MER22_te_92, MER22_te_108, MER22_te_112...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-944.png)<!-- -->

    ## Rows: 553 Columns: 42
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (41): MER28_te_34, MER28_te_35, MER28_te_39, MER28_te_50, MER28_te_51, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-945.png)<!-- -->

    ## Rows: 553 Columns: 28
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (27): MER2B_te_19, MER2B_te_20, MER2B_te_45, MER2B_te_46, MER2B_te_47, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-946.png)<!-- -->

    ## Rows: 553 Columns: 11
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (10): MER3_te_4, MER3_te_28, MER3_te_51, MER3_te_52, MER3_te_63, MER3_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-947.png)<!-- -->

    ## Rows: 553 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (22): MER30_te_37, MER30_te_88, MER30_te_89, MER30_te_104, MER30_te_117,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-948.png)<!-- -->

    ## Rows: 553 Columns: 29
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (28): MER30B_te_11, MER30B_te_20, MER30B_te_21, MER30B_te_25, MER30B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-949.png)<!-- -->

    ## Rows: 553 Columns: 53
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (52): MER31_te_15, MER31_te_23, MER31_te_32, MER31_te_33, MER31_te_37, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-950.png)<!-- -->

    ## Rows: 553 Columns: 26
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (25): MER31A_te_14, MER31A_te_21, MER31A_te_40, MER31A_te_133, MER31A_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-951.png)<!-- -->

    ## Rows: 553 Columns: 48
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (47): MER31B_te_9, MER31B_te_11, MER31B_te_13, MER31B_te_42, MER31B_te_5...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-952.png)<!-- -->

    ## Rows: 553 Columns: 17
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (16): MER33_te_4, MER33_te_16, MER33_te_17, MER33_te_78, MER33_te_97, ME...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-953.png)<!-- -->

    ## Rows: 553 Columns: 66
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (65): MER34_te_52, MER34_te_79, MER34_te_89, MER34_te_99, MER34_te_109, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-954.png)<!-- -->

    ## Rows: 553 Columns: 61
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (60): MER34A_te_72, MER34A_te_88, MER34A_te_95, MER34A_te_98, MER34A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-955.png)<!-- -->

    ## Rows: 553 Columns: 50
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (49): MER34A1_te_142, MER34A1_te_143, MER34A1_te_154, MER34A1_te_180, ME...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-956.png)<!-- -->

    ## Rows: 553 Columns: 27
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (26): MER34B_te_36, MER34B_te_42, MER34B_te_52, MER34B_te_66, MER34B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-957.png)<!-- -->

    ## Rows: 553 Columns: 22
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (21): MER34C_te_80, MER34C_te_132, MER34C_te_133, MER34C_te_186, MER34C_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-958.png)<!-- -->

    ## Rows: 553 Columns: 33
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (32): MER34C2_te_70, MER34C2_te_72, MER34C2_te_107, MER34C2_te_146, MER3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-959.png)<!-- -->

    ## Rows: 553 Columns: 43
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (42): MER34D_te_7, MER34D_te_10, MER34D_te_14, MER34D_te_16, MER34D_te_1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-960.png)<!-- -->

    ## Rows: 553 Columns: 40
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (39): MER39_te_18, MER39_te_19, MER39_te_64, MER39_te_139, MER39_te_140,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-961.png)<!-- -->

    ## Rows: 553 Columns: 24
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (23): MER39B_te_43, MER39B_te_114, MER39B_te_115, MER39B_te_154, MER39B_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-962.png)<!-- -->

    ## Rows: 553 Columns: 47
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (46): MER41A_te_9, MER41A_te_11, MER41A_te_36, MER41A_te_56, MER41A_te_6...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-963.png)<!-- -->

    ## Rows: 553 Columns: 49
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (48): MER41B_te_10, MER41B_te_11, MER41B_te_56, MER41B_te_87, MER41B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-964.png)<!-- -->

    ## Rows: 553 Columns: 54
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (53): MER41C_te_11, MER41C_te_14, MER41C_te_15, MER41C_te_24, MER41C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-965.png)<!-- -->

    ## Rows: 553 Columns: 69
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (68): MER41D_te_11, MER41D_te_32, MER41D_te_94, MER41D_te_95, MER41D_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-966.png)<!-- -->

    ## Rows: 553 Columns: 61
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (60): MER41E_te_12, MER41E_te_13, MER41E_te_26, MER41E_te_27, MER41E_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-967.png)<!-- -->

    ## Rows: 553 Columns: 64
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (63): MER41F_te_27, MER41F_te_33, MER41F_te_38, MER41F_te_40, MER41F_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-968.png)<!-- -->

    ## Rows: 553 Columns: 97
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (96): MER41G_te_6, MER41G_te_23, MER41G_te_26, MER41G_te_27, MER41G_te_6...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-969.png)<!-- -->

    ## Rows: 553 Columns: 424
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (423): MER41I_te_28, MER41I_te_34, MER41I_te_42, MER41I_te_43, MER41I_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-970.png)<!-- -->

    ## Rows: 553 Columns: 26
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (25): MER44A_te_9, MER44A_te_19, MER44A_te_20, MER44A_te_21, MER44A_te_2...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-971.png)<!-- -->

    ## Rows: 553 Columns: 50
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (49): MER44B_te_95, MER44B_te_120, MER44B_te_121, MER44B_te_122, MER44B_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-972.png)<!-- -->

    ## Rows: 553 Columns: 24
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (23): MER45_te_30, MER45_te_35, MER45_te_36, MER45_te_37, MER45_te_45, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-973.png)<!-- -->

    ## Rows: 553 Columns: 96
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (95): MER45B_te_264, MER45B_te_266, MER45B_te_272, MER45B_te_274, MER45B...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-974.png)<!-- -->

    ## Rows: 553 Columns: 125
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (124): MER45R_te_513, MER45R_te_515, MER45R_te_516, MER45R_te_522, MER45...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-975.png)<!-- -->

    ## Rows: 553 Columns: 67
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (66): MER47B_te_35, MER47B_te_36, MER47B_te_47, MER47B_te_48, MER47B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-976.png)<!-- -->

    ## Rows: 553 Columns: 49
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (48): MER48_te_16, MER48_te_31, MER48_te_33, MER48_te_35, MER48_te_36, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-977.png)<!-- -->

    ## Rows: 553 Columns: 72
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (71): MER49_te_21, MER49_te_26, MER49_te_31, MER49_te_32, MER49_te_34, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-978.png)<!-- -->

    ## Rows: 553 Columns: 60
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (59): MER4A_te_21, MER4A_te_22, MER4A_te_30, MER4A_te_33, MER4A_te_70, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-979.png)<!-- -->

    ## Rows: 553 Columns: 34
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (33): MER4A1_LTR_te_21, MER4A1_LTR_te_22, MER4A1_LTR_te_30, MER4A1_LTR_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-980.png)<!-- -->

    ## Rows: 553 Columns: 35
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (34): MER4A1_te_30, MER4A1_te_107, MER4A1_te_127, MER4A1_te_154, MER4A1_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-981.png)<!-- -->

    ## Rows: 553 Columns: 60
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (59): MER4B_te_36, MER4B_te_37, MER4B_te_57, MER4B_te_58, MER4B_te_120, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-982.png)<!-- -->

    ## Rows: 553 Columns: 600
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (599): MER4BI_te_6, MER4BI_te_10, MER4BI_te_11, MER4BI_te_20, MER4BI_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-983.png)<!-- -->

    ## Rows: 553 Columns: 24
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (23): MER4C_te_22, MER4C_te_23, MER4C_te_39, MER4C_te_49, MER4C_te_84, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-984.png)<!-- -->

    ## Rows: 553 Columns: 73
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (72): MER4CL34_te_15, MER4CL34_te_21, MER4CL34_te_22, MER4CL34_te_76, ME...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-985.png)<!-- -->

    ## Rows: 553 Columns: 66
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (65): MER4D_LTR_te_35, MER4D_LTR_te_36, MER4D_LTR_te_37, MER4D_LTR_te_58...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-986.png)<!-- -->

    ## Rows: 553 Columns: 43
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (42): MER4D_te_29, MER4D_te_30, MER4D_te_35, MER4D_te_36, MER4D_te_48, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-987.png)<!-- -->

    ## Rows: 553 Columns: 32
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (31): MER4D1_te_37, MER4D1_te_57, MER4D1_te_58, MER4D1_te_78, MER4D1_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-988.png)<!-- -->

    ## Rows: 553 Columns: 35
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (34): MER4E_te_43, MER4E_te_55, MER4E_te_105, MER4E_te_111, MER4E_te_117...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-989.png)<!-- -->

    ## Rows: 553 Columns: 36
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (35): MER4E1_te_76, MER4E1_te_118, MER4E1_te_119, MER4E1_te_135, MER4E1_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-990.png)<!-- -->

    ## Rows: 553 Columns: 432
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (431): MER4I_te_31, MER4I_te_60, MER4I_te_61, MER4I_te_77, MER4I_te_99, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-991.png)<!-- -->

    ## Rows: 553 Columns: 76
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (75): MER50_te_6, MER50_te_20, MER50_te_24, MER50_te_58, MER50_te_75, ME...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-992.png)<!-- -->

    ## Rows: 553 Columns: 81
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (80): MER50B_te_24, MER50B_te_25, MER50B_te_43, MER50B_te_58, MER50B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-993.png)<!-- -->

    ## Rows: 553 Columns: 153
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (152): MER50C_te_6, MER50C_te_9, MER50C_te_15, MER50C_te_16, MER50C_te_2...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-994.png)<!-- -->

    ## Rows: 553 Columns: 945
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (944): MER50I_te_12, MER50I_te_13, MER50I_te_14, MER50I_te_21, MER50I_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-995.png)<!-- -->

    ## Rows: 553 Columns: 59
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (58): MER51A_te_42, MER51A_te_43, MER51A_te_44, MER51A_te_49, MER51A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-996.png)<!-- -->

    ## Rows: 553 Columns: 32
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (31): MER51B_te_42, MER51B_te_43, MER51B_te_61, MER51B_te_87, MER51B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-997.png)<!-- -->

    ## Rows: 553 Columns: 56
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (55): MER51C_te_6, MER51C_te_10, MER51C_te_33, MER51C_te_47, MER51C_te_7...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-998.png)<!-- -->

    ## Rows: 553 Columns: 96
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (95): MER51D_te_13, MER51D_te_16, MER51D_te_17, MER51D_te_25, MER51D_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-999.png)<!-- -->

    ## Rows: 553 Columns: 39
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (38): MER51E_te_45, MER51E_te_87, MER51E_te_95, MER51E_te_108, MER51E_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1000.png)<!-- -->

    ## Rows: 553 Columns: 714
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (713): MER51I_te_9, MER51I_te_10, MER51I_te_51, MER51I_te_62, MER51I_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1001.png)<!-- -->

    ## Rows: 553 Columns: 313
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (312): MER52A_te_19, MER52A_te_27, MER52A_te_32, MER52A_te_37, MER52A_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1002.png)<!-- -->

    ## Rows: 553 Columns: 830
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (829): MER52AI_te_5, MER52AI_te_17, MER52AI_te_18, MER52AI_te_26, MER52A...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1003.png)<!-- -->

    ## Rows: 553 Columns: 292
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (291): MER52B_te_10, MER52B_te_11, MER52B_te_14, MER52B_te_19, MER52B_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1004.png)<!-- -->

    ## Rows: 553 Columns: 152
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (151): MER52C_te_25, MER52C_te_28, MER52C_te_45, MER52C_te_46, MER52C_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1005.png)<!-- -->

    ## Rows: 553 Columns: 294
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (293): MER52D_te_31, MER52D_te_32, MER52D_te_37, MER52D_te_40, MER52D_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1006.png)<!-- -->

    ## Rows: 553 Columns: 154
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (153): MER54_te_32, MER54_te_33, MER54_te_53, MER54_te_56, MER54_te_83, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1007.png)<!-- -->

    ## Rows: 553 Columns: 113
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (112): MER54A_te_12, MER54A_te_15, MER54A_te_30, MER54A_te_32, MER54A_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1008.png)<!-- -->

    ## Rows: 553 Columns: 79
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (78): MER54B_te_30, MER54B_te_32, MER54B_te_33, MER54B_te_50, MER54B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1009.png)<!-- -->

    ## Rows: 553 Columns: 523
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (522): MER57A_I_te_6, MER57A_I_te_7, MER57A_I_te_22, MER57A_I_te_23, MER...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1010.png)<!-- -->

    ## Rows: 553 Columns: 39
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (38): MER57A1_te_10, MER57A1_te_37, MER57A1_te_38, MER57A1_te_65, MER57A...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1011.png)<!-- -->

    ## Rows: 553 Columns: 46
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (45): MER57B1_te_49, MER57B1_te_50, MER57B1_te_63, MER57B1_te_73, MER57B...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1012.png)<!-- -->

    ## Rows: 553 Columns: 48
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (47): MER57B2_te_49, MER57B2_te_50, MER57B2_te_66, MER57B2_te_67, MER57B...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1013.png)<!-- -->

    ## Rows: 553 Columns: 36
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (35): MER57C1_te_19, MER57C1_te_24, MER57C1_te_32, MER57C1_te_33, MER57C...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1014.png)<!-- -->

    ## Rows: 553 Columns: 24
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (23): MER57C2_te_24, MER57C2_te_36, MER57C2_te_37, MER57C2_te_38, MER57C...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1015.png)<!-- -->

    ## Rows: 553 Columns: 33
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (32): MER57D_te_45, MER57D_te_46, MER57D_te_73, MER57D_te_89, MER57D_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1016.png)<!-- -->

    ## Rows: 553 Columns: 59
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (58): MER57E1_te_24, MER57E1_te_31, MER57E1_te_37, MER57E1_te_38, MER57E...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1017.png)<!-- -->

    ## Rows: 553 Columns: 60
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (59): MER57E3_te_15, MER57E3_te_22, MER57E3_te_36, MER57E3_te_37, MER57E...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1018.png)<!-- -->

    ## Rows: 553 Columns: 36
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (35): MER57F_te_41, MER57F_te_55, MER57F_te_80, MER57F_te_107, MER57F_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1019.png)<!-- -->

    ## Rows: 553 Columns: 592
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (591): MER57I_te_6, MER57I_te_15, MER57I_te_23, MER57I_te_27, MER57I_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1020.png)<!-- -->

    ## Rows: 553 Columns: 48
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (47): MER58B_te_23, MER58B_te_24, MER58B_te_27, MER58B_te_36, MER58B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1021.png)<!-- -->

    ## Rows: 553 Columns: 36
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (35): MER58D_te_99, MER58D_te_103, MER58D_te_115, MER58D_te_121, MER58D_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1022.png)<!-- -->

    ## Rows: 553 Columns: 9
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): familyname_position
    ## dbl (8): MER5A_te_25, MER5A_te_26, MER5A_te_70, MER5A_te_76, MER5A_te_77, ME...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1023.png)<!-- -->

    ## Rows: 553 Columns: 5
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): familyname_position
    ## dbl (4): MER5A1_te_22, MER5A1_te_23, MER5A1_te_74, MER5A1_te_143
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1024.png)<!-- -->

    ## Rows: 553 Columns: 16
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (15): MER5B_te_5, MER5B_te_15, MER5B_te_16, MER5B_te_25, MER5B_te_30, ME...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1025.png)<!-- -->

    ## Rows: 553 Columns: 20
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (19): MER5C_te_68, MER5C_te_70, MER5C_te_71, MER5C_te_78, MER5C_te_164, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1026.png)<!-- -->

    ## Rows: 553 Columns: 67
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (66): MER6_te_200, MER6_te_203, MER6_te_206, MER6_te_208, MER6_te_212, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1027.png)<!-- -->

    ## Rows: 553 Columns: 46
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (45): MER61A_te_19, MER61A_te_28, MER61A_te_29, MER61A_te_50, MER61A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1028.png)<!-- -->

    ## Rows: 553 Columns: 57
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (56): MER61B_te_16, MER61B_te_32, MER61B_te_33, MER61B_te_35, MER61B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1029.png)<!-- -->

    ## Rows: 553 Columns: 47
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (46): MER61C_te_26, MER61C_te_35, MER61C_te_36, MER61C_te_46, MER61C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1030.png)<!-- -->

    ## Rows: 553 Columns: 58
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (57): MER61D_te_51, MER61D_te_52, MER61D_te_62, MER61D_te_63, MER61D_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1031.png)<!-- -->

    ## Rows: 553 Columns: 56
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (55): MER61E_te_24, MER61E_te_25, MER61E_te_53, MER61E_te_54, MER61E_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1032.png)<!-- -->

    ## Rows: 553 Columns: 68
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (67): MER61F_te_37, MER61F_te_42, MER61F_te_50, MER61F_te_56, MER61F_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1033.png)<!-- -->

    ## Rows: 553 Columns: 439
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (438): MER61I_te_16, MER61I_te_26, MER61I_te_31, MER61I_te_32, MER61I_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1034.png)<!-- -->

    ## Rows: 553 Columns: 39
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (38): MER65A_te_7, MER65A_te_8, MER65A_te_11, MER65A_te_25, MER65A_te_34...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1035.png)<!-- -->

    ## Rows: 553 Columns: 57
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (56): MER65C_te_14, MER65C_te_18, MER65C_te_19, MER65C_te_53, MER65C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1036.png)<!-- -->

    ## Rows: 553 Columns: 66
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (65): MER65D_te_44, MER65D_te_48, MER65D_te_54, MER65D_te_56, MER65D_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1037.png)<!-- -->

    ## Rows: 553 Columns: 968
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (967): MER66_I_te_100, MER66_I_te_110, MER66_I_te_178, MER66_I_te_201, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1038.png)<!-- -->

    ## Rows: 553 Columns: 68
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (67): MER66A_te_5, MER66A_te_19, MER66A_te_20, MER66A_te_21, MER66A_te_2...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1039.png)<!-- -->

    ## Rows: 553 Columns: 62
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (61): MER66B_te_20, MER66B_te_21, MER66B_te_26, MER66B_te_27, MER66B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1040.png)<!-- -->

    ## Rows: 553 Columns: 45
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (44): MER66C_te_14, MER66C_te_15, MER66C_te_77, MER66C_te_85, MER66C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1041.png)<!-- -->

    ## Rows: 553 Columns: 53
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (52): MER66D_te_31, MER66D_te_32, MER66D_te_34, MER66D_te_36, MER66D_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1042.png)<!-- -->

    ## Rows: 553 Columns: 59
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (58): MER67A_te_49, MER67A_te_51, MER67A_te_57, MER67A_te_58, MER67A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1043.png)<!-- -->

    ## Rows: 553 Columns: 68
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (67): MER67B_te_19, MER67B_te_20, MER67B_te_22, MER67B_te_24, MER67B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1044.png)<!-- -->

    ## Rows: 553 Columns: 37
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (36): MER67C_te_9, MER67C_te_25, MER67C_te_48, MER67C_te_68, MER67C_te_7...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1045.png)<!-- -->

    ## Rows: 553 Columns: 56
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (55): MER67D_te_19, MER67D_te_20, MER67D_te_51, MER67D_te_52, MER67D_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1046.png)<!-- -->

    ## Rows: 553 Columns: 42
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (41): MER68A_te_60, MER68A_te_100, MER68A_te_101, MER68A_te_123, MER68A_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1047.png)<!-- -->

    ## Rows: 553 Columns: 45
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (44): MER68B_te_40, MER68B_te_43, MER68B_te_44, MER68B_te_47, MER68B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1048.png)<!-- -->

    ## Rows: 553 Columns: 49
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (48): MER6A_te_156, MER6A_te_164, MER6A_te_170, MER6A_te_172, MER6A_te_1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1049.png)<!-- -->

    ## Rows: 553 Columns: 46
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (45): MER6B_te_21, MER6B_te_22, MER6B_te_23, MER6B_te_25, MER6B_te_26, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1050.png)<!-- -->

    ## Rows: 553 Columns: 40
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (39): MER6C_te_9, MER6C_te_14, MER6C_te_19, MER6C_te_20, MER6C_te_23, ME...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1051.png)<!-- -->

    ## Rows: 553 Columns: 72
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (71): MER70A_te_32, MER70A_te_43, MER70A_te_48, MER70A_te_54, MER70A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1052.png)<!-- -->

    ## Rows: 553 Columns: 81
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (80): MER70C_te_11, MER70C_te_12, MER70C_te_20, MER70C_te_28, MER70C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1053.png)<!-- -->

    ## Rows: 553 Columns: 68
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (67): MER72_te_20, MER72_te_44, MER72_te_49, MER72_te_50, MER72_te_83, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1054.png)<!-- -->

    ## Rows: 553 Columns: 83
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (82): MER72B_te_20, MER72B_te_28, MER72B_te_33, MER72B_te_35, MER72B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1055.png)<!-- -->

    ## Rows: 553 Columns: 68
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (67): MER73_te_33, MER73_te_34, MER73_te_63, MER73_te_77, MER73_te_78, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1056.png)<!-- -->

    ## Rows: 553 Columns: 86
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (85): MER74_te_71, MER74_te_77, MER74_te_78, MER74_te_79, MER74_te_92, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1057.png)<!-- -->

    ## Rows: 553 Columns: 55
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (54): MER74A_te_9, MER74A_te_10, MER74A_te_15, MER74A_te_33, MER74A_te_7...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1058.png)<!-- -->

    ## Rows: 553 Columns: 81
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (80): MER74C_te_51, MER74C_te_52, MER74C_te_76, MER74C_te_77, MER74C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1059.png)<!-- -->

    ## Rows: 553 Columns: 61
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (60): MER75_te_10, MER75_te_11, MER75_te_19, MER75_te_20, MER75_te_29, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1060.png)<!-- -->

    ## Rows: 553 Columns: 11
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (10): MER75B_te_10, MER75B_te_11, MER75B_te_19, MER75B_te_20, MER75B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1061.png)<!-- -->

    ## Rows: 553 Columns: 56
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (55): MER76_te_4, MER76_te_20, MER76_te_31, MER76_te_32, MER76_te_50, ME...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1062.png)<!-- -->

    ## Rows: 553 Columns: 77
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (76): MER77_te_34, MER77_te_37, MER77_te_64, MER77_te_65, MER77_te_77, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1063.png)<!-- -->

    ## Rows: 553 Columns: 78
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (77): MER82_te_46, MER82_te_47, MER82_te_48, MER82_te_49, MER82_te_51, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1064.png)<!-- -->

    ## Rows: 553 Columns: 49
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (48): MER83_te_37, MER83_te_38, MER83_te_42, MER83_te_45, MER83_te_57, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1065.png)<!-- -->

    ## Rows: 553 Columns: 35
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (34): MER83B_te_8, MER83B_te_29, MER83B_te_38, MER83B_te_47, MER83B_te_5...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1066.png)<!-- -->

    ## Rows: 553 Columns: 687
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (686): MER83BI_te_72, MER83BI_te_73, MER83BI_te_74, MER83BI_te_90, MER83...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1067.png)<!-- -->

    ## Rows: 553 Columns: 39
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (38): MER83C_te_12, MER83C_te_33, MER83C_te_34, MER83C_te_44, MER83C_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1068.png)<!-- -->

    ## Rows: 553 Columns: 65
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (64): MER84_te_16, MER84_te_34, MER84_te_38, MER84_te_44, MER84_te_45, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1069.png)<!-- -->

    ## Rows: 553 Columns: 737
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (736): MER84I_te_9, MER84I_te_13, MER84I_te_16, MER84I_te_25, MER84I_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1070.png)<!-- -->

    ## Rows: 553 Columns: 47
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (46): MER87B_te_76, MER87B_te_77, MER87B_te_110, MER87B_te_138, MER87B_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1071.png)<!-- -->

    ## Rows: 553 Columns: 81
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (80): MER88_te_24, MER88_te_25, MER88_te_32, MER88_te_59, MER88_te_67, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1072.png)<!-- -->

    ## Rows: 553 Columns: 37
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (36): MER89_te_21, MER89_te_22, MER89_te_44, MER89_te_46, MER89_te_81, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1073.png)<!-- -->

    ## Rows: 553 Columns: 49
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (48): MER9_te_5, MER9_te_79, MER9_te_88, MER9_te_90, MER9_te_97, MER9_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1074.png)<!-- -->

    ## Rows: 553 Columns: 86
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (85): MER90a_LTR_te_16, MER90a_LTR_te_17, MER90a_LTR_te_36, MER90a_LTR_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1075.png)<!-- -->

    ## Rows: 553 Columns: 67
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (66): MER92B_te_14, MER92B_te_32, MER92B_te_62, MER92B_te_68, MER92B_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1076.png)<!-- -->

    ## Rows: 553 Columns: 56
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (55): MER95_te_41, MER95_te_43, MER95_te_54, MER95_te_55, MER95_te_58, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1077.png)<!-- -->

    ## Rows: 553 Columns: 24
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (23): MER96B_te_121, MER96B_te_128, MER96B_te_129, MER96B_te_130, MER96B...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1078.png)<!-- -->

    ## Rows: 553 Columns: 45
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (44): MER9B_te_20, MER9B_te_87, MER9B_te_96, MER9B_te_98, MER9B_te_100, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1079.png)<!-- -->

    ## Rows: 553 Columns: 12
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (11): MIR_te_17, MIR_te_21, MIR_te_27, MIR_te_55, MIR_te_56, MIR_te_62, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1080.png)<!-- -->

    ## Rows: 553 Columns: 9
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): familyname_position
    ## dbl (8): MIR3_te_70, MIR3_te_71, MIR3_te_159, MIR3_te_186, MIR3_te_187, MIR3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1081.png)<!-- -->

    ## Rows: 553 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (22): MIRb_te_48, MIRb_te_49, MIRb_te_53, MIRb_te_62, MIRb_te_63, MIRb_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1082.png)<!-- -->

    ## Rows: 553 Columns: 145
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (144): MLT-int_te_1, MLT-int_te_12, MLT-int_te_13, MLT-int_te_57, MLT-in...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1083.png)<!-- -->

    ## Rows: 553 Columns: 133
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (132): MLT1_I_te_12, MLT1_I_te_13, MLT1_I_te_68, MLT1_I_te_84, MLT1_I_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1084.png)<!-- -->

    ## Rows: 553 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (22): MLT1A0_te_47, MLT1A0_te_48, MLT1A0_te_54, MLT1A0_te_73, MLT1A0_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1085.png)<!-- -->

    ## Rows: 553 Columns: 26
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (25): MLT1A1_te_20, MLT1A1_te_55, MLT1A1_te_56, MLT1A1_te_63, MLT1A1_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1086.png)<!-- -->

    ## Rows: 553 Columns: 13
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (12): MLT1B_te_91, MLT1B_te_213, MLT1B_te_229, MLT1B_te_230, MLT1B_te_23...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1087.png)<!-- -->

    ## Rows: 553 Columns: 21
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (20): MLT1C_te_92, MLT1C_te_93, MLT1C_te_197, MLT1C_te_198, MLT1C_te_212...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1088.png)<!-- -->

    ## Rows: 553 Columns: 24
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (23): MLT1C1_te_26, MLT1C1_te_27, MLT1C1_te_69, MLT1C1_te_127, MLT1C1_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1089.png)<!-- -->

    ## Rows: 553 Columns: 18
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (17): MLT1C2_te_39, MLT1C2_te_53, MLT1C2_te_54, MLT1C2_te_220, MLT1C2_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1090.png)<!-- -->

    ## Rows: 553 Columns: 31
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (30): MLT1D_te_11, MLT1D_te_38, MLT1D_te_39, MLT1D_te_50, MLT1D_te_51, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1091.png)<!-- -->

    ## Rows: 553 Columns: 55
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (54): MLT1E_te_19, MLT1E_te_29, MLT1E_te_38, MLT1E_te_39, MLT1E_te_41, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1092.png)<!-- -->

    ## Rows: 553 Columns: 77
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (76): MLT1E1_te_26, MLT1E1_te_39, MLT1E1_te_55, MLT1E1_te_56, MLT1E1_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1093.png)<!-- -->

    ## Rows: 553 Columns: 46
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (45): MLT1E1A_te_61, MLT1E1A_te_64, MLT1E1A_te_76, MLT1E1A_te_81, MLT1E1...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1094.png)<!-- -->

    ## Rows: 553 Columns: 49
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (48): MLT1E2_te_37, MLT1E2_te_67, MLT1E2_te_68, MLT1E2_te_86, MLT1E2_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1095.png)<!-- -->

    ## Rows: 553 Columns: 186
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (185): MLT1F_I_te_73, MLT1F_I_te_79, MLT1F_I_te_102, MLT1F_I_te_104, MLT...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1096.png)<!-- -->

    ## Rows: 553 Columns: 51
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (50): MLT1F_te_27, MLT1F_te_31, MLT1F_te_34, MLT1F_te_40, MLT1F_te_50, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1097.png)<!-- -->

    ## Rows: 553 Columns: 41
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (40): MLT1F1_te_39, MLT1F1_te_44, MLT1F1_te_49, MLT1F1_te_52, MLT1F1_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1098.png)<!-- -->

    ## Rows: 553 Columns: 36
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (35): MLT1F2_te_14, MLT1F2_te_125, MLT1F2_te_175, MLT1F2_te_180, MLT1F2_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1099.png)<!-- -->

    ## Rows: 553 Columns: 83
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (82): MLT1G1_te_10, MLT1G1_te_30, MLT1G1_te_31, MLT1G1_te_35, MLT1G1_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1100.png)<!-- -->

    ## Rows: 553 Columns: 110
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (109): MLT1G2_te_90, MLT1G2_te_100, MLT1G2_te_102, MLT1G2_te_104, MLT1G2...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1101.png)<!-- -->

    ## Rows: 553 Columns: 28
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (27): MLT1G3_te_192, MLT1G3_te_207, MLT1G3_te_238, MLT1G3_te_239, MLT1G3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1102.png)<!-- -->

    ## Rows: 553 Columns: 57
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (56): MLT1H_te_28, MLT1H_te_30, MLT1H_te_33, MLT1H_te_56, MLT1H_te_98, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1103.png)<!-- -->

    ## Rows: 553 Columns: 60
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (59): MLT1H1_te_29, MLT1H1_te_34, MLT1H1_te_35, MLT1H1_te_38, MLT1H1_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1104.png)<!-- -->

    ## Rows: 553 Columns: 41
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (40): MLT1H2_te_27, MLT1H2_te_31, MLT1H2_te_34, MLT1H2_te_40, MLT1H2_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1105.png)<!-- -->

    ## Rows: 553 Columns: 8
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): familyname_position
    ## dbl (7): MLT1I_te_143, MLT1I_te_144, MLT1I_te_156, MLT1I_te_159, MLT1I_te_16...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1106.png)<!-- -->

    ## Rows: 553 Columns: 19
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (18): MLT1J2_te_50, MLT1J2_te_101, MLT1J2_te_146, MLT1J2_te_148, MLT1J2_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1107.png)<!-- -->

    ## Rows: 553 Columns: 28
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (27): MLT2A1_te_88, MLT2A1_te_128, MLT2A1_te_152, MLT2A1_te_170, MLT2A1_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1108.png)<!-- -->

    ## Rows: 553 Columns: 41
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (40): MLT2A2_te_15, MLT2A2_te_56, MLT2A2_te_112, MLT2A2_te_119, MLT2A2_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1109.png)<!-- -->

    ## Rows: 553 Columns: 44
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (43): MLT2B2_te_60, MLT2B2_te_97, MLT2B2_te_98, MLT2B2_te_101, MLT2B2_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1110.png)<!-- -->

    ## Rows: 553 Columns: 61
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (60): MLT2B3_te_11, MLT2B3_te_17, MLT2B3_te_81, MLT2B3_te_101, MLT2B3_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1111.png)<!-- -->

    ## Rows: 553 Columns: 32
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (31): MLT2B4_te_39, MLT2B4_te_40, MLT2B4_te_60, MLT2B4_te_79, MLT2B4_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1112.png)<!-- -->

    ## Rows: 553 Columns: 40
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (39): MLT2C2_te_7, MLT2C2_te_37, MLT2C2_te_43, MLT2C2_te_44, MLT2C2_te_4...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1113.png)<!-- -->

    ## Rows: 553 Columns: 37
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (36): MLT2D_te_18, MLT2D_te_20, MLT2D_te_25, MLT2D_te_33, MLT2D_te_37, M...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1114.png)<!-- -->

    ## Rows: 553 Columns: 22
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (21): MSR1_te_47, MSR1_te_54, MSR1_te_60, MSR1_te_65, MSR1_te_67, MSR1_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1115.png)<!-- -->

    ## Rows: 553 Columns: 160
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (159): MST_I_te_23, MST_I_te_57, MST_I_te_58, MST_I_te_59, MST_I_te_76, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1116.png)<!-- -->

    ## Rows: 553 Columns: 31
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (30): MSTA_te_14, MSTA_te_41, MSTA_te_57, MSTA_te_67, MSTA_te_139, MSTA_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1117.png)<!-- -->

    ## Rows: 553 Columns: 24
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (23): MSTA1_te_58, MSTA1_te_143, MSTA1_te_145, MSTA1_te_146, MSTA1_te_14...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1118.png)<!-- -->

    ## Rows: 553 Columns: 25
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (24): MSTA2_te_63, MSTA2_te_69, MSTA2_te_70, MSTA2_te_145, MSTA2_te_155,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1119.png)<!-- -->

    ## Rows: 553 Columns: 46
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (45): MSTB_te_48, MSTB_te_51, MSTB_te_56, MSTB_te_90, MSTB_te_100, MSTB_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1120.png)<!-- -->

    ## Rows: 553 Columns: 26
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (25): MSTB1_te_32, MSTB1_te_56, MSTB1_te_100, MSTB1_te_101, MSTB1_te_133...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1121.png)<!-- -->

    ## Rows: 553 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (22): MSTC_te_16, MSTC_te_45, MSTC_te_98, MSTC_te_150, MSTC_te_182, MSTC...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1122.png)<!-- -->

    ## Rows: 553 Columns: 42
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (41): MSTD_te_69, MSTD_te_78, MSTD_te_83, MSTD_te_85, MSTD_te_133, MSTD_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1123.png)<!-- -->

    ## Rows: 553 Columns: 19
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (18): ORSL_te_13, ORSL_te_14, ORSL_te_22, ORSL_te_31, ORSL_te_50, ORSL_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1124.png)<!-- -->

    ## Rows: 553 Columns: 32
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (31): PABL_A_te_24, PABL_A_te_25, PABL_A_te_51, PABL_A_te_88, PABL_A_te_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1125.png)<!-- -->

    ## Rows: 553 Columns: 332
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (331): PABL_AI_te_11, PABL_AI_te_13, PABL_AI_te_21, PABL_AI_te_28, PABL_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1126.png)<!-- -->

    ## Rows: 553 Columns: 44
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (43): PABL_B_te_47, PABL_B_te_98, PABL_B_te_109, PABL_B_te_158, PABL_B_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1127.png)<!-- -->

    ## Rows: 553 Columns: 1139
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr    (1): familyname_position
    ## dbl (1138): PRIMA4_I_te_12, PRIMA4_I_te_13, PRIMA4_I_te_14, PRIMA4_I_te_31, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1128.png)<!-- -->

    ## Rows: 553 Columns: 34
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (33): PRIMA4_LTR_te_19, PRIMA4_LTR_te_39, PRIMA4_LTR_te_75, PRIMA4_LTR_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1129.png)<!-- -->

    ## Rows: 553 Columns: 774
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (773): PRIMA41_te_10, PRIMA41_te_11, PRIMA41_te_17, PRIMA41_te_18, PRIMA...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1130.png)<!-- -->

    ## Rows: 553 Columns: 60
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (59): PrimLTR79_te_14, PrimLTR79_te_24, PrimLTR79_te_32, PrimLTR79_te_33...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1131.png)<!-- -->

    ## Rows: 553 Columns: 93
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (92): PTR5_te_14, PTR5_te_43, PTR5_te_48, PTR5_te_51, PTR5_te_52, PTR5_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1132.png)<!-- -->

    ## Rows: 553 Columns: 241
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (240): RICKSHA_0_te_30, RICKSHA_0_te_34, RICKSHA_0_te_35, RICKSHA_0_te_3...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1133.png)<!-- -->

    ## Rows: 553 Columns: 194
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (193): Ricksha_a_te_17, Ricksha_a_te_19, Ricksha_a_te_20, Ricksha_a_te_2...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1134.png)<!-- -->

    ## Rows: 553 Columns: 232
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (231): RICKSHA_te_18, RICKSHA_te_19, RICKSHA_te_20, RICKSHA_te_25, RICKS...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1135.png)<!-- -->

    ## Rows: 553 Columns: 42
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (41): SATR1_te_21, SATR1_te_30, SATR1_te_35, SATR1_te_77, SATR1_te_106, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1136.png)<!-- -->

    ## Rows: 553 Columns: 47
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (46): SATR2_te_22, SATR2_te_25, SATR2_te_33, SATR2_te_35, SATR2_te_41, S...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1137.png)<!-- -->

    ## Rows: 553 Columns: 24
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (23): SN5_te_6, SN5_te_36, SN5_te_41, SN5_te_44, SN5_te_51, SN5_te_140, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1138.png)<!-- -->

    ## Rows: 553 Columns: 46
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (45): SVA_A_te_92, SVA_A_te_155, SVA_A_te_166, SVA_A_te_251, SVA_A_te_26...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1139.png)<!-- -->

    ## Rows: 553 Columns: 24
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (23): SVA2_te_2, SVA2_te_13, SVA2_te_32, SVA2_te_36, SVA2_te_40, SVA2_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1140.png)<!-- -->

    ## Rows: 553 Columns: 177
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (176): TAR1_te_1, TAR1_te_27, TAR1_te_44, TAR1_te_49, TAR1_te_53, TAR1_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1141.png)<!-- -->

    ## Rows: 553 Columns: 125
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (124): THE1_I_te_26, THE1_I_te_27, THE1_I_te_28, THE1_I_te_34, THE1_I_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1142.png)<!-- -->

    ## Rows: 553 Columns: 27
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (26): THE1A_te_64, THE1A_te_65, THE1A_te_110, THE1A_te_111, THE1A_te_120...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1143.png)<!-- -->

    ## Rows: 553 Columns: 34
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (33): THE1B_te_64, THE1B_te_65, THE1B_te_69, THE1B_te_70, THE1B_te_82, T...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1144.png)<!-- -->

    ## Rows: 553 Columns: 33
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (32): THE1C_te_64, THE1C_te_65, THE1C_te_69, THE1C_te_70, THE1C_te_82, T...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1145.png)<!-- -->

    ## Rows: 553 Columns: 27
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (26): THE1D_te_40, THE1D_te_65, THE1D_te_66, THE1D_te_71, THE1D_te_78, T...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1146.png)<!-- -->

    ## Rows: 553 Columns: 15
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (14): THER1_te_75, THER1_te_84, THER1_te_90, THER1_te_96, THER1_te_98, T...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1147.png)<!-- -->

    ## Rows: 553 Columns: 129
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (128): TIGGER1_te_12, TIGGER1_te_13, TIGGER1_te_22, TIGGER1_te_23, TIGGE...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1148.png)<!-- -->

    ## Rows: 553 Columns: 164
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (163): TIGGER2_te_380, TIGGER2_te_393, TIGGER2_te_399, TIGGER2_te_411, T...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1149.png)<!-- -->

    ## Rows: 553 Columns: 65
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (64): Tigger2b_Pri_te_391, Tigger2b_Pri_te_405, Tigger2b_Pri_te_407, Tig...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## `summarise()` has grouped output by 'location'. You can override using the `.groups` argument.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1150.png)<!-- -->

    ## Rows: 553 Columns: 93
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (92): Tigger3b_te_22, Tigger3b_te_23, Tigger3b_te_27, Tigger3b_te_28, Ti...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1151.png)<!-- -->

    ## Rows: 553 Columns: 43
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (42): Tigger3d_te_23, Tigger3d_te_27, Tigger3d_te_30, Tigger3d_te_34, Ti...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1152.png)<!-- -->

    ## Rows: 553 Columns: 29
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (28): Tigger4a_te_21, Tigger4a_te_22, Tigger4a_te_29, Tigger4a_te_30, Ti...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1153.png)<!-- -->

    ## Rows: 553 Columns: 34
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (33): TIGGER5_A_te_116, TIGGER5_A_te_117, TIGGER5_A_te_125, TIGGER5_A_te...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1154.png)<!-- -->

    ## Rows: 553 Columns: 36
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (35): TIGGER5A_te_5, TIGGER5A_te_6, TIGGER5A_te_12, TIGGER5A_te_13, TIGG...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1155.png)<!-- -->

    ## Rows: 553 Columns: 199
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (198): TIGGER6A_te_22, TIGGER6A_te_27, TIGGER6A_te_28, TIGGER6A_te_52, T...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1156.png)<!-- -->

    ## Rows: 553 Columns: 295
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (294): TIGGER6B_te_30, TIGGER6B_te_37, TIGGER6B_te_50, TIGGER6B_te_53, T...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1157.png)<!-- -->

    ## Rows: 553 Columns: 243
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (242): TIGGER7_te_291, TIGGER7_te_292, TIGGER7_te_295, TIGGER7_te_296, T...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1158.png)<!-- -->

    ## Rows: 553 Columns: 110
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (109): Tigger9b_te_12, Tigger9b_te_16, Tigger9b_te_19, Tigger9b_te_20, T...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1159.png)<!-- -->

    ## Rows: 553 Columns: 25
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (24): ZOMBI_A_te_21, ZOMBI_A_te_22, ZOMBI_A_te_29, ZOMBI_A_te_30, ZOMBI_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1160.png)<!-- -->

    ## Rows: 553 Columns: 16
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): familyname_position
    ## dbl (15): ZOMBI_B_te_29, ZOMBI_B_te_30, ZOMBI_B_te_67, ZOMBI_B_te_116, ZOMBI...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1161.png)<!-- -->

    ## Rows: 553 Columns: 166
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr   (1): familyname_position
    ## dbl (165): ZOMBI_te_29, ZOMBI_te_30, ZOMBI_te_54, ZOMBI_te_67, ZOMBI_te_68, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](private-snps_files/figure-gfm/unnamed-chunk-5-1162.png)<!-- -->![](private-snps_files/figure-gfm/unnamed-chunk-5-1163.png)<!-- -->
