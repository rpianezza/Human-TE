HGDP - Controls
================

Script 5. This scripts contains some controls on the data quality, such
as checking for **batch effects** and **reads length**.

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

The two following functions have been previously described in scripts 1
and 4. The first plot the distribution of the **copynumbers** of a TE
divided by **continent**. The second plot the **mean copynumber** of
every **population** for a particular TE on the **world map**.

These functions are slightly modified in order to add `batch` among the
arguments and check for batch effect.

``` r
plotTEfamily <- function(data, famname, bat, binwidht, x_title, y_title, x_numbers, y_numbers){
filtered <- filter(data, familyname==famname, batch==bat)
ggplot(data = filtered, mapping = aes(x = copynumber, fill = Country)) +
  geom_histogram(binwidth = binwidht) +
  facet_wrap(~sex)+
  ggtitle(paste0(filtered$batch, ' - ',famname)) + theme(plot.title = element_text(size = 8, hjust = 0.5)) +
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

``` r
plot_map <- function(data, famname, batch){
TE <- filter(data, familyname == famname)
world_map = map_data("world")

ggplot() +
  geom_map(
    data = world_map, map = world_map,
    aes(long, lat, map_id = region),
    color = "white", fill = "lightgray", size = 0) +
  geom_point(
    data = TE, aes(longitude, latitude, color = copynumber, size = count)
  ) + scale_colour_gradient(low = "green", high = "red") + theme(legend.position="top") + theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~sex) + ggtitle(paste0(batch, ' - ',famname))}
```

## Batch effect

In our dataset, the `batch` column refers to the computer on which the
pipeline was ran for each sample, thus `flo` for Florian’s PC and `ro`
for Robert’s PC. This could introduce some biases, even though the
pipeline are the same, with the exact same version for each of the used
programs. To check if this is the case and if the bias is significant, I
did the same plots previously described in scripts 1-2-4 but dividing
the dataset for the two **batches**. I choose the TE `L1ME5` which was
the most interesting.

``` r
(ro_L1ME5_plot<-plotTEfamily(HGDPcutoff, 'L1ME5', 'ro', 1, 'y', 'y', 'y', 'y'))
```

![](5_HGDP_controls_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
(flo_L1ME5_plot<-plotTEfamily(HGDPcutoff, 'L1ME5', 'flo', 1, 'y', 'y', 'y', 'y'))
```

![](5_HGDP_controls_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
coordinates <- read_tsv("/Users/rpianezza/TE/summary-HGDP/HGDP_populationcoordinates.txt", col_names = c("Pop", "region", "latitude", "longitude"))
```

    ## Rows: 54 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): Pop, region
    ## dbl (2): latitude, longitude
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
coord <- select(coordinates, Pop, latitude, longitude)

ro_TE <- filter(HGDPcutoff, type=='te', batch=='ro')
flo_TE <- filter(HGDPcutoff, type=='te', batch=='flo')

ro_by_pop <- group_by(ro_TE, Pop, familyname, sex) %>% summarise(copynumber = mean(copynumber), count=n())
```

    ## `summarise()` has grouped output by 'Pop', 'familyname'. You can override using
    ## the `.groups` argument.

``` r
ro_data <- inner_join(x = coord, y = ro_by_pop, by = "Pop")

flo_by_pop <- group_by(flo_TE, Pop, familyname, sex) %>% summarise(copynumber = mean(copynumber), count=n())
```

    ## `summarise()` has grouped output by 'Pop', 'familyname'. You can override using
    ## the `.groups` argument.

``` r
flo_data <- inner_join(x = coord, y = flo_by_pop, by = "Pop")

(ro_L1ME5_map<-plot_map(ro_data, 'L1ME5', 'ro'))
```

    ## Warning: Ignoring unknown aesthetics: x, y

![](5_HGDP_controls_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
(flo_L1ME5_map<-plot_map(flo_data, 'L1ME5', 'flo'))
```

    ## Warning: Ignoring unknown aesthetics: x, y

![](5_HGDP_controls_files/figure-gfm/unnamed-chunk-5-2.png)<!-- --> I do
not notice significant differences among the two batches, in neither of
the two figures.

## Read length

Another possible bias could be introduced by different **read lengths**
of the samples. To check this, I found the original dataset here
<https://www.ebi.ac.uk/ena/browser/view/PRJEB6463> on the ENA website.
Then I selected the columns containing the number of reads for each
sample (`read_count`) and the total base number sequenced in each sample
(`base_count`), and finally I downloaded the report in TSV format. Here
I import the report in R as `HGDP_report` and I add the column
`read_lenght`, calculated as the ratio between `base_count` and
`read_count`.

``` r
(HGDP_report <- read_tsv('/Users/rpianezza/TE/ric-documentation-Rmd/other-files/filereport_read_run_PRJEB6463_tsv.txt') %>% mutate(read_length = base_count/read_count))
```

    ## Rows: 7942 Columns: 7
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (5): study_accession, sample_accession, experiment_accession, run_access...
    ## dbl (2): read_count, base_count
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

    ## # A tibble: 7,942 × 8
    ##    study_accession sample_acce…¹ exper…² run_a…³ read_…⁴ base_…⁵ sampl…⁶ read_…⁷
    ##    <chr>           <chr>         <chr>   <chr>     <dbl>   <dbl> <chr>     <dbl>
    ##  1 PRJEB6463       SAMEA2580809  ERX132… ERR125…  8.40e7 1.27e10 HGDP00…     151
    ##  2 PRJEB6463       SAMEA2580910  ERX132… ERR125…  6.83e7 1.03e10 HGDP00…     151
    ##  3 PRJEB6463       SAMEA2580921  ERX132… ERR125…  7.51e7 1.13e10 HGDP00…     151
    ##  4 PRJEB6463       SAMEA2580929  ERX132… ERR125…  7.36e7 1.11e10 HGDP00…     151
    ##  5 PRJEB6463       SAMEA2580937  ERX132… ERR125…  7.38e7 1.11e10 HGDP00…     151
    ##  6 PRJEB6463       SAMEA2580945  ERX132… ERR125…  8.66e7 1.31e10 HGDP01…     151
    ##  7 PRJEB6463       SAMEA2580957  ERX132… ERR125…  7.08e7 1.07e10 HGDP01…     151
    ##  8 PRJEB6463       SAMEA2580969  ERX132… ERR125…  7.14e7 1.08e10 HGDP01…     151
    ##  9 PRJEB6463       SAMEA2580981  ERX132… ERR125…  6.67e7 1.01e10 HGDP01…     151
    ## 10 PRJEB6463       SAMEA2580992  ERX132… ERR125…  7.36e7 1.11e10 HGDP01…     151
    ## # … with 7,932 more rows, and abbreviated variable names ¹​sample_accession,
    ## #   ²​experiment_accession, ³​run_accession, ⁴​read_count, ⁵​base_count,
    ## #   ⁶​sample_alias, ⁷​read_length

First I wanted to check if the dataset used for my data analysis was
still complete respect to the original one, so I compared the number of
individuals in the two datasets. It’s `828` in both, confirming the
completness of the used dataset.

``` r
length(unique(HGDP_report$sample_accession)) == length(unique(HGDPcutoff$ID))
```

    ## [1] TRUE

Then I checked for the mean of read lenghts among all the samples.
Surprisingly it is not `151`, as I was expecting looking at the dataset,
but a bit lower.

``` r
(mean_rl = mean(HGDP_report$read_length))
```

    ## [1] 150.9875

Thus, I checked for the samples with `read_lenght` different from the
expected `151`.

``` r
(not151 <- filter(HGDP_report, !(read_length==151)))
```

    ## # A tibble: 14 × 8
    ##    study_accession sample_acce…¹ exper…² run_a…³ read_…⁴ base_…⁵ sampl…⁶ read_…⁷
    ##    <chr>           <chr>         <chr>   <chr>     <dbl>   <dbl> <chr>     <dbl>
    ##  1 PRJEB6463       SAMEA2580910  ERX136… ERR128…  7.15e7 1.03e10 HGDP00…    144.
    ##  2 PRJEB6463       SAMEA2580945  ERX136… ERR128…  8.97e7 1.30e10 HGDP01…    144.
    ##  3 PRJEB6463       SAMEA2581005  ERX136… ERR128…  7.81e7 1.13e10 HGDP01…    144.
    ##  4 PRJEB6463       SAMEA2580922  ERX136… ERR128…  7.10e7 1.02e10 HGDP00…    144.
    ##  5 PRJEB6463       SAMEA2580930  ERX136… ERR128…  1.38e8 1.95e10 HGDP00…    141.
    ##  6 PRJEB6463       SAMEA2580960  ERX136… ERR128…  9.43e7 1.35e10 HGDP01…    143.
    ##  7 PRJEB6463       SAMEA2580972  ERX136… ERR128…  6.93e7 9.97e 9 HGDP01…    144.
    ##  8 PRJEB6463       SAMEA2580948  ERX136… ERR128…  7.94e7 1.14e10 HGDP01…    144.
    ##  9 PRJEB6463       SAMEA2580985  ERX136… ERR128…  8.24e7 1.19e10 HGDP01…    144.
    ## 10 PRJEB6463       SAMEA2580965  ERX136… ERR128…  8.05e7 1.15e10 HGDP01…    143.
    ## 11 PRJEB6463       SAMEA2580927  ERX136… ERR128…  8.08e7 1.17e10 HGDP00…    144.
    ## 12 PRJEB6463       SAMEA2581017  ERX136… ERR128…  8.17e7 1.18e10 HGDP01…    144.
    ## 13 PRJEB6463       SAMEA2581004  ERX136… ERR128…  7.50e7 1.09e10 HGDP01…    145.
    ## 14 PRJEB6463       SAMEA2581018  ERX136… ERR128…  7.25e7 1.05e10 HGDP01…    145.
    ## # … with abbreviated variable names ¹​sample_accession, ²​experiment_accession,
    ## #   ³​run_accession, ⁴​read_count, ⁵​base_count, ⁶​sample_alias, ⁷​read_length

We found 14/828 samples for which read length is not `151`, but in a
range between `141` and `146`. I think that these outliers are not
introducing a significant bias, since the read length is very similar to
the other samples. To further investigate this, I also checked from
which populations these samples come from.

``` r
out_rl_names <- c(not151$sample_alias)
out_rl <- filter(HGDPcutoff, type == "te", ID %in% out_rl_names)
(unique(out_rl$Pop))
```

    ## [1] "Han"         "Tujia"       "Miao"        "NorthernHan" "Dai"        
    ## [6] "Lahu"        "She"         "Naxi"        "Tu"

The samples are distributed among 9 different populations, confirming
the non significance of these outliers on the final results.
