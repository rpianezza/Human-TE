Pipeline installation
================

## Installation

To install our pipeline on my local computer (vetgrid28), I did the
following:

I first create a **conda environment** where I installed all the
softwares needed in their correct versions. The environment location is:
/usr/local/Caskroom/miniconda/base/envs/human-te

    conda create --name human-te python==3.8.5
    conda activate human-te

After activating the env, I installed the following softwares:

    conda install -c bioconda bwa=0.7.17
    conda install -c conda-forge curl=7.71.1
    conda install -c bioconda mosdepth=0.3.1

Unfortunately, I encountered lot of problems in installing **samtools
version 1.10**, thus I kept the local version of the program. Later we
checked if the pipeline was producing the same output on the two
computers even with a different samtools version.

I then checked where each software was installed and substituted the
path of each software into the script **run.map.sh**.

    which md5
    which samtools
    which bwa
    which mosdepth
    which python
    which java
    which gzip

## Batch effect control

Here I analysed with the pipeline the same genome (cram file) on the two
computers (flo and ric) and checked if the two outputs were identical.

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
ric_mq10.mapstat <- read_tsv("/Users/rpianezza/Library/CloudStorage/OneDrive-Personal/PHD/test-2pc/ric/SAMEA3302822.mq10.mapstat")
```

    ## Warning: One or more parsing issues, see `problems()` for details

    ## Rows: 1710 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): summary, all_reads, 290500572
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
flo_mq10.mapstat <- read_tsv("/Users/rpianezza/Library/CloudStorage/OneDrive-Personal/PHD/test-2pc/flo/SAMEA3302822.mq10.mapstat")
```

    ## Warning: One or more parsing issues, see `problems()` for details

    ## Rows: 1710 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): summary, all_reads, 290500572
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
all.equal(ric_mq10.mapstat,flo_mq10.mapstat)
```

    ## [1] TRUE

``` r
ric_mq0.mapstat <- read_tsv("/Users/rpianezza/Library/CloudStorage/OneDrive-Personal/PHD/test-2pc/ric/SAMEA3302822.mq0.mapstat")
```

    ## Warning: One or more parsing issues, see `problems()` for details

    ## Rows: 1710 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): summary, all_reads, 290500572
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
flo_mq0.mapstat <- read_tsv("/Users/rpianezza/Library/CloudStorage/OneDrive-Personal/PHD/test-2pc/flo/SAMEA3302822.mq0.mapstat")
```

    ## Warning: One or more parsing issues, see `problems()` for details

    ## Rows: 1710 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): summary, all_reads, 290500572
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
all.equal(ric_mq0.mapstat,flo_mq0.mapstat)
```

    ## [1] TRUE

``` r
ric.mosdepth.global.dist <- read_tsv("/Users/rpianezza/Library/CloudStorage/OneDrive-Personal/PHD/test-2pc/ric/SAMEA3302822.mosdepth.global.dist.txt")
```

    ## Rows: 939035 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): LTR65_te
    ## dbl (2): 1142, 0.00
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
flo.mosdepth.global.dist <- read_tsv("/Users/rpianezza/Library/CloudStorage/OneDrive-Personal/PHD/test-2pc/flo/SAMEA3302822.mosdepth.global.dist.txt")
```

    ## Rows: 939035 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): LTR65_te
    ## dbl (2): 1142, 0.00
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
all.equal(ric.mosdepth.global.dist,flo.mosdepth.global.dist)
```

    ## [1] TRUE

``` r
ric.mosdepth.sum <- read_tsv("/Users/rpianezza/Library/CloudStorage/OneDrive-Personal/PHD/test-2pc/ric/SAMEA3302822.mosdepth.summary.txt")
```

    ## Rows: 1668 Columns: 6
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): chrom
    ## dbl (5): length, bases, mean, min, max
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
flo.mosdepth.sum <- read_tsv("/Users/rpianezza/Library/CloudStorage/OneDrive-Personal/PHD/test-2pc/flo/SAMEA3302822.mosdepth.summary.txt")
```

    ## Rows: 1668 Columns: 6
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): chrom
    ## dbl (5): length, bases, mean, min, max
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
all.equal(ric.mosdepth.sum,flo.mosdepth.sum)
```

    ## [1] TRUE

``` r
ric.sync <- read_tsv("/Users/rpianezza/Library/CloudStorage/OneDrive-Personal/PHD/test-2pc/ric/SAMEA3302822.mq10.sync.gz")
```

    ## Rows: 3124685 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): LTR65_te, t, 0:60:0:0:0:0
    ## dbl (1): 1
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
flo.sync <- read_tsv("/Users/rpianezza/Library/CloudStorage/OneDrive-Personal/PHD/test-2pc/flo/SAMEA3302822.mq10.sync.gz")
```

    ## Rows: 3124685 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): LTR65_te, t, 0:60:0:0:0:0
    ## dbl (1): 1
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
all.equal(ric.sync,flo.sync)
```

    ## [1] TRUE

``` r
ric.per_base <- read_tsv("/Users/rpianezza/Library/CloudStorage/OneDrive-Personal/PHD/test-2pc/ric/SAMEA3302822.per-base.bed.gz")
```

    ## Rows: 1763145 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): LTR65_te
    ## dbl (3): 0, 1, 60
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
flo.per_base <- read_tsv("/Users/rpianezza/Library/CloudStorage/OneDrive-Personal/PHD/test-2pc/flo/SAMEA3302822.per-base.bed.gz")
```

    ## Rows: 1763145 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): LTR65_te
    ## dbl (3): 0, 1, 60
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
all.equal(ric.per_base,flo.per_base)
```

    ## [1] TRUE

All the output files are equal. Note that this process has been did two
times, with another cram file (SAMEA3302884.cram), to further ensure the
pipeline robustness among the two computers. The result was the same.
