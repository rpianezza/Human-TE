SGDP - Data download and summary file creation
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
#(report <- read_tsv("/Volumes/Temp1/rpianezza/SGDP/filereport_read_run_PRJEB9586_tsv.txt", skip=1, col_names = c("run", "sample", "md5", "link")) %>% separate(link, into = c("fastq1", "fastq2"), sep = ";") %>% separate(md5, into = c("md5_1", "md5_2"), sep = ";") %>% select(!(run)) %>% group_by(sample) %>% filter(row_number(sample) == 1))

#write_tsv(report, "/Volumes/Temp1/rpianezza/SGDP/SGDP_report.tsv")

#distinct(report, sample)
#distinct(report, md5_2)
#group_by(report, sample) %>% filter(n()>1)
```

I first went to the repository where the **cram** files are stored
(<https://www.internationalgenome.org/data-portal/data-collection/sgdp>)
and I downloaded the list of files. 276 samples are present in the
repository.

Here I read the .tvs file downloaded and I create 3 files:

- `ftp_and_metadata`: it contains all the data that can be useful for us
- `ftp`: it contains only the url for downloading the files using
  **wget**
- `samples`: it contains only the **sample names**. Note that the sample
  names that are written in the cram files names (es. SAMEA3302884) are
  different from the sample names reported on the website (es. iran11).
  Here we create a file containing the sample names of the cram files.

``` r
ftp <- read_tsv("/Volumes/Temp1/rpianezza/SGDP/igsr_SGDP_files.tsv", skip=1, col_names = c("url","md5","data_collection","data_type","analysis_group","sample","population","data_reuse_policy")) %>% select(!(c(data_collection, data_type, analysis_group, data_reuse_policy)))
```

    ## Rows: 276 Columns: 8
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (8): url, md5, data_collection, data_type, analysis_group, sample, popul...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
ftp_only <- select(ftp, url)
sample_only <- select(ftp, url) %>% separate(url, into = c("ftp:", "/", "repo", "vol1", "ftp", "data_coll", "sgdd", "data", "pop", "sample", "alig"), sep="/") %>% select(sample)
```

    ## Warning: Expected 11 pieces. Additional pieces discarded in 276 rows [1, 2, 3,
    ## 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].

``` r
sample_md5 <- select(ftp, md5, url) %>% separate(url, into = c("ftp:", "/", "repo", "vol1", "ftp", "data_coll", "sgdd", "data", "pop", "sample", "alig"), sep="/") %>% select(sample, md5)
```

    ## Warning: Expected 11 pieces. Additional pieces discarded in 276 rows [1, 2, 3,
    ## 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].

``` r
write_tsv(ftp, "/Volumes/Temp1/rpianezza/SGDP/SGDP_ftp_and_metadata.tsv")
write_tsv(sample_md5, "/Volumes/Temp1/rpianezza/SGDP/SGDP_md5.tsv")
write_tsv(ftp_only, "/Volumes/Temp1/rpianezza/SGDP/SGDP_ftp.tsv")
write_tsv(sample_only, "/Volumes/Temp1/rpianezza/SGDP/SGDP_samples.tsv")
```

I then **split** the file `ftp` into small files containing 5 url once,
and then I download all the files contained in each file using **wget**.

    split -l 5 /Volumes/Temp2/riccardo/SGDP/data/SGDP_ftp.tsv /Volumes/Temp2/riccardo/SGDP/data/ftp_part_

    wget -i /Volumes/Temp2/riccardo/SGDP/data/SGDP_ftp_parts/ftp_part_aa -P /Volumes/Temp2/riccardo/SGDP/data/

I split also the file `samples` into small files (4 samples each) for
the **analysis pipeline**. These files can be directly used as input
files in the pipeline, if previously renamed removing all the character
after the sample code from the name (to match the sample codes in the
`sample` file).

    split -l 4 /Volumes/Temp2/riccardo/SGDP/data/SGDP_samples.tsv /Volumes/Temp2/riccardo/SGDP/data/samples_part_

    bash /Volumes/Temp2/human_TEs/human-te-dynamics-svn/dev/pipe-v1/run-map.sh  /Volumes/Temp2/riccardo/SGDP/data/samples_part_aa /Volumes/Temp2/riccardo/SGDP/data/ 

Here I take info from different files to create the `metadata` file
which contains all the useful information for our analysis. The files
can be downloaded from
<https://www.internationalgenome.org/data-portal/data-collection/sgdp>
(samples, populations, data).

``` r
sample_ID <- select(ftp, url, sample, population) %>% separate(url, into = c("ftp:", "/", "repo", "vol1", "ftp", "data_coll", "sgdd", "data", "pop", "biosample", "alig"), sep="/") %>% select(sample, biosample, population) %>% separate(population, into = c("pop", "description"), sep=" in") %>% select(!(description))
```

    ## Warning: Expected 11 pieces. Additional pieces discarded in 276 rows [1, 2, 3,
    ## 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].

``` r
coordinates <- read_tsv("/Volumes/Temp1/rpianezza/SGDP/coordinates.tsv", skip=1, col_names = c("nan", "population_code","pop","pop_description","latitude","longitude","nan", "country")) %>% select("pop","latitude","longitude","country") %>% separate(country, into = c("country", "sgdp"), sep="\\(") %>% select(!(sgdp))
```

    ## New names:
    ## Rows: 129 Columns: 11
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "\t" chr
    ## (7): nan...1, population_code, pop, pop_description, country, X9, X11 dbl (3):
    ## latitude, longitude, X10 lgl (1): nan...7
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `nan` -> `nan...1`
    ## • `nan` -> `nan...7`

``` r
sex <- read_tsv("/Volumes/Temp1/rpianezza/SGDP/igsr-SGDP_samples.tsv", skip=1, col_names = c("sample",  "sex",  "biosample")) %>% select(sex, sample)
```

    ## Rows: 276 Columns: 9
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (9): sample, sex, biosample, X4, X5, X6, X7, X8, X9
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
sample_and_coord <- inner_join(sample_ID, coordinates, by="pop")
(metadata <- inner_join(sample_and_coord, sex, by="sample") %>% arrange(pop) %>% relocate(sex, .before=pop) %>% relocate(country, .before=latitude))
```

    ## # A tibble: 276 × 7
    ##    sample    biosample    sex    pop       country               latit…¹ longi…²
    ##    <chr>     <chr>        <chr>  <chr>     <chr>                   <dbl>   <dbl>
    ##  1 abh107    SAMEA3302655 male   Abkhasian "West Eurasia "          43      41  
    ##  2 abh100    SAMEA3302693 male   Abkhasian "West Eurasia "          43      41  
    ##  3 HGDP01401 SAMEA3302728 female Adygei    "West Eurasia "          44      39  
    ##  4 HGDP01402 SAMEA3302656 male   Adygei    "West Eurasia "          44      39  
    ##  5 ALB212    SAMEA3302780 female Albanian  "West Eurasia "          41.3    19.8
    ##  6 Ale20     SAMEA3302642 female Aleut     "Central Asia and Si…    55.2   166  
    ##  7 Ale22     SAMEA3302773 male   Aleut     "Central Asia and Si…    55.2   166  
    ##  8 altai363p SAMEA3302753 male   Altaian   "Central Asia and Si…    51.9    86  
    ##  9 NA13607   SAMEA3302722 male   Ami       "East Asia "             22.8   121. 
    ## 10 NA13616   SAMEA3302623 male   Ami       "East Asia "             22.8   121. 
    ## # … with 266 more rows, and abbreviated variable names ¹​latitude, ²​longitude

``` r
write_tsv(metadata, "/Volumes/Temp1/rpianezza/SGDP/SGDP_metadata.tsv")
```

From **sample_part_ca** I start the pipeline on my local PC (ric). The
corresponding ftp file is **ftp_part_bp**, where the first sample in
**sample_part_ca** is stored (SAMEA3302678). Note that the files are not
synchronized, thus the ftp file contains other 3 links before the one we
are starting with on this pc.
