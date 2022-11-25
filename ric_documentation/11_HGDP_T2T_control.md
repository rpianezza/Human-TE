Telomer-to-telomer (T2T) Repeat Masker output analysis
================

Idea: run RepeatMasker (RM) on the T2T human genome assembly published
in 2022, then creat artificial short reads from the same genome and
compare the TE copynumbers obtained with the two methods.

First, with this command I run RM, using our customized RepBase
reference library as `- lib`.

    RepeatMasker -gccalc -s -cutoff 200 -no_is -nolow -norna -gff -u -pa 20 -lib /Volumes/Temp2/human_TEs/human-te-dynamics-svn/refg/reflibrary_humans_v6.2.fasta /Volumes/Temp2/riccardo/T2T/genome_assemblies_genome_fasta/ncbi-genomes-2022-11-18/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna

    # RepeatMasker version open-4.0.7
    # Search Engine: NCBI/RMBLAST [ 2.2.27+ ]

Explanation of the RepeatMasker output file:

- `SWscore` = Smith-Waterman score of the match, usually complexity
  adjusted
- `perc_div` = % substitutions in matching region compared to the
  consensus
- `perc_del` = % of bases opposite a gap in the query sequence (deleted
  bp)
- `perc_ins` = % of bases opposite a gap in the repeat consensus
  (inserted bp)
- `query_sequence` = name of query sequence
- `position_in_query_begin` = starting position of match in query
  sequence
- `position_in_query_end` = ending position of match in query sequence
- `position_in_query_left` = no. of bases in query sequence past the
  ending position of match
- `C` = match is with the Complement of the consensus sequence in the
  database
- `matching_repeat` = name of the matching interspersed repeat
- `repeat_class/family` = the class of the repeat
- `position_in_repeat_begin` = starting position of match in database
  sequence (using top-strand numbering)
- `position_in_repeat_end` = ending position of match in database
  sequence
- `position_in_repeat_left` = no. of bases in the repeat consensus
  sequence prior to beginning of the match (so 0 means that the match
  extended all the way to the end of the repeat consensus sequence)
- `ID` = estimated unique transposon (es. two segments of the same
  transposon could be separated by another insertion, thus these two
  sequences have the same ID)
- An asterisk (\*) in the final column indicates that there is a
  higher-scoring match whose domain partly (\<80%) includes the domain
  of this match.

Command to remove multiple spaces from the RM output and make it
readable in R:

    less path/rm.out | sed 's/  */ /g' | cut -c2- | > output

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
library(ggpubr)

(RM <- read_delim("/Users/rpianezza/TE/T2T/RM_complete/GCF_009914755.1_T2T-CHM13v2.0_genomic_mod.fna.out", delim = " ",skip = 3, col_names = c("SWscore", "perc_div", "perc_del", "perc_ins", "query_sequence", "position_in_query_begin", "position_in_query_end", "position_in_query_left",  "+/c", "matching_repeat", "repeat_class/family", "position_in_repeat_left", "position_in_repeat_begin", "position_in_repeat_end", "ID", "other_match")))
```

    ## Rows: 8353919 Columns: 16
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: " "
    ## chr (9): SWscore, query_sequence, position_in_query_left, +/c, matching_repe...
    ## dbl (7): perc_div, perc_del, perc_ins, position_in_query_begin, position_in_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

    ## # A tibble: 8,353,919 × 16
    ##    SWscore perc_…¹ perc_…² perc_…³ query…⁴ posit…⁵ posit…⁶ posit…⁷ `+/c` match…⁸
    ##    <chr>     <dbl>   <dbl>   <dbl> <chr>     <dbl>   <dbl> <chr>   <chr> <chr>  
    ##  1 234         4.8     0       0   NC_060…    1885    1926 (24838… C     TAR1_te
    ##  2 250         2.4     0       2.4 NC_060…    1897    1939 (24838… C     TAR1_te
    ##  3 5109       11.7     1.1     0.8 NC_060…    2710    3803 (24838… C     TAR1_te
    ##  4 1236       20.6     2.2     3.1 NC_060…    3484    3843 (24838… C     TAR1_te
    ##  5 964        25.2     4.3     2.1 NC_060…    3542    3910 (24838… C     TAR1_te
    ##  6 902        25.5     4.4     2.7 NC_060…    3602    3968 (24838… C     TAR1_te
    ##  7 3879       13.6     5.3     2.2 NC_060…    3663    4398 (24838… C     TAR1_te
    ##  8 990        26.2    12.4     1   NC_060…    4083    4533 (24838… C     LTR60B…
    ##  9 552        18       5.7     0   NC_060…    4534    4655 (24838… C     LTR60_…
    ## 10 1013       24.2     5       1.1 NC_060…    4718    5139 (24838… +     L1ME_O…
    ## # … with 8,353,909 more rows, 6 more variables: `repeat_class/family` <chr>,
    ## #   position_in_repeat_left <chr>, position_in_repeat_begin <dbl>,
    ## #   position_in_repeat_end <chr>, ID <dbl>, other_match <chr>, and abbreviated
    ## #   variable names ¹​perc_div, ²​perc_del, ³​perc_ins, ⁴​query_sequence,
    ## #   ⁵​position_in_query_begin, ⁶​position_in_query_end, ⁷​position_in_query_left,
    ## #   ⁸​matching_repeat

Here I filter the RM output putting a threshold on the maximum
**divergence** percentage among the query and the reference, as well as
a minimum **length** for the match. I also remove the brackets present
in some data to indicate that the sequence aligned to the complement of
the reference, which create some problems with the data type and which
are not interesting for this analysis.

As output of this chunk, I create a .csv file to analyze in Python using
the script *from_RM_to_copynumber*, which estimates the copynumber of
each sequence from the ref library found in the RM run. This is done by
the script by simply adding the length of each match of a sequence and
then dividing by the sequence length.

``` r
(RM_cutoff <- filter(RM, perc_div < 15, position_in_query_end-position_in_query_begin>150) %>% replace_na(list(other_match = "-")) %>% filter(other_match == "-") %>% mutate(position_in_repeat_begin = str_replace(position_in_repeat_begin, "\\(", "")) %>% mutate(position_in_repeat_begin = str_replace(position_in_repeat_begin, "\\)", "")) %>% mutate(position_in_repeat_left = str_replace(position_in_repeat_left, "\\(", "")) %>% mutate(position_in_repeat_left = str_replace(position_in_repeat_left, "\\)", "")) %>% mutate(position_in_query_left = str_replace(position_in_query_left, "\\(", "")) %>% mutate(position_in_query_left = str_replace(position_in_query_left, "\\)", "")) %>% mutate(position_in_repeat_end = str_replace(position_in_repeat_end, "\\(", "")) %>% mutate(position_in_repeat_end = str_replace(position_in_repeat_end, "\\)", "")) %>% mutate(position_in_query_begin = str_replace(position_in_query_begin, "\\(", "")) %>% mutate(position_in_query_begin = str_replace(position_in_query_begin, "\\)", "")) %>% mutate(position_in_query_end = str_replace(position_in_query_end, "\\(", "")) %>% mutate(position_in_query_end = str_replace(position_in_query_end, "\\)", "")) %>% arrange(matching_repeat))
```

    ## # A tibble: 1,169,025 × 16
    ##    SWscore perc_…¹ perc_…² perc_…³ query…⁴ posit…⁵ posit…⁶ posit…⁷ `+/c` match…⁸
    ##    <chr>     <dbl>   <dbl>   <dbl> <chr>   <chr>   <chr>   <chr>   <chr> <chr>  
    ##  1 4251        1.9     0       0.2 NC_060… 127766… 127769… 120618… +     6kbHsa…
    ##  2 5953        5       0.4     0.3 NC_060… 127770… 127772… 120614… C     6kbHsa…
    ##  3 201        11.3     6.3     5.6 NC_060… 208871… 208871… 395155… +     6kbHsa…
    ##  4 1105       14.6     1       0   NC_060… 693439… 693441… 173352… C     6kbHsa…
    ##  5 2985        6.8     1.2     0   NC_060… 132718… 132721… 109975… +     6kbHsa…
    ##  6 3406        7.4     1       0.4 NC_060… 132721… 132726… 109970… +     6kbHsa…
    ##  7 2922       10.5     3.6     3   NC_060… 132727… 132727… 109968… +     6kbHsa…
    ##  8 0558        4.8     1.3     0   NC_060… 161737… 161740… 809562… C     6kbHsa…
    ##  9 6574        2.4     0.2     0.1 NC_060… 957535… 957556… 105350… +     6kbHsa…
    ## 10 2956        3.7     0.3     0.1 NC_060… 957549… 957591… 105346… +     6kbHsa…
    ## # … with 1,169,015 more rows, 6 more variables: `repeat_class/family` <chr>,
    ## #   position_in_repeat_left <chr>, position_in_repeat_begin <chr>,
    ## #   position_in_repeat_end <chr>, ID <dbl>, other_match <chr>, and abbreviated
    ## #   variable names ¹​perc_div, ²​perc_del, ³​perc_ins, ⁴​query_sequence,
    ## #   ⁵​position_in_query_begin, ⁶​position_in_query_end, ⁷​position_in_query_left,
    ## #   ⁸​matching_repeat

``` r
write_csv(RM_cutoff, "/Users/rpianezza/TE/T2T/RM_cutoff.csv")
```

Here I read the output of the Python script.

``` r
RM_copynumbers <- read_csv("/Users/rpianezza/TE/T2T/RM_complete/copynumber_RM.csv", col_names = c("sequence", "copynumber"))
```

    ## Rows: 1467 Columns: 2
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (1): sequence
    ## dbl (1): copynumber
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
RM_scg <- filter(RM_copynumbers, str_detect(sequence, ".scg"))
(RM_scg_outliers <- filter(RM_scg, copynumber > 1.5))
```

    ## # A tibble: 6 × 2
    ##   sequence                      copynumber
    ##   <chr>                              <dbl>
    ## 1 chr1:2629757-2632979_scg            1.53
    ## 2 chr14:105475310-105478776_scg       2.92
    ## 3 chr16:765731-768798_scg             2.73
    ## 4 chr16:792523-795720_scg             2.06
    ## 5 chr16:89905008-89909174_scg         2.07
    ## 6 chr7:95409506-95411684_scg         10.2

``` r
RM_krab <- filter(RM_copynumbers, str_detect(sequence, ".krab"))
(RM_krab_outliers <- filter(RM_krab, copynumber > 2.5))
```

    ## # A tibble: 10 × 2
    ##    sequence        copynumber
    ##    <chr>                <dbl>
    ##  1 a_CDK8_20_krab        4.93
    ##  2 a_EWSR1_17_krab       5.49
    ##  3 a_FOXD3_0_krab        3.27
    ##  4 a_FOXO1_0_krab        2.88
    ##  5 a_PCBP2_29_krab       5.80
    ##  6 a_SRSF2_10_krab       2.92
    ##  7 a_YWHAB_9_krab        2.58
    ##  8 a_ZNF248_3_krab       4.48
    ##  9 a_ZNF302_4_krab       3.14
    ## 10 a_ZNF674_3_krab       4.88

Here I select only the TE, removing KRAB and SCG from the dataset. The
resulting CSV file can be used for other analyses.

``` r
(RM_cutoff_te <- filter(RM_cutoff, str_detect(matching_repeat, ".te")) %>% mutate(position_in_query_begin = as.integer(position_in_query_begin)) %>% arrange(matching_repeat, position_in_query_begin))
```

    ## # A tibble: 1,168,220 × 16
    ##    SWscore perc_…¹ perc_…² perc_…³ query…⁴ posit…⁵ posit…⁶ posit…⁷ `+/c` match…⁸
    ##    <chr>     <dbl>   <dbl>   <dbl> <chr>     <int> <chr>   <chr>   <chr> <chr>  
    ##  1 8749        8.3     0.6     0.9 NC_060…    2070 4815    101156… C     6kbHsa…
    ##  2 6479       14.8     1.4     0.3 NC_060…    4286 5535    101155… C     6kbHsa…
    ##  3 2599        5.1     0.4     0.5 NC_060…    8384 12657   101148… C     6kbHsa…
    ##  4 2470        5.2     0.3     0.6 NC_060…   14657 18937   101142… C     6kbHsa…
    ##  5 2530        5.1     0.3     0.6 NC_060…   20945 25225   101136… C     6kbHsa…
    ##  6 2606        5       0.3     0.6 NC_060…   27225 31505   101129… C     6kbHsa…
    ##  7 2610        5       0.3     0.6 NC_060…   33505 37785   101123… C     6kbHsa…
    ##  8 2510        5.2     0.3     0.6 NC_060…   39785 44065   101117… C     6kbHsa…
    ##  9 2436        5.2     0.3     0.6 NC_060…   46065 50345   101111… C     6kbHsa…
    ## 10 2620        5       0.3     0.6 NC_060…   52345 56625   101104… C     6kbHsa…
    ## # … with 1,168,210 more rows, 6 more variables: `repeat_class/family` <chr>,
    ## #   position_in_repeat_left <chr>, position_in_repeat_begin <chr>,
    ## #   position_in_repeat_end <chr>, ID <dbl>, other_match <chr>, and abbreviated
    ## #   variable names ¹​perc_div, ²​perc_del, ³​perc_ins, ⁴​query_sequence,
    ## #   ⁵​position_in_query_begin, ⁶​position_in_query_end, ⁷​position_in_query_left,
    ## #   ⁸​matching_repeat

``` r
write_csv(RM_cutoff, "/Users/rpianezza/TE/T2T/RM_cutoff_te.csv")
```

Using this command and the Python2 script *create-reads-for-human.py*, I
create artificial short reads from the T2T genome.

    python /Users/rpianezza/TE/human-te-dynamics-svn/scripts/create-reads-for-human.py --fasta /Users/rpianezza/TE/T2T/genome_assemblies_genome_fasta/ncbi-genomes-2022-11-22/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna --coverage 30 --read-length 150 --output /Users/rpianezza/TE/T2T/T2T-artificial-reads/reads.fastq.gz --method uniform

The artificial reads are then processed into the pipeline to estimate
copynumbers. Here I read the main output file and I create a comparison
between the two methods copynumber estimates: RepeatMasker (`RM`) and
the pipeline run on artificial short reads (`pipeline`).

``` r
artificial <- read_delim("/Users/rpianezza/TE/T2T/T2T-artificial-reads/T2T_artificial_reads.mq10.mapstat", delim = "\t", skip = 8, col_names = c("type", "familyname", "length", "reads", "copynumber"))
```

    ## Rows: 1703 Columns: 5
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): type, familyname
    ## dbl (3): length, reads, copynumber
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
art_subset <- select(artificial, familyname, copynumber) %>% arrange(familyname)
RM_subset <- RM_copynumbers %>% mutate(sequence = str_replace(sequence, "_te", "")) %>% mutate(sequence = str_replace(sequence, "_krab", "")) %>% mutate(sequence = str_replace(sequence, "_scg", "")) %>% mutate(sequence = str_replace(sequence, "_scgx", "")) %>% rename(familyname = sequence)

cn_comparison <- inner_join(art_subset, RM_subset, by = "familyname") %>% rename(pipeline_copynumber = copynumber.x) %>% rename(RM_copynumber = copynumber.y)

ggplot(cn_comparison, aes(x=log(pipeline_copynumber), y=log(RM_copynumber))) +
  geom_point(size=1) +
  geom_smooth(method="lm",color="grey", se=F) +
  ylab("RepeatMasker copynumber estimates (log)") + xlab("Pipeline copynumber estimates (log)") +
  stat_regline_equation(label.y = 13, aes(label = ..rr.label..), size=5)
```

    ## `geom_smooth()` using formula 'y ~ x'

    ## Warning: Removed 3 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 3 rows containing non-finite values (stat_regline_equation).

![](11_HGDP_T2T_control_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Here I check if there is a correlation between the number of reads (not
normalized) aligned by the pipeline and the number of hits of RM for
each sequence in the reference library.

``` r
(RM_reads <- group_by(RM, matching_repeat) %>% count() %>% mutate(matching_repeat = str_replace(matching_repeat, "_te", "")) %>% mutate(matching_repeat = str_replace(matching_repeat, "_krab", "")) %>% mutate(matching_repeat = str_replace(matching_repeat, "_scg", "")) %>% mutate(matching_repeat = str_replace(matching_repeat, "_scgx", "")) %>% rename(familyname = matching_repeat))
```

    ## # A tibble: 1,703 × 2
    ## # Groups:   familyname [1,703]
    ##    familyname         n
    ##    <chr>          <int>
    ##  1 6kbHsap        20241
    ##  2 a_AC067968.1_2    29
    ##  3 a_AC092835.1_5   132
    ##  4 a_CDK8_20         36
    ##  5 a_CHD3_3          19
    ##  6 a_CLTC_6          18
    ##  7 a_DCT_12          77
    ##  8 a_DNAJA2_11       53
    ##  9 a_E2F1_2         199
    ## 10 a_EGR1_4          54
    ## # … with 1,693 more rows

``` r
(pipeline_reads <- select(artificial, familyname, reads) %>% arrange(familyname))
```

    ## # A tibble: 1,703 × 2
    ##    familyname       reads
    ##    <chr>            <dbl>
    ##  1 6kbHsap        344075.
    ##  2 a_AC067968.1_2    312.
    ##  3 a_AC092835.1_5    296.
    ##  4 a_CDK8_20         487.
    ##  5 a_CHD3_3          232.
    ##  6 a_CLTC_6          276.
    ##  7 a_DCT_12          496.
    ##  8 a_DNAJA2_11       318.
    ##  9 a_E2F1_2          299.
    ## 10 a_EGR1_4          329.
    ## # … with 1,693 more rows

``` r
reads_comparison <- inner_join(pipeline_reads, RM_reads, by = "familyname") %>% rename(pipeline_reads = reads, RM_reads = n)

ggplot(reads_comparison, aes(x=log(pipeline_reads), y=log(RM_reads))) +
  geom_point(size=1) +
  geom_smooth(method="lm",color="grey", se=F) +
  ylab("RepeatMasker read number (log)") + xlab("Pipeline read number (log)") +
  stat_regline_equation(label.y = 13, aes(label = ..rr.label..), size=5)
```

    ## `geom_smooth()` using formula 'y ~ x'

    ## Warning: Removed 40 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 40 rows containing non-finite values (stat_regline_equation).

![](11_HGDP_T2T_control_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Here I calculate the mean copynumber for each `familyname` among all the
HGDP samples and I look for a correlation between these values and the
copynumber estimated by the pipeline in the T2T genome, to see if the
values are in line with the previous analysis.

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
names(HGDPcutoff)<-c("ID","pop","sex","country","type","familyname","length","reads","copynumber","batch")

mean_cn <- group_by(HGDPcutoff, familyname) %>% summarise(min = min(copynumber), mean = mean(copynumber), max = max(copynumber))

(T2T_HGDP <- inner_join(mean_cn, art_subset, by="familyname") %>% rename(T2T = copynumber))
```

    ## # A tibble: 1,684 × 5
    ##    familyname         min    mean    max     T2T
    ##    <chr>            <dbl>   <dbl>  <dbl>   <dbl>
    ##  1 6kbHsap        178.    305.    434.   285.   
    ##  2 a_AC067968.1_2   0.953   1.31    1.63   1.29 
    ##  3 a_AC092835.1_5   0.941   1.22    1.51   1.17 
    ##  4 a_CDK8_20        1.63    2.07    2.59   1.99 
    ##  5 a_CHD3_3         0.605   1.02    1.39   1.06 
    ##  6 a_CLTC_6         0.736   1.05    1.40   0.992
    ##  7 a_DCT_12         0.789   0.984   1.61   0.994
    ##  8 a_DNAJA2_11      0.774   1.05    1.42   0.993
    ##  9 a_E2F1_2         0.812   1.02    1.25   0.993
    ## 10 a_EGR1_4         0.967   1.25    1.57   1.24 
    ## # … with 1,674 more rows

``` r
ggplot(T2T_HGDP, aes(x=log(mean), y=log(T2T))) +
  geom_point(size=1) +
  geom_smooth(method="lm",color="grey")+
  ylab("HGDP mean copynumber (log)") + xlab("T2T copynumber (log)") +
  stat_regline_equation(label.y = 13, aes(label = ..rr.label..), size=5)
```

    ## `geom_smooth()` using formula 'y ~ x'

    ## Warning: Removed 21 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 21 rows containing non-finite values (stat_regline_equation).

![](11_HGDP_T2T_control_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->
