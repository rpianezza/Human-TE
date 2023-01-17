HGDP - SNPs analysis
================

All the files are present in the folders
<https://sourceforge.net/p/human-te-dynamics/data/HEAD/tree/raw-data/>.

`HGDP01382-Adygei.mq10.sync.gz` was NOT present in the folders and was
manually added from the computer **vetgrid27**.

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
library(ggpubr)
library(janitor)
```

    ## 
    ## Attaching package: 'janitor'
    ## 
    ## The following objects are masked from 'package:stats':
    ## 
    ##     chisq.test, fisher.test

``` r
HGDP <- read_tsv("/Volumes/Temp1/rpianezza/TE/summary-HGDP/HGDP_cutoff_classified.tsv", col_names = c( "ID","pop","sex","country","type","familyname","length","reads","copynumber","batch", "superfamily", "shared_with"), skip=1)
```

    ## Rows: 1396835 Columns: 12
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (9): ID, pop, sex, country, type, familyname, batch, superfamily, shared...
    ## dbl (3): length, reads, copynumber
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
sync_files <- read_tsv("/Volumes/Temp1/rpianezza/TE/SNP/file_list", col_names = "ID") %>% separate(ID, into = c("ID", "desc"), sep="-") %>% select(ID) %>% distinct(ID)
```

    ## Rows: 828 Columns: 1
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): ID
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
ID <- distinct(HGDP, ID)
identical(ID, sync_files)
```

    ## [1] TRUE

    cd /Volumes/Temp1/rpianezza/TE/SNP/HGDP_data
    gunzip *gz

`HGDP00902-Russian.mq10.sync.gz` and `HGDP01202-Mandenka.mq10.sync` were
corrupted, so I changed them with their version from the folder
*OLD_flo_vetlinux*.

Here I modify the files by adding a 5th column containing the sample
name, to differentiate file between each other prior to merging.

    mkdir modified_files &&
    for i in *.sync; do
      awk 'NR==1{ sub(/\.mq10.sync$/, "", FILENAME) } { $5=FILENAME }1' "$i" |
        column -t > "./modified_files/$i"
    done

Merge files. After creating the file I realized that the columns are
divided by multiple spaces. To remove the multiple spaces leaving only
one space between them, I ran the second commmand. These commands are
ideally merged together with a pipe.

    cat /path/to/files/*.sync >> sync_summary.sync

    less /Volumes/Temp1/rpianezza/TE/SNP/HGDP_data/sync_summary.sync | sed 's/  */ /g' | cut -c2- | > /Volumes/Temp1/rpianezza/TE/SNP/HGDP_data/sync_summary2.0.sync 

After using `gzip` on the summary file, I sent all the zipped files to
Robert’s server with:

    # The summary file:
    rsync -v -e ssh /Volumes/Temp1/rpianezza/TE/SNP/HGDP_data/sync_summary2.0.sync.gz student@roco.local:/home/student/sync_summary2.0.sync.gz

    # The single files:
    rsync -v -e ssh /Volumes/Temp1/rpianezza/TE/SNP/HGDP_data/modified_files/*.gz student@roco.local:/mnt/curana/riccardo/all-files/

Robert processed the files and created a matrix using a custom script
“sync2matrix.py”.

    python3 /home/robert/dev/human-te-dynamics-svn/scripts/sync2matrix.py --fai /home/robert/dev/human-te-dynamics-svn/refg/reflibrary_humans_v6.2.fasta.fai --sync summary.sync.gz |gzip -c > allte.matrix.gz

I take the file and slice extract the data for the TEs of interest,
creating a single file for each TE.

    scp student@roco.local:/mnt/curana/riccardo/allte.matrix.gz /Volumes/Temp1/rpianezza/TE/SNP

    head -1 /Volumes/Temp1/rpianezza/TE/SNP/allte.matrix > /Volumes/Temp1/rpianezza/TE/SNP/header.matrix
    # Manually added: "familyname" and "position" at the beginning of "header.matrix" using vi.

    # For all the investigated TEs (L1PB1, L1PA16, L1PA3, MER2, MLT2A1, L1PA7_5):
    less /Volumes/Temp1/rpianezza/TE/SNP/allte.matrix| grep 'MER2_te' > /Volumes/Temp1/rpianezza/TE/SNP/MER2.matrix

    cat /Volumes/Temp1/rpianezza/TE/SNP/header.matrix /Volumes/Temp1/rpianezza/TE/SNP/MER2.matrix >> /Volumes/Temp1/rpianezza/TE/SNP/matrixes/MER2.matrix

The single TE matrix is further processed with another custom script
“frequency_matrix.ipynb” which does:

- 1)  Remove all the bases of the TE for which we don’t have at least
      1000 reads covering that position in each sample. In other words,
      we remove all the bases with **coverage \< 1000**.

- 2)  Remove all the **non-polymorphic** positions. I filtered out the
      positions for which the major allele was \> 99%.

- 3)  Convert the raw counts for each base into **frequencies**.

This function creates the PCA, taking as inputs:

- **TE**: the name of the TE (es. “L1PB1”).
- **freq_matrix**: the path of the frequency matrix created with the
  Python script.
- **metadata**: the tibble containing “ID”, “sex” and “country” of the
  analysed data.

``` r
SNP_PCA <- function(TE, freq_matrix, metadata){
  matrix <- read_delim(freq_matrix) %>% select(!("...1")) 
  
  positions <- select(matrix, position) %>% pull
  bases <- paste(positions, "_", sep="")
  
  reverted <- as.data.frame(t(matrix)) %>% as_tibble() %>% filter(!row_number()==1) %>% filter(!row_number()==1)
  colnames(reverted) <- bases
  
  ID_cont_sex <- select(metadata, ID, sex, country, pop) %>% distinct()
  males <- reverted %>% add_column(.before = 1, ID=ID_cont_sex$ID, sex=ID_cont_sex$sex) %>% filter(sex=="male") %>% select(!(c(ID, sex)))
  females <- reverted %>% add_column(.before = 1, ID=ID_cont_sex$ID, sex=ID_cont_sex$sex) %>% filter(sex=="female")%>% select(!(c(ID, sex)))
  
  f_separated=females
  m_separated=males
  
  for (pos in bases){
  f_separated <- f_separated %>% separate(pos, into = c(paste0(pos,"A"), paste0(pos,"T"), paste0(pos,"C"), paste0(pos,"G"), paste0(pos,"N"), paste0(pos,"D")), sep=":")
  
  m_separated <- m_separated %>% separate(pos, into = c(paste0(pos,"A"), paste0(pos,"T"), paste0(pos,"C"), paste0(pos,"G"), paste0(pos,"N"), paste0(pos,"D")), sep=":")
  }
  
  f_separated <- f_separated %>% type_convert()
  f_all_col <- colnames(f_separated)
  m_separated <- m_separated %>% type_convert()
  m_all_col <- colnames(m_separated)
  
  f_name <- f_separated %>% dplyr::summarise(across(c(f_all_col), sd, na.rm = TRUE)) %>% select_if(~ any(. > 0))
  f_names <- colnames(f_name)
  m_name <- m_separated %>% dplyr::summarise(across(c(m_all_col), sd, na.rm = TRUE)) %>% select_if(~ any(. > 0))
  m_names <- colnames(m_name)
  
  f_pca_data <- select(f_separated, c(f_names))
f_pca_result <- prcomp(f_pca_data, center = TRUE, scale = TRUE)
  m_pca_data <- select(m_separated, c(m_names))
  m_pca_result <- prcomp(m_pca_data, center = TRUE, scale = TRUE)
 
  f_var_explained <- f_pca_result$sdev^2/sum(f_pca_result$sdev^2)
  m_var_explained <- m_pca_result$sdev^2/sum(m_pca_result$sdev^2)
   
  f_metadata <- ID_cont_sex %>% filter(sex=="female")
 m_metadata <- ID_cont_sex %>% filter(sex=="male")
 
 f_pca_data
   
#   f <- f_pca_result$x %>% as_tibble() %>% add_column(.before = 1, ID=f_metadata$ID, sex=f_metadata$sex, pop=f_metadata$pop, country=f_metadata$country) %>% as.data.frame() %>%
#   ggplot(aes(x=PC1,y=PC2, color=country)) + geom_point(size=2) +
#   labs(x=paste0("PC1: ",round(f_var_explained[1]*100,1),"%"),
#       y=paste0("PC2: ",round(f_var_explained[2]*100,1),"%")) +
#  theme(legend.position="bottom") + ggtitle(paste0(TE, " - females")) +
#  theme(plot.title = element_text(hjust = 0.5))
#   
#   m <- m_pca_result$x %>% as_tibble() %>% add_column(.before = 1, ID=m_metadata$ID, sex=m_metadata$sex, pop=m_metadata$pop, country=m_metadata$country) %>% as.data.frame() %>%
#   ggplot(aes(x=PC1,y=PC2, color=country)) + geom_point(size=2) +
#   labs(x=paste0("PC1: ",round(m_var_explained[1]*100,1),"%"),
#        y=paste0("PC2: ",round(m_var_explained[2]*100,1),"%")) +
#   theme(legend.position="bottom") + ggtitle(paste0(TE, " - males")) +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# ggarrange(f, m, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom", align = "hv", font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))
  }
```

``` r
#SNP_PCA("L1PA3", "/Volumes/Temp1/rpianezza/TE/SNP/frequency_matrixes/L1PA3.matrix", HGDP)
```

``` r
#SNP_PCA("L1PB1", "/Volumes/Temp1/rpianezza/TE/SNP/frequency_matrixes/L1PB1.matrix", HGDP)
```

``` r
#SNP_PCA("L1PA7_5", "/Volumes/Temp1/rpianezza/TE/SNP/frequency_matrixes/L1PA7_5.matrix", HGDP)
```

``` r
#SNP_PCA("MLT2A1", "/Volumes/Temp1/rpianezza/TE/SNP/frequency_matrixes/MLT2A1.matrix", HGDP)
```

``` r
#SNP_PCA("L1PA16", "/Volumes/Temp1/rpianezza/TE/SNP/frequency_matrixes/L1PA16.matrix", HGDP)
```

``` r
#SNP_PCA("MER2", "/Volumes/Temp1/rpianezza/TE/SNP/frequency_matrixes/MER2.matrix", HGDP)
```

``` r
#SNP_PCA("ALU", "/Volumes/Temp1/rpianezza/TE/SNP/frequency_matrixes/ALU.matrix", HGDP)
#SNP_PCA("L2", "/Volumes/Temp1/rpianezza/TE/SNP/frequency_matrixes/L2.matrix", HGDP)
#SNP_PCA("L1PA4", "/Volumes/Temp1/rpianezza/TE/SNP/frequency_matrixes/L1PA4.matrix", HGDP)
#SNP_PCA("L1PA6", "/Volumes/Temp1/rpianezza/TE/SNP/frequency_matrixes/L1PA6.matrix", HGDP)
#SNP_PCA("L1PREC2", "/Volumes/Temp1/rpianezza/TE/SNP/frequency_matrixes/L1PREC2.matrix", HGDP)
#SNP_PCA("SVA_A", "/Volumes/Temp1/rpianezza/TE/SNP/frequency_matrixes/SVA_A.matrix", HGDP)
#SNP_PCA("ALR_", "/Volumes/Temp1/rpianezza/TE/SNP/frequency_matrixes/ALR_.matrix", HGDP)
#SNP_PCA("ALR1", "/Volumes/Temp1/rpianezza/TE/SNP/frequency_matrixes/ALR1.matrix", HGDP)
```

ALU SWAsia outliers: Kalash, Brahui, Balochi

## Y chromosme haplotypes

``` r
hap1 <- read_csv("/Volumes/Temp1/rpianezza/TE/Y-chromosome/males-Y-chromosome.txt")
```

    ## Rows: 553 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (3): ID, haplotype, hap_group
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
id <- HGDP %>% filter(sex=="male") %>% select(ID) %>% distinct()

hap2 <- inner_join(hap1, id, by="ID")

(metadata_hap <- inner_join(hap2, HGDP, by="ID"))
```

    ## # A tibble: 932,910 × 14
    ##    ID     haplo…¹ hap_g…² pop   sex   country type  famil…³ length reads copyn…⁴
    ##    <chr>  <chr>   <chr>   <chr> <chr> <chr>   <chr> <chr>    <dbl> <dbl>   <dbl>
    ##  1 HGDP0… R1a     7       Brah… male  Centra… scg   chr1:9…   5136 1105.   0.861
    ##  2 HGDP0… R1a     7       Brah… male  Centra… scg   chr1:1…   3064  832.   1.09 
    ##  3 HGDP0… R1a     7       Brah… male  Centra… scg   chr1:1…   3239  901.   1.11 
    ##  4 HGDP0… R1a     7       Brah… male  Centra… scg   chr1:1…   4035 1102.   1.09 
    ##  5 HGDP0… R1a     7       Brah… male  Centra… scg   chr1:1…   2500  733.   1.17 
    ##  6 HGDP0… R1a     7       Brah… male  Centra… scg   chr1:1…   2599  580.   0.894
    ##  7 HGDP0… R1a     7       Brah… male  Centra… scg   chr1:1…   2124  477.   0.899
    ##  8 HGDP0… R1a     7       Brah… male  Centra… scg   chr1:2…   6284 1527.   0.973
    ##  9 HGDP0… R1a     7       Brah… male  Centra… scg   chr1:2…   3222  889.   1.10 
    ## 10 HGDP0… R1a     7       Brah… male  Centra… scg   chr1:3…   3698  868.   0.940
    ## # … with 932,900 more rows, 3 more variables: batch <chr>, superfamily <chr>,
    ## #   shared_with <chr>, and abbreviated variable names ¹​haplotype, ²​hap_group,
    ## #   ³​familyname, ⁴​copynumber

``` r
SNP_PCA_hap <- function(TE, freq_matrix, metadata, subset){
  samples <- metadata %>% select(ID) %>% distinct() %>% pull()
  matrix <- read_delim(freq_matrix) %>% select(!("...1"))
  
  positions <- select(matrix, position) %>% pull
  bases <- paste(positions, "_", sep="")
  males <- select(subset, ID, country, pop, haplotype, hap_group) %>% distinct()
    
  reverted <- as.data.frame(t(matrix)) %>% as_tibble() %>% filter(!row_number()==1) %>% filter(!row_number()==1)
  colnames(reverted) <- bases
  
  separated=reverted
  for (pos in bases){
  separated <- separated %>% separate(pos, into = c(paste0(pos,"A"), paste0(pos,"T"), paste0(pos,"C"), paste0(pos,"_G"), paste0(pos,"N"), paste0(pos,"D")), sep=":")
  }
  
  md <- select(metadata, ID, sex, country, pop) %>% distinct()
  
  separated <- separated %>% type_convert() %>% add_column(.before = 1, ID=md$ID) %>% filter(ID %in% subset$ID) %>% select(!(ID))
  all_col <- colnames(separated)
  
  name <- separated %>% dplyr::summarise(across(c(all_col), sd, na.rm = TRUE)) %>% select_if(~ any(. > 0))
  names <- colnames(name)
  
  pca_data <- select(separated, c(names))
  pca_result <- prcomp(pca_data, center = TRUE, scale = TRUE)
  var_explained <- pca_result$sdev^2/sum(pca_result$sdev^2)

  ppp <- pca_result$x %>% as_tibble() %>% add_column(.before = 1, ID=males$ID, pop=males$pop, country=males$country, hap=males$haplotype, hap_group=males$hap_group) %>% as.data.frame() %>%
  ggplot(aes(x=PC1,y=PC2, color=hap_group, label=hap)) + geom_point(size=2) + geom_text(hjust=0, vjust=0, size=2) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="bottom") + ggtitle(paste0(TE, " - males")) +
  theme(plot.title = element_text(hjust = 0.5))
 ppp
   }
```

``` r
#SNP_PCA_hap("L1PB1", "/Volumes/Temp1/rpianezza/TE/SNP/frequency_matrixes/L1PB1.matrix", HGDP, metadata_hap)
#SNP_PCA_hap("L1PA3", "/Volumes/Temp1/rpianezza/TE/SNP/frequency_matrixes/L1PA3.matrix", HGDP, metadata_hap)
#SNP_PCA_hap("L1PA16", "/Volumes/Temp1/rpianezza/TE/SNP/frequency_matrixes/L1PA16.matrix", HGDP, metadata_hap)
#SNP_PCA_hap("L1PA7_5", "/Volumes/Temp1/rpianezza/TE/SNP/frequency_matrixes/L1PA7_5.matrix", HGDP, metadata_hap)
#SNP_PCA_hap("MER2", "/Volumes/Temp1/rpianezza/TE/SNP/frequency_matrixes/MER2.matrix", HGDP, metadata_hap)
#SNP_PCA_hap("MLT2A1", "/Volumes/Temp1/rpianezza/TE/SNP/frequency_matrixes/MLT2A1.matrix", HGDP, metadata_hap)
```

``` r
SNP_step <- function(freq_matrix, metadata){
matrix <- read_tsv(freq_matrix)
samples <- metadata %>% select(ID) %>% distinct() %>% pull()
matrix <- matrix %>% add_column(.before = 1, ID=samples)
#write_tsv(matrix, "/Volumes/Temp1/rpianezza/TE/SNP/try04/IDinverted.all.08.5000x.matrix.tsv")
}

SNP_step("/Volumes/Temp1/rpianezza/TE/SNP/try04/inverted.all.08.5000x.matrix.tsv", HGDP)
```

    ## Rows: 828 Columns: 11580
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (11580): HERV9_te_78, HERV9_te_79, HERV9_te_89, HERV9_te_90, HERV9_te_12...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
SNP_all <- function(freq_matrix, metadata){

matrix <- read_csv(freq_matrix)
  
f_metadata <- metadata %>% filter(sex=="female") %>% select(ID, sex, country, pop) %>% distinct()
m_metadata <- metadata %>% filter(sex=="male") %>% select(ID, sex, country, pop) %>% distinct()
males_matrix <- filter(matrix, ID %in% m_metadata$ID) %>% select(!(ID))
females_matrix <- filter(matrix, ID %in% f_metadata$ID) %>% select(!(ID))

f_pca_data <- females_matrix %>% select_if(negate(function(col) sd(col)==0))
m_pca_data <- males_matrix %>% select_if(negate(function(col) sd(col)==0))

f_pca_result <- prcomp(f_pca_data, center = TRUE, scale = TRUE)
m_pca_result <- prcomp(m_pca_data, center = TRUE, scale = TRUE)
  
f_var_explained <- f_pca_result$sdev^2/sum(f_pca_result$sdev^2)
m_var_explained <- m_pca_result$sdev^2/sum(m_pca_result$sdev^2)
   
f <- f_pca_result$x %>% as_tibble() %>% add_column(.before = 1, ID=f_metadata$ID, sex=f_metadata$sex, pop=f_metadata$pop, country=f_metadata$country) %>% as.data.frame() %>%
 ggplot(aes(x=PC1,y=PC2, color=country)) + geom_point(size=2) +
labs(x=paste0("PC1: ",round(f_var_explained[1]*100,1),"%"),
        y=paste0("PC2: ",round(f_var_explained[2]*100,1),"%")) + ggtitle("All the repetitive sequences - females") +
 theme(plot.title = element_text(hjust = 0.5))
   
  m <- m_pca_result$x %>% as_tibble() %>% add_column(.before = 1, ID=m_metadata$ID, sex=m_metadata$sex, pop=m_metadata$pop, country=m_metadata$country) %>% as.data.frame() %>%
  ggplot(aes(x=PC1,y=PC2, color=country)) + geom_point(size=2) +
   labs(x=paste0("PC1: ",round(m_var_explained[1]*100,1),"%"),
        y=paste0("PC2: ",round(m_var_explained[2]*100,1),"%")) + ggtitle("All the repetitive sequences - males") +
 theme(plot.title = element_text(hjust = 0.5))
 
ggarrange(f, m, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom", align = "hv", font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))
  }
```

``` r
SNP_all("/Volumes/Temp1/rpianezza/TE/SNP/try04/separed.all.08.5000x.matrix.tsv", HGDP)
```

    ## Rows: 828 Columns: 69481
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr     (1): ID
    ## dbl (69480): HERV9_te_78A, HERV9_te_78T, HERV9_te_78C, HERV9_te_78G, HERV9_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](13_HGDP_SNP_files/figure-gfm/all-1.png)<!-- -->

``` r
SNP_all_hap <- function(freq_matrix, metadata){

matrix <- read_csv(freq_matrix)

m_metadata <- metadata %>% filter(sex=="male") %>% select(ID, sex, country, pop, haplotype) %>% distinct()
males_matrix <- filter(matrix, ID %in% m_metadata$ID) %>% select(!(ID))
m_pca_data <- males_matrix %>% select_if(negate(function(col) sd(col)==0))
m_pca_result <- prcomp(m_pca_data, center = TRUE, scale = TRUE)
m_var_explained <- m_pca_result$sdev^2/sum(m_pca_result$sdev^2)
   
m_pca_result$x %>% as_tibble() %>% add_column(.before = 1, ID=m_metadata$ID, sex=m_metadata$sex, pop=m_metadata$pop, country=m_metadata$country, hap=m_metadata$haplotype) %>% as.data.frame() %>%
  ggplot(aes(x=PC1,y=PC2, color=hap
             )) + geom_point(size=2) +
   labs(x=paste0("PC1: ",round(m_var_explained[1]*100,1),"%"),
        y=paste0("PC2: ",round(m_var_explained[2]*100,1),"%")) + ggtitle("All the repetitive sequences - males") +
 theme(plot.title = element_text(hjust = 0.5))
  }
```

``` r
SNP_all_hap("/Volumes/Temp1/rpianezza/TE/SNP/try04/separed.all.08.5000x.matrix.tsv", metadata_hap)
```

    ## Rows: 828 Columns: 69481
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr     (1): ID
    ## dbl (69480): HERV9_te_78A, HERV9_te_78T, HERV9_te_78C, HERV9_te_78G, HERV9_t...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

![](13_HGDP_SNP_files/figure-gfm/all-hap-1.png)<!-- -->
