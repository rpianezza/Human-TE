SGDP - SNPs analysis
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
SGDP <- read_tsv("/Volumes/Temp1/rpianezza/SGDP/summary/SGDP_classified.tsv", col_names = c( "ID","sex","pop","country","type","familyname","length","reads","copynumber","batch", "superfamily", "shared_with"), skip=1)
```

    ## Rows: 470028 Columns: 12
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (9): ID, sex, pop, country, type, familyname, batch, superfamily, shared...
    ## dbl (3): length, reads, copynumber
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
SNP_step <- function(freq_matrix, metadata){
matrix <- read_tsv(freq_matrix)
samples <- metadata %>% select(ID) %>% distinct() %>% pull()
matrix <- matrix %>% add_column(.before = 1, ID=samples)
#write_tsv(matrix, "/Volumes/Temp1/rpianezza/SGDP/SNP/IDinverted.all.08.5000x.matrix.tsv")
}

SNP_step("/Volumes/Temp1/rpianezza/SGDP/SNP/inverted.all.08.5000x.matrix.tsv", SGDP)
```

    ## Rows: 276 Columns: 11048
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (11048): HERV9_te_79, HERV9_te_89, HERV9_te_90, HERV9_te_123, HERV9_te_1...
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
   
f <- f_pca_result$x %>% as_tibble() %>% add_column(.before = 1, ID=f_metadata$ID, sex=f_metadata$sex, pop=f_metadata$pop, country=f_metadata$country) %>% mutate(country= relevel(country, ref=c("Africa", "America", "Central Asia and Siberia", "East Asia", "South Asia", "West Eurasia", "Oceania"))) %>% as.data.frame() %>%
 ggplot(aes(x=PC1,y=PC2, color=country)) + geom_point(size=2) +
labs(x=paste0("PC1: ",round(f_var_explained[1]*100,1),"%"),
        y=paste0("PC2: ",round(f_var_explained[2]*100,1),"%")) + ggtitle("All the repetitive sequences - females") +
 theme(plot.title = element_text(hjust = 0.5))
   
  m <- m_pca_result$x %>% as_tibble() %>% add_column(.before = 1, ID=m_metadata$ID, sex=m_metadata$sex, pop=m_metadata$pop, country=m_metadata$country) %>% mutate(country= relevel(country, ref=c("Africa", "America", "Central Asia and Siberia", "East Asia", "South Asia", "West Eurasia", "Oceania"))) %>% as.data.frame() %>%
  ggplot(aes(x=PC1,y=PC2, color=country)) + geom_point(size=2) +
labs(x=paste0("PC1: ",round(m_var_explained[1]*100,1),"%"),
        y=paste0("PC2: ",round(m_var_explained[2]*100,1),"%")) + ggtitle("All the repetitive sequences - males") +
 theme(plot.title = element_text(hjust = 0.5))
 
ggarrange(f, m, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom", align = "hv", font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))
  }
```

`{r all}SNP_all("/Volumes/Temp1/rpianezza/SGDP/SNP/separed.all.08.5000x.matrix.tsv", SGDP)`

``` r
#SNP_all_hap <- function(freq_matrix, metadata){

#matrix <- read_csv(freq_matrix)

#m_metadata <- metadata %>% filter(sex=="male") %>% select(ID, sex, country, pop, haplotype) %>% distinct()
#males_matrix <- filter(matrix, ID %in% m_metadata$ID) %>% select(!(ID))
#m_pca_data <- males_matrix %>% select_if(negate(function(col) sd(col)==0))
#m_pca_result <- prcomp(m_pca_data, center = TRUE, scale = TRUE)
#m_var_explained <- m_pca_result$sdev^2/sum(m_pca_result$sdev^2)
   
#m_pca_result$x %>% as_tibble() %>% add_column(.before = 1, ID=m_metadata$ID, sex=m_metadata$sex, pop=m_metadata$pop, country=m_metadata$country, hap=m_metadata$haplotype) %>% as.data.frame() %>%
  #ggplot(aes(x=PC1,y=PC2, color=hap
            # )) + geom_point(size=2) +
   #labs(x=paste0("PC1: ",round(m_var_explained[1]*100,1),"%"),
       # y=paste0("PC2: ",round(m_var_explained[2]*100,1),"%")) + ggtitle("All the repetitive sequences - males") +
# theme(plot.title = element_text(hjust = 0.5))
 # }
```

``` r
#SNP_all_hap("/Volumes/Temp1/rpianezza/TE/SNP/try04/separed.all.08.5000x.matrix.tsv", metadata_hap)
```
