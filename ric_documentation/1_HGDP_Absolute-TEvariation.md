HGDP - Analyzing the TEs with the highest absolute difference between
their minimum and maximum abundance estimate
================

This is the first script written by me, Riccardo. This script works with
the version of the HGDP dataset created in Script 2 by Florian. We
perform all analyses separately for males and females, as we established
that there are significant differences between the sexes in Script 5
from Florian.

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
library("ggpubr")
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

females <- filter(HGDPcutoff, sex == 'female')
males <- filter(HGDPcutoff, sex == 'male')
```

The idea of the first part of this code is to plot the **general
distribution of TE abundance** estimates in a boxplot. First I will
perform the analysis on the female subset, later on the male subset
using the same settings.

Looking at this plot for all 965 TEs seems excessive (and not
informative), so I created a subset. First, I identify the **minimum**,
**mean** and **maximum** value for each TE sequence. This code chunk
creates a data frame containing respectively the name of the TE and its
minimum, mean and maximum value in the female subset of the data set.
Everything containing minimum, mean and maximum value for each TE is
named from now `MMM`.

``` r
f_MMM <- filter(HGDPcutoff, type == "te", sex == "female") %>% group_by(familyname) %>% summarise(min = min(copynumber), mean = mean(copynumber), max = max(copynumber))

m_MMM <- filter(HGDPcutoff, type == "te", sex == "male") %>% group_by(familyname) %>% summarise(min = min(copynumber), mean = mean(copynumber), max = max(copynumber))
```

The criteria to select TEs was an **absolute value** of maximum copy
number vs minimum copy number as a threshold. In this way, the idea is
to keep only informative observations in this subset. For this analysis,
I set the threshold at `200`.

In the two sections below I create the two subsets for females and
males, I order the subsets by descent absolute difference and plot them
as boxplots.

## Females

``` r
f_outliers_names <- mutate(f_MMM, diff = max-min) %>% filter(diff>200 & diff<Inf)

f_outliers <- filter(HGDPcutoff, familyname %in% f_outliers_names$familyname, type == "te", sex == "female")
f_outliers <- f_outliers[order(f_outliers$copynumber,decreasing=T),]
f_outliers$familyname<-factor(f_outliers$familyname,levels=unique(f_outliers$familyname))

ggplot(f_outliers, aes(x=familyname, y=log(copynumber))) + geom_boxplot(notch=F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
```

![](1_HGDP_Absolute-TEvariation-copy-optimized_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## Males

``` r
m_outliers_names <- mutate(m_MMM, diff = max-min) %>% filter(diff>200 & diff<Inf)

m_outliers <- filter(HGDPcutoff, familyname %in% m_outliers_names$familyname, type == "te", sex == "male")
m_outliers <- m_outliers[order(m_outliers$copynumber,decreasing=T),]
m_outliers$familyname<-factor(m_outliers$familyname,levels=unique(m_outliers$familyname))

ggplot(m_outliers, aes(x=familyname, y=log(copynumber))) + geom_boxplot(notch=F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
```

![](1_HGDP_Absolute-TEvariation-copy-optimized_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

## Details for each relevant TE

In this section I go in details analyzing the **copynumber geographical
distribution** for the TE families previously selected in the subsets.

First, I analysed the difference between the two sex-related subsets. I
found that the **females** dataset (`fTEoutlier_order`) has 40 different
`familyname`, while the **males** has 44. While there are no sequence
families present in the females and not in the males, there are 4
sequence families present in the males and not in the females:

- `L1MC1`
- `L1PB2`
- `L1PB4`
- `MER22` (a satellite)

In other words, these sequences have high variance particularly on the
**Y chromosome**. I analysed those separately.

``` r
f_familynames <- f_outliers$familyname
m_familynames <- m_outliers$familyname

# TEs present in the female subset and not in the males
length(unique(f_familynames))
```

    ## [1] 40

``` r
setdiff(f_familynames,m_familynames)
```

    ## character(0)

``` r
# TEs present in the male subset and not in the females
length(unique(m_familynames))
```

    ## [1] 44

``` r
setdiff(m_familynames,f_familynames)
```

    ## [1] "L1MC1" "L1PB2" "L1PB4" "MER22"

I wrote this function to easily plot each TE family. The 8 arguments
are:

- `data` = Dataset where are the data to be plotted
- `sex` = provide ‘Males’ or ‘Females’
- `famname`= provide the TE `familyname`
- `binwidht` = widht of the bin of the histogram

All the other arguments to provide are ‘y’ for yes and ‘n’ for no,
depending if you want the following plot features.

- `x_title`= Title of the x axis
- `y_title`= Title of the y axis
- `x_numbers`= Ticks and scale of the x axis
- `y_numbers`= Ticks and scale of the y axis

``` r
plotTEfamily <- function(data, sex, famname, binwidht, x_title, y_title, x_numbers, y_numbers){
filtered <- filter(data, familyname==famname)
ggplot(data = filtered, mapping = aes(x = copynumber, fill = Country)) +
  geom_histogram(binwidth = binwidht) + 
  ggtitle(paste0(sex,' - ',famname)) + theme(plot.title = element_text(size = 8, hjust = 0.5)) +
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

### L1 family

To analyse all the TE families part of the **L1 superfamily**, I first
create two subsubsets containing the useful data. Only the observation
with `familyname` starting with `L1` were retained.

``` r
fL1<-f_outliers %>% filter(str_detect(f_outliers$familyname, "^L1"))
mL1<-m_outliers %>% filter(str_detect(m_outliers$familyname, "^L1"))
```

Subsequently, I created the plots for each L1 family contained in the
data. In total, I had 16 L1 families, for a total of 32 plots (16 for
males, 16 for females). Here is the list of the investigated L1
families, representing the most variant families in absolute terms:

- `L1PA4`
- `L1`
- `L1PREC1`
- `L1PA16`
- `L1PA6`
- `L1PA7_5`
- `L1HS`
- `L1PA7`
- `L1PB2c`
- `L1PB1`
- `L1PA10`
- `L1PA8`
- `L1PREC2`
- `L1PA3`
- `L1PA15`
- `L1P_MA2`

``` r
# L1PA4
mL1PA4_plot<-plotTEfamily(m_outliers, 'Males', 'L1PA4', 15, 'n', 'n', 'n', 'n')
fL1PA4_plot<-plotTEfamily(f_outliers, 'Females', 'L1PA4', 15, 'n', 'n', 'n', 'n')

# L1
mL1_plot<-plotTEfamily(m_outliers, 'Males', 'L1', 15, 'n', 'n', 'n', 'n')
fL1_plot<-plotTEfamily(f_outliers, 'Females', 'L1', 15, 'n', 'n', 'n', 'n')

# L1PREC1
mL1PREC1_plot<-plotTEfamily(m_outliers, 'Males', 'L1PREC1', 15, 'n', 'n', 'n', 'n')
fL1PREC1_plot<-plotTEfamily(f_outliers, 'Females', 'L1PREC1', 15, 'n', 'n', 'n', 'n')

# L1PA16
mL1PA16_plot<-plotTEfamily(m_outliers, 'Males', 'L1PA16', 15, 'n', 'n', 'n', 'n')
fL1PA16_plot<-plotTEfamily(f_outliers, 'Females', 'L1PA16', 15, 'n', 'n', 'n', 'n')

# L1PA6
mL1PA6_plot<-plotTEfamily(m_outliers, 'Males', 'L1PA6', 15, 'n', 'n', 'n', 'n')
fL1PA6_plot<-plotTEfamily(f_outliers, 'Females', 'L1PA6', 15, 'n', 'n', 'n', 'n')

# L1PA7_5
mL1PA7_5_plot<-plotTEfamily(m_outliers, 'Males', 'L1PA7_5', 15, 'n', 'n', 'n', 'n')
fL1PA7_5_plot<-plotTEfamily(f_outliers, 'Females', 'L1PA7_5', 15, 'n', 'n', 'n', 'n')

# L1HS
mL1HS_plot<-plotTEfamily(m_outliers, 'Males', 'L1HS', 15, 'n', 'n', 'n', 'n')
fL1HS_plot<-plotTEfamily(f_outliers, 'Females', 'L1HS', 15, 'n', 'n', 'n', 'n')

# L1PA7
mL1PA7_plot<-plotTEfamily(m_outliers, 'Males', 'L1PA7', 15, 'n', 'n', 'n', 'n')
fL1PA7_plot<-plotTEfamily(f_outliers, 'Females', 'L1PA7', 15, 'n', 'n', 'n', 'n')

# L1PB2c
mL1PB2c_plot<-plotTEfamily(m_outliers, 'Males', 'L1PB2c', 15, 'n', 'n', 'n', 'n')
fL1PB2c_plot<-plotTEfamily(f_outliers, 'Females', 'L1PB2c', 15, 'n', 'n', 'n', 'n')

# L1PB1
mL1PB1_plot<-plotTEfamily(m_outliers, 'Males', 'L1PB1', 15, 'n', 'n', 'n', 'n')
fL1PB1_plot<-plotTEfamily(f_outliers, 'Females', 'L1PB1', 15, 'n', 'n', 'n', 'n')

# L1PA10
mL1PA10_plot<-plotTEfamily(m_outliers, 'Males', 'L1PA10', 15, 'n', 'n', 'n', 'n')
fL1PA10_plot<-plotTEfamily(f_outliers, 'Females', 'L1PA10', 15, 'n', 'n', 'n', 'n')

# L1PA8
mL1PA8_plot<-plotTEfamily(m_outliers, 'Males', 'L1PA8', 15, 'n', 'n', 'n', 'n')
fL1PA8_plot<-plotTEfamily(f_outliers, 'Females', 'L1PA8', 15, 'n', 'n', 'n', 'n')

# L1PREC2
mL1PREC2_plot<-plotTEfamily(m_outliers, 'Males', 'L1PREC2', 15, 'n', 'n', 'n', 'n')
fL1PREC2_plot<-plotTEfamily(f_outliers, 'Females', 'L1PREC2', 15, 'n', 'n', 'n', 'n')

# L1PA3
mL1PA3_plot<-plotTEfamily(m_outliers, 'Males', 'L1PA3', 15, 'n', 'n', 'n', 'n')
fL1PA3_plot<-plotTEfamily(f_outliers, 'Females', 'L1PA3', 15, 'n', 'n', 'n', 'n')

# L1PA15
mL1PA15_plot<-plotTEfamily(m_outliers, 'Males', 'L1PA15', 15, 'n', 'n', 'n', 'n')
fL1PA15_plot<-plotTEfamily(f_outliers, 'Females', 'L1PA15', 15, 'n', 'n', 'n', 'n')

# L1P_MA2
mL1P_MA2_plot<-plotTEfamily(m_outliers, 'Males', 'L1P_MA2', 15, 'n', 'n', 'n', 'n')
fL1P_MA2_plot<-plotTEfamily(f_outliers, 'Females', 'L1P_MA2', 15, 'n', 'n', 'n', 'n')
```

Eventually, I create the two final figures, one for each sex. They
contain 16 plots each, one for each L1 family present in the analysis.

``` r
m_L1_figure <- ggarrange(mL1PA4_plot, mL1_plot, mL1PREC1_plot, mL1PA16_plot, mL1PA6_plot, mL1PA7_5_plot, mL1HS_plot, mL1PA7_plot, mL1PB2c_plot, mL1PB1_plot, mL1PA10_plot, mL1PA8_plot, mL1PREC2_plot, mL1PA3_plot, mL1PA15_plot, mL1P_MA2_plot, ncol = 4, nrow = 4, common.legend = TRUE, legend = "top", align = "hv", font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

m_L1_final <- annotate_figure(m_L1_figure, left = text_grob("Count", color = "black", rot = 90), bottom = text_grob("Copynumber", color = "black"), top = text_grob("Males", color = "black", size = 20), fig.lab = "")


f_L1_figure <- ggarrange(fL1PA4_plot, fL1_plot, fL1PREC1_plot, fL1PA16_plot, fL1PA6_plot, fL1PA7_5_plot, fL1HS_plot, fL1PA7_plot, fL1PB2c_plot, fL1PB1_plot, fL1PA10_plot, fL1PA8_plot, fL1PREC2_plot, fL1PA3_plot, fL1PA15_plot, fL1P_MA2_plot, ncol = 4, nrow = 4, common.legend = TRUE, legend = "top", align = "hv", font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f_L1_final <- annotate_figure(f_L1_figure, left = text_grob("Count", color = "black", rot = 90), bottom = text_grob("Copynumber", color = "black"), top = text_grob("Females", color = "black", size = 20), fig.lab = "")

f_L1_final
```

![](1_HGDP_Absolute-TEvariation-copy-optimized_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
m_L1_final
```

![](1_HGDP_Absolute-TEvariation-copy-optimized_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

What do we understand from these plots? I try to summarize my ideas:

- Many L1 distributions (10/16) are closer to a **normal** distribution
  (`L1`, `L1PREC`, `L1PA16`, `L1PA6`, `L1HS`, `L1PB2c`, `L1PA10`,
  `L1PREC2`, `L1PA15`, `L1P_MA2`).
- All the other families (`L1PA4`, `L1PA7_5`, `L1PA7`, `L1PB1`, `L1PA8`,
  `L1PA3`) shows a distribution closer to a **bimodal** (6/16). Some are
  more pronounced, some are less evident. These distributions are more
  evident in the `male` figure (remember that the number of observations
  in the male plot is much higher, and thus more significant).
- In all the bimodal distributed families, the **africans** are strongly
  biased for the first peak. There is a strong exception, which is
  `L1PA7_5`, which shows almost all african observation in the second
  peak (1/6).
- Even in the normal distributed families, we can notice that
  **africans** are almost never in the second part of the distribution,
  with the exceptions of `L1` and `L1HS` (2/10).
- The outliers with high copynumber are almost always **europeans**.

My final interpretation of this figures is that the **bimodal**
distribution of some families can be explained with a **TE invasion**.
The first peak is composed by individual with low copynumber, thus the
TE invasion did not spread into these populations. The second peak
contains individuals with high copynumber, indicating a spread of the TE
in these populations. These invasions are particularly clear for
`L1PA3`, `L1PA7_5`, `L1PB1` and `L1PA7_5`.

The fact that the **africans** are almost always on the left side of the
distributions is in line with my previous observation (see script 1):
this can be explained by the fact that during the OOA bottleneck
multiple TE invasions were triggered. For this reason, I think it can be
interesting to look more into details also the L1 families that shows
different distribution of africans (`L1PA7_5`, `L1`, `L1HS`).

We also see a **sex-specific** pattern, very clear in the L1s with
bimodal distribution: for the females, almost all the low copynumber
individuals are **africans**, while the low copynumber males are more
mixed, despite the african males being all in the low copynumber part of
the distribution. This pattern is present also in the `L1PA7_5`
distribution, but reversed: the africans are in the second peak (high
copynumber).

#### L1 sequences present only in the males subset

There are 3 L1 families which are only present in the male subset:
`L1MC1`, `L1PB2`, `L1PB4`. Thus, they not match the condition of
“absolute difference \> 200” to be included in the female dataset. In
order to compare male and female distributions of these sequences, I
take the female data directly from the main dataset, not filtered.

``` r
# L1MC1
mL1MC1_plot<-plotTEfamily(m_outliers, 'Males', 'L1MC1', 10, 'n', 'n', 'n', 'n')
fL1MC1_plot<-plotTEfamily(females, 'Females', 'L1MC1', 10, 'n', 'n', 'n', 'n')

# L1PB2
mL1PB2_plot<-plotTEfamily(m_outliers, 'Males', 'L1PB2', 10, 'n', 'n', 'n', 'n')
fL1PB2_plot<-plotTEfamily(females, 'Females', 'L1PB2', 10, 'n', 'n', 'n', 'n')

# L1PB4
mL1PB4_plot<-plotTEfamily(m_outliers, 'Males', 'L1PB4', 10, 'n', 'n', 'n', 'n')
fL1PB4_plot<-plotTEfamily(females, 'Females', 'L1PB4', 10, 'n', 'n', 'n', 'n')
```

``` r
m_L1_Y_figure <- ggarrange(mL1MC1_plot, mL1PB2_plot, mL1PB4_plot, ncol = 1, nrow = 3, common.legend = TRUE, legend = "right", align = "hv", font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

m_L1_Y_final <- annotate_figure(m_L1_Y_figure, left = text_grob("Count", color = "black", rot = 90), bottom = text_grob("Copynumber", color = "black"), top = text_grob("Males", color = "black", size = 20), fig.lab = "")


f_L1_Y_figure <- ggarrange(fL1MC1_plot, fL1PB2_plot, fL1PB4_plot, ncol = 1, nrow = 3, common.legend = TRUE, legend = "right", align = "hv", font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f_L1_Y_final <- annotate_figure(f_L1_Y_figure, left = text_grob("Count", color = "black", rot = 90), bottom = text_grob("Copynumber", color = "black"), top = text_grob("Females", color = "black", size = 20), fig.lab = "")

f_L1_Y_final
```

![](1_HGDP_Absolute-TEvariation-copy-optimized_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
m_L1_Y_final
```

![](1_HGDP_Absolute-TEvariation-copy-optimized_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

- First, I notice that, even if these sequences were not included in the
  female subset, their absolute variance is very close to 200.
- I also notice 3/3 distributions that looks intermediate bewteen a
  **normal** and a **bimodal**.
- The **geographical pattern** is again very clear, with **africans** on
  the left side of the plots and **europeans** that usually represent
  the high copynumber outliers (sometimes with **middle east**
  individuals).

### Non-LTR Retrotransposons

#### ALU

The most variant sequences found in the dataset is, not suprisingly,
`ALU`. This is a SINE TE, the most abundant in the human genome.

``` r
(mALU_plot<-plotTEfamily(males, 'Males', 'ALU', 750, 'y', 'y', 'y', 'y'))
```

![](1_HGDP_Absolute-TEvariation-copy-optimized_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
(fALU_plot<-plotTEfamily(females, 'Females', 'ALU', 1000, 'y', 'y', 'y', 'y'))
```

![](1_HGDP_Absolute-TEvariation-copy-optimized_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->

Also for `ALU`, we see the same **bimodal** pattern already present in
some L1 families. We also see a **sex-specific** pattern: for the
females, almost all the low copynumber individuals are **africans**,
while the low copynumber males are more mixed, despite the african males
being all in the low copynumber part of the distribution.

#### SVA_A

Only one element of the SVA superfamily appears among the TEs with
highest absolute variance: `SVA_A`.

``` r
(mSVA_A_plot<-plotTEfamily(males, 'Males', 'SVA_A', 8, 'y', 'y', 'y', 'y'))
```

![](1_HGDP_Absolute-TEvariation-copy-optimized_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
(fSVA_A_plot<-plotTEfamily(females, 'Females', 'SVA_A', 10, 'y', 'y', 'y', 'y'))
```

![](1_HGDP_Absolute-TEvariation-copy-optimized_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->

Also here we see a **bimodal** distribution with **africans** on the
left side of the distribution. I also notice the same pattern
highlighted in the script 1 for L1ME1: while the **females** with low
copynumber are almost all africans, the **males** with low copynumber
are mixed, but the great part of male africans are still showing low
copynumbers (left side).

### DNA transposon

The only **DNA transposons** present in the dataset are `MER2` and
`TIGGER1`.

``` r
(mMER2_plot<-plotTEfamily(males, 'Males', 'MER2', 8, 'y', 'y', 'y', 'y'))
```

![](1_HGDP_Absolute-TEvariation-copy-optimized_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
(fMER2_plot<-plotTEfamily(females, 'Females', 'MER2', 8, 'y', 'y', 'y', 'y'))
```

![](1_HGDP_Absolute-TEvariation-copy-optimized_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->

Again, **bimodal** distribution with the same **africans-OOU** and
**males-females** pattern.

``` r
(mTIGGER1_plot<-plotTEfamily(males, 'Males', 'TIGGER1', 10, 'y', 'y', 'y', 'y'))
```

![](1_HGDP_Absolute-TEvariation-copy-optimized_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
(fTIGGER1_plot<-plotTEfamily(females, 'Females', 'TIGGER1', 10, 'y', 'y', 'y', 'y'))
```

![](1_HGDP_Absolute-TEvariation-copy-optimized_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->

I do not notice anything interesting in this pattern. It makes sense to
me because we expect DNA transposons to be much more stable than
retrotransposon during recent time. It could be an “ancient” transposon
with no more activity?

### Endogenous retroviruses

`THE1B`, `THE1_I`, `THE1C`, `MSTA`, `THE1A`, `THE1D`, `MLT2A2`, `MLT1B`,
`MLT2A1` are all **ERV3**.

``` r
# THE1B
mTHE1B_plot<-plotTEfamily(males, 'Males', 'THE1B', 25, 'n', 'n', 'n', 'n')
fTHE1B_plot<-plotTEfamily(females, 'Females', 'THE1B', 25, 'n', 'n', 'n', 'n')

# THE1_I
mTHE1_I_plot<-plotTEfamily(males, 'Males', 'THE1_I', 25, 'n', 'n', 'n', 'n')
fTHE1_I_plot<-plotTEfamily(females, 'Females', 'THE1_I', 25, 'n', 'n', 'n', 'n')

# THE1C
mTHE1C_plot<-plotTEfamily(males, 'Males', 'THE1C', 25, 'n', 'n', 'n', 'n')
fTHE1C_plot<-plotTEfamily(females, 'Females', 'THE1C', 25, 'n', 'n', 'n', 'n')

# MSTA
mMSTA_plot<-plotTEfamily(males, 'Males', 'MSTA', 25, 'n', 'n', 'n', 'n')
fMSTA_plot<-plotTEfamily(females, 'Females', 'MSTA', 25, 'n', 'n', 'n', 'n')

# THE1A
mTHE1A_plot<-plotTEfamily(males, 'Males', 'THE1A', 25, 'n', 'n', 'n', 'n')
fTHE1A_plot<-plotTEfamily(females, 'Females', 'THE1A', 25, 'n', 'n', 'n', 'n')

#THE1D
mTHE1D_plot<-plotTEfamily(males, 'Males', 'THE1D', 25, 'n', 'n', 'n', 'n')
fTHE1D_plot<-plotTEfamily(females, 'Females', 'THE1D', 25, 'n', 'n', 'n', 'n')

# MLT2A2
mMLT2A2_plot<-plotTEfamily(males, 'Males', 'MLT2A2', 10, 'n', 'n', 'n', 'n')
fMLT2A2_plot<-plotTEfamily(females, 'Females', 'MLT2A2', 10, 'n', 'n', 'n', 'n')

# MLT1B
mMLT1B_plot<-plotTEfamily(males, 'Males', 'MLT1B', 10, 'n', 'n', 'n', 'n')
fMLT1B_plot<-plotTEfamily(females, 'Females', 'MLT1B', 10, 'n', 'n', 'n', 'n')

# MLT2A1
mMLT2A1_plot<-plotTEfamily(m_outliers, 'Males', 'MLT2A1', 10, 'n', 'n', 'n', 'n')
fMLT2A1_plot<-plotTEfamily(females, 'Females', 'MLT2A1', 10, 'n', 'n', 'n', 'n')
```

``` r
m_ERV3_figure <- ggarrange(mTHE1B_plot, mTHE1_I_plot, mTHE1C_plot, mMSTA_plot, mTHE1A_plot, mTHE1D_plot, mMLT2A2_plot, mMLT1B_plot, mMLT2A1_plot, ncol = 3, nrow = 3, common.legend = TRUE, legend = "top", align = "hv", font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

m_ERV3_final <- annotate_figure(m_ERV3_figure, left = text_grob("Count", color = "black", rot = 90), bottom = text_grob("Copynumber", color = "black"), top = text_grob("Males", color = "black", size = 20), fig.lab = "")


f_ERV3_figure <- ggarrange(fTHE1B_plot, fTHE1_I_plot, fTHE1C_plot, fMSTA_plot, fTHE1A_plot, fTHE1D_plot, fMLT2A2_plot, fMLT1B_plot, fMLT2A1_plot, ncol = 3, nrow = 3, common.legend = TRUE, legend = "top", align = "hv", font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f_ERV3_final <- annotate_figure(f_ERV3_figure, left = text_grob("Count", color = "black", rot = 90), bottom = text_grob("Copynumber", color = "black"), top = text_grob("Females", color = "black", size = 20), fig.lab = "")

f_ERV3_final
```

![](1_HGDP_Absolute-TEvariation-copy-optimized_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
m_ERV3_final
```

![](1_HGDP_Absolute-TEvariation-copy-optimized_files/figure-gfm/unnamed-chunk-17-2.png)<!-- -->

We see 8/9 distributions close to a **normal**, and 1/9 **bimodal**
(`MLT2A1`). In the `MLT2A1` distribution we do notice the usual pattern
previously described (africans with low copynumber, males with low
copynumber more mixed than the females).

### Satellites

All the sequences present in the dataset with the most variant sequences
that were not analysed yet are **satellites**:`ALR1`, `ALR_`, `HSATII`,
`ALR`, `ALRb`, `HSATI`, `ALR2`, `ALRa_`, `CER`, `LSAU`, `6kbHsap`. I’m
not sure about the significance and importance of this results, but I
performed the analysis anyway.

I also added the only satellite present only in the **male** subset and
not in the female one, which is `MER22`.

``` r
# ALR1
mALR1_plot<-plotTEfamily(males, 'Males', 'ALR1', 1000, 'n', 'n', 'n', 'n')
fALR1_plot<-plotTEfamily(females, 'Females', 'ALR1', 1000, 'n', 'n', 'n', 'n')

# ALR_
mALR__plot<-plotTEfamily(males, 'Males', 'ALR_', 1000, 'n', 'n', 'n', 'n')
fALR__plot<-plotTEfamily(females, 'Females', 'ALR_', 1000, 'n', 'n', 'n', 'n')

# HSATII
mHSATII_plot<-plotTEfamily(males, 'Males', 'HSATII', 1000, 'n', 'n', 'n', 'n')
fHSATII_plot<-plotTEfamily(females, 'Females', 'HSATII', 1000, 'n', 'n', 'n', 'n')

# ALR
mALR_plot<-plotTEfamily(males, 'Males', 'ALR', 1000, 'n', 'n', 'n', 'n')
fALR_plot<-plotTEfamily(females, 'Females', 'ALR', 1000, 'n', 'n', 'n', 'n')

# ALRb
mALRb_plot<-plotTEfamily(males, 'Males', 'ALRb', 1000, 'n', 'n', 'n', 'n')
fALRb_plot<-plotTEfamily(females, 'Females', 'ALRb', 1000, 'n', 'n', 'n', 'n')

# HSATI
mHSATI_plot<-plotTEfamily(males, 'Males', 'HSATI', 100, 'n', 'n', 'n', 'n')
fHSATI_plot<-plotTEfamily(females, 'Females', 'HSATI', 100, 'n', 'n', 'n', 'n')

# ALR2
mALR2_plot<-plotTEfamily(males, 'Males', 'ALR2', 50, 'n', 'n', 'n', 'n')
fALR2_plot<-plotTEfamily(females, 'Females', 'ALR2', 50, 'n', 'n', 'n', 'n')

# ALRa_
mALRa__plot<-plotTEfamily(males, 'Males', 'ALRa_', 25, 'n', 'n', 'n', 'n')
fALRa__plot<-plotTEfamily(females, 'Females', 'ALRa_', 25, 'n', 'n', 'n', 'n')

# CER
mCER_plot<-plotTEfamily(males, 'Males', 'CER', 10, 'n', 'n', 'n', 'n')
fCER_plot<-plotTEfamily(females, 'Females', 'CER', 10, 'n', 'n', 'n', 'n')

# LSAU
mLSAU_plot<-plotTEfamily(males, 'Males', 'LSAU', 10, 'n', 'n', 'n', 'n')
fLSAU_plot<-plotTEfamily(females, 'Females', 'LSAU', 10, 'n', 'n', 'n', 'n')

# 6kbHsap
m6kbHsap_plot<-plotTEfamily(males, 'Males', '6kbHsap', 20, 'n', 'n', 'n', 'n')
f6kbHsap_plot<-plotTEfamily(females, 'Females', '6kbHsap', 20, 'n', 'n', 'n', 'n')

# MER22
mMER22_plot<-plotTEfamily(males, 'Males', 'MER22', 10, 'n', 'n', 'n', 'n')
fMER22_plot<-plotTEfamily(females, 'Females', 'MER22', 10, 'n', 'n', 'n', 'n')
```

``` r
m_SAT_figure <- ggarrange(mALR1_plot, mALR__plot, mHSATII_plot, mALR_plot, mALRb_plot, mHSATI_plot, mALR2_plot, mALRa__plot, mCER_plot, mLSAU_plot, m6kbHsap_plot, mMER22_plot, ncol = 3, nrow = 4, common.legend = TRUE, legend = "top", align = "hv", font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

m_SAT_final <- annotate_figure(m_SAT_figure, left = text_grob("Count", color = "black", rot = 90), bottom = text_grob("Copynumber", color = "black"), top = text_grob("Males", color = "black", size = 20), fig.lab = "")


f_SAT_figure <- ggarrange(fALR1_plot, fALR__plot, fHSATII_plot, fALR_plot, fALRb_plot, fHSATI_plot, fALR2_plot, fALRa__plot, fCER_plot, fLSAU_plot, f6kbHsap_plot, fMER22_plot, ncol = 3, nrow = 4, common.legend = TRUE, legend = "top", align = "hv", font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

f_SAT_final <- annotate_figure(f_SAT_figure, left = text_grob("Count", color = "black", rot = 90), bottom = text_grob("Copynumber", color = "black"), top = text_grob("Females", color = "black", size = 20), fig.lab = "")

f_SAT_final
```

![](1_HGDP_Absolute-TEvariation-copy-optimized_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
m_SAT_final
```

![](1_HGDP_Absolute-TEvariation-copy-optimized_files/figure-gfm/unnamed-chunk-19-2.png)<!-- -->

I notice only **normal** distributions. Is there something interesting
here?
