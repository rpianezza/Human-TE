HGDP - Summary
================

The previous 8 scripts highlight some TE which are more interesting in
our analysis, but exclude some of them after different controls. Here I
want to answer a simple but crucial question for the future paper: which
are the most interesting TEs found in the HGDP project? These TEs
deserve a more detailed analysis (**phylogenetic trees**).

In the **scripts 1-2-3** I looked for the most variable TEs copynumbers
across all the genomes, both in absolute and relative terms. In the same
analysis, I plot the copynumber distribution of each selected TE,
dividing the samples for geographical origin. Here is the TEs which
remained interesting after the plotting:

- **LINE-1 family**: lot of members of this family showed high variance
  as well as geography-specific pattern of copynumber distribution. Here
  is the list of the selected L1s: `L1PA4`, `L1`, `L1PREC1`, `L1PA16`,
  `L1PA6`, `L1PA7_5`, `L1PB2c`, `L1PB1`, `L1PA10`, `L1PREC2`, `L1PA15`,
  `L1P_MA2`, `L1MC1`, `L1PB2`, `L1PB4`.
  - I removed `L1ME5` for the extremely **uneven coverage** as shown in
    **script 5**.
  - I removed `L1PA3`, `L1PA7`, `L1PA8` and `L1HS` which are showing
    some percentage of **cross-mapped reads**, like shown in **script
    8**. Note that the level of cross-mapping are not so high to
    invalidate our results, but since the goal of this script is to
    select just few TEs, this is one of the “filtering” steps.
  - From the remaining TEs, we can select:
    - `L1PB1`: the L1 with the more pronounced bimodal distribution,
      showing the classic pattern of copynumber distribution (low in
      Africans, high in Eurasia).
    - `L1PA7_5`: bimodal distribution but reverse pattern of copynumber
      distribution (low in Eurasia, high in Africans).
- **SINEs**: `ALU` and `SVA_A`, both showing the classic patterns of
  copynumber distributions and both known to be part of some of the most
  active TEs in humans.
- `MER2`: a **DNA TE** showing the classic patterns.
- `MLT2A1`: an **endogenous retrovirus** showing the classic patterns.
- As a **control**, I would like to add to the list `L2`, known to be
  ancient, inactive and which shows low variance among individuals in
  our dataset.

Summarizing, the final list of TEs to be investigated in more details
is:

- `L1PB1`
- `L1PA7_5`
- `ALU`
- `SVA_A`
- `MER2`
- `MLT2A1`
- `L2` (control)
