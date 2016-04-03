## Simulation code, simulated data sets, and analysis for *A Genealogical Look at Shared Ancestry on the X Chromosome*

This repository contains the code and other materials for our paper **A
Genealogical Look at Shared Ancestry on the X Chromosome*. A preprint of our
manuscript is [available on
bioRxiv](http://biorxiv.org/content/early/2016/04/03/046912).


## Contents

 - `Makefile`: main makefile for the repository; runs simulations to generate data, runs analysis for plots.

 - `data/`: directory containing simulated datasets used in the paper.

 - `ancestrysim/`: directory containing genealogical simulation code. This
   includes the prototype Python routine `ancestralsim.py` and the finalized
   C/Python `ancestrysim2.py` version.

  - `paper/`: empty directory which was home to the manuscript, now contains
    `paper/images` which contains images produced by `plots.R`.

  - `plots.R`: R code to run analyses and generate paper's figures.

  - `dist-functions.R`: distribution functions, test code, and other functions used by `plots.R`.

  - `jsviz/`: code for visualizations made with Javascript/d3 for paper.

  - `images/`: other images and image requisites like `.dot` files.

## Session Info

Here's the `sessionInfo()` from R, after `plots.R` was run:

```R
R version 3.2.2 (2015-08-14)
Platform: x86_64-apple-darwin14.5.0 (64-bit)
Running under: OS X 10.10.5 (Yosemite)

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods
[8] base

other attached packages:
 [1] colorspace_1.2-6   RColorBrewer_1.1-2 readr_0.2.2        tidyr_0.4.0
 [5] wesanderson_0.3.2  reshape2_1.4.1     purrr_0.2.0        dplyr_0.4.3
 [9] scales_0.3.0       ggplot2_2.0.0      devtools_1.9.1

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.2     digest_0.6.9    assertthat_0.1  plyr_1.8.3
 [5] R6_2.1.1        gtable_0.1.2    DBI_0.3.1       magrittr_1.5
 [9] stringi_1.0-1   lazyeval_0.1.10 tools_3.2.2     stringr_1.0.0
[13] munsell_0.4.2   parallel_3.2.2  memoise_0.2.1k
```
