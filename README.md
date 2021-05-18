
# enet

<!-- badges: start -->
[![R-CMD-check](https://github.com/yufree/enet/workflows/R-CMD-check/badge.svg)](https://github.com/yufree/enet/actions)
<!-- badges: end -->

The goal of enet is to perform network analysis between exposome and metabolites. It also include function to perform network analysis among metabolites and exposures.

## Installation

You can install the released version of enet from GitHub with:

``` r
remotes::install_github("yufree/enet")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(enet)
data(expo)
date(meta)
re <- getgk(exp(meta),expo)
```

