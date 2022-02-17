
# enet

<!-- badges: start -->
[![R-CMD-check](https://github.com/yufree/enet/workflows/R-CMD-check/badge.svg)](https://github.com/yufree/enet/actions)
<!-- badges: end -->

The goal of enet is to perform network analysis between exposome and metabolites. It also include function to perform network analysis among metabolites and exposures. At current stage, this package can be used to find gatekeepers.

## Installation

You can install the released version of enet from GitHub with:

``` r
remotes::install_github("yufree/enet")
```

## Example

This is a basic example to find gatekeepers

``` r
library(enet)
data(expo)
data(meta)
re <- getgk(exp(meta),expo)
```

## Citation

- Yu, M.; Teitelbaum, S. L.; Dolios, G.; Dang, L.-H. T.; Tu, P.; Wolff, M. S.; Petrick, L. M. Molecular Gatekeeper Discovery: Workflow for Linking Multiple Exposure Biomarkers to Metabolomics. Environ. Sci. Technol. 2022. https://doi.org/10.1021/acs.est.1c04039.

