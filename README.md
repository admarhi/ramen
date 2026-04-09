
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ramen

<!-- badges: start -->

[![R-CMD-check](https://github.com/admarhi/ramen/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/admarhi/ramen/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/admarhi/ramen/branch/main/graph/badge.svg)](https://codecov.io/gh/admarhi/ramen?branch=main)
<!-- badges: end -->

**ramen** (**R**econstruction and **A**lignment of **M**icrobial
**E**xchange **N**etworks) compares microbial communities by their
metabolic exchange networks rather than by species composition. Given
species-metabolite-flux data (e.g. from community flux balance
analysis), ramen reconstructs directed metabolite-to-metabolite networks
and quantifies how functionally similar different communities are – even
when they share no species at all.

## Quick start

``` r
library(ramen)
data("misosoup24")

## Build two consortium objects
cm1 <- ConsortiumMetabolism(
    misosoup24[[1]], name = names(misosoup24)[1],
    species_col = "species",
    metabolite_col = "metabolites",
    flux_col = "fluxes"
)
cm2 <- ConsortiumMetabolism(
    misosoup24[[2]], name = names(misosoup24)[2],
    species_col = "species",
    metabolite_col = "metabolites",
    flux_col = "fluxes"
)

## Align them
cma <- align(cm1, cm2)
scores(cma)
#> $FOS
#> [1] 0.7634409
#> 
#> $jaccard
#> [1] 0.3397129
#> 
#> $brayCurtis
#> [1] 0.3115158
#> 
#> $redundancyOverlap
#> [1] 0.3397129
#> 
#> $coverageQuery
#> [1] 0.7634409
#> 
#> $coverageReference
#> [1] 0.3796791
```

For the full workflow – from data import through alignment, functional
groups, and visualization – see `vignette("ramen", package = "ramen")`.

## Installation

Install the development version from GitHub:

``` r
pak::pkg_install("admarhi/ramen")
# or
devtools::install_github("admarhi/ramen")
```

## Core workflow

1.  **Import** species-metabolite-flux data (from MiSoSoup YAML, CSV, or
    any data.frame)
2.  **Build** `ConsortiumMetabolism` objects (one per community)
3.  **Group** them into a `ConsortiumMetabolismSet` (computes pairwise
    overlap and clustering)
4.  **Align** communities pairwise or across the full set using five
    complementary metrics (FOS, Jaccard, Bray-Curtis, Redundancy
    Overlap, MAAS)
5.  **Visualise** with heatmaps, network plots, dendrograms, and
    functional group analyses

All three classes inherit from Bioconductor’s
`TreeSummarizedExperiment`, so standard accessor functions (`assay()`,
`colData()`, `dim()`, etc.) work out of the box.

## Vignettes

- `vignette("ramen")` – full introduction and analysis pipeline
- `vignette("alignment")` – alignment metrics, coverage ratios, and
  permutation p-values
- `vignette("visualisation")` – gallery of all plot types
