# ramen

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
    misosoup24[[1]], name = names(misosoup24)[1]
)
cm2 <- ConsortiumMetabolism(
    misosoup24[[2]], name = names(misosoup24)[2]
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
groups, and visualization – see
[`vignette("ramen", package = "ramen")`](https://admarhi.github.io/ramen/articles/ramen.md).

## Compare many communities at once

Bundle a handful of consortia into a `ConsortiumMetabolismSet`, align
them all-against-all, and plot the resulting similarity matrix:

``` r

sel <- names(misosoup24)[seq_len(8)]
cms <- do.call(
    ConsortiumMetabolismSet,
    c(
        lapply(sel, function(n) {
            ConsortiumMetabolism(misosoup24[[n]], name = n)
        }),
        list(name = "misosoup24[1:8]", verbose = FALSE)
    )
)
plot(align(cms))
```

![Heatmap of pairwise FOS scores across eight MiSoSoup consortia, with
rows and columns reordered by hierarchical
clustering.](reference/figures/README-hero-1.png)

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
`colData()`, [`dim()`](https://rdrr.io/r/base/dim.html), etc.) work out
of the box.

## Vignettes

- [`vignette("ramen")`](https://admarhi.github.io/ramen/articles/ramen.md)
  – full introduction and analysis pipeline
- [`vignette("alignment")`](https://admarhi.github.io/ramen/articles/alignment.md)
  – alignment metrics, coverage ratios, and permutation p-values
- [`vignette("visualisation")`](https://admarhi.github.io/ramen/articles/visualisation.md)
  – gallery of all plot types
