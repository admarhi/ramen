# Visualisation

## Introduction

`ramen` provides
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods for all
three core classes: `ConsortiumMetabolism`, `ConsortiumMetabolismSet`,
and `ConsortiumMetabolismAlignment`. This vignette demonstrates each
plot type with synthetic and example data.

``` r
library(ramen)
```

## Preparing example data

We create four synthetic consortia, combine them into a CMS, and compute
both pairwise and multiple alignments.

``` r
cm_alpha <- synCM("Alpha", n_species = 5, max_met = 10, seed = 42)
cm_beta <- synCM("Beta", n_species = 5, max_met = 10, seed = 43)
cm_gamma <- synCM("Gamma", n_species = 4, max_met = 8, seed = 44)
cm_delta <- synCM("Delta", n_species = 6, max_met = 10, seed = 45)

cms <- ConsortiumMetabolismSet(
    cm_alpha, cm_beta, cm_gamma, cm_delta,
    name = "Demo"
)

cma_pair <- align(cm_alpha, cm_beta)
cma_mult <- align(cms)
```

## ConsortiumMetabolism plots

`plot(CM)` renders a directed metabolic flow network using `igraph`.
Nodes are coloured by role: **lightblue** = source (only outgoing),
**salmon** = sink (only incoming), **yellow** = intermediate (both). The
`type` argument selects which assay matrix to use for edge weights.

### Binary network

The default view shows the structure of the metabolic network where any
edge present has weight 1.

``` r
plot(cm_alpha, type = "Binary")
```

![Binary metabolic network for
Alpha.](visualisation_files/figure-html/plot-cm-binary-1.png)

Binary metabolic network for Alpha.

### Number of species per pathway

Edge weights reflect how many species catalyse each pathway.

``` r
plot(cm_alpha, type = "nSpecies")
```

![nSpecies network for
Alpha.](visualisation_files/figure-html/plot-cm-nspecies-1.png)

nSpecies network for Alpha.

### Consumption and production

Consumption and production views weight edges by summed flux magnitudes.

``` r
plot(cm_alpha, type = "Consumption")
```

![Consumption-weighted
network.](visualisation_files/figure-html/plot-cm-consumption-1.png)

Consumption-weighted network.

``` r
plot(cm_alpha, type = "Production")
```

![Production-weighted
network.](visualisation_files/figure-html/plot-cm-production-1.png)

Production-weighted network.

### Effective consumption and production

The effective diversity measures indicate how evenly species contribute
to each pathway (Shannon-based effective number of species).

``` r
plot(cm_alpha, type = "EffectiveConsumption")
```

![Effective consumption
network.](visualisation_files/figure-html/plot-cm-effcons-1.png)

Effective consumption network.

``` r
plot(cm_alpha, type = "EffectiveProduction")
```

![Effective production
network.](visualisation_files/figure-html/plot-cm-effprod-1.png)

Effective production network.

## ConsortiumMetabolismSet plots

### Dendrogram

`plot(CMS)` draws the hierarchical clustering dendrogram with numbered
internal nodes. These node IDs are used by
[`extractCluster()`](https://admarhi.github.io/ramen/reference/extractCluster.md)
to pull sub-clusters.

``` r
plot(cms)
```

![CMS dendrogram with numbered
nodes.](visualisation_files/figure-html/plot-cms-dend-1.png)

CMS dendrogram with numbered nodes.

### Dendrogram with custom label colours

You can supply a tibble mapping leaf labels to colours.

``` r
colour_map <- data.frame(
    label = c("Alpha", "Beta", "Gamma", "Delta"),
    colour = c("steelblue", "firebrick",
               "forestgreen", "darkorange")
)
plot(cms, label_colours = colour_map)
```

![CMS dendrogram with coloured
labels.](visualisation_files/figure-html/plot-cms-colours-1.png)

CMS dendrogram with coloured labels.

## ConsortiumMetabolismAlignment plots

The CMA [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method
supports three types: `"heatmap"`, `"network"`, and `"scores"`. The
default depends on the alignment type: `"network"` for pairwise,
`"heatmap"` for multiple.

### Heatmap (multiple alignment)

The heatmap shows pairwise FOS similarities with dendrogram-based
ordering. Values range from 0 (no overlap) to 1 (identical).

``` r
plot(cma_mult, type = "heatmap")
```

![Similarity
heatmap.](visualisation_files/figure-html/plot-cma-heatmap-1.png)

Similarity heatmap.

### Network (pairwise alignment)

The network view visualises the alignment result as a directed
metabolite flow graph. Edge colours indicate the source: **green** =
shared pathways, **blue** = unique to query (Alpha), **red** = unique to
reference (Beta).

``` r
plot(cma_pair, type = "network")
```

![Pairwise alignment
network.](visualisation_files/figure-html/plot-cma-network-1.png)

Pairwise alignment network.

### Scores (pairwise)

The score bar chart displays all computed similarity metrics for a
pairwise alignment.

``` r
plot(cma_pair, type = "scores")
```

![Pairwise alignment
scores.](visualisation_files/figure-html/plot-cma-scores-pair-1.png)

Pairwise alignment scores.

### Scores (multiple alignment)

For a multiple alignment, the scores plot shows summary statistics
across all pairwise comparisons.

``` r
plot(cma_mult, type = "scores")
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_col()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_text()`).
```

![Multiple alignment
scores.](visualisation_files/figure-html/plot-cma-scores-mult-1.png)

Multiple alignment scores.

## Using `plotDirectedFlow()` directly

For fine-grained control over the network layout, you can call
[`plotDirectedFlow()`](https://admarhi.github.io/ramen/reference/plotDirectedFlow.md)
directly on an igraph object. This is the function underlying all CM and
CMA network plots.

``` r
g <- igraph::graph_from_adjacency_matrix(
    SummarizedExperiment::assay(cm_beta, "nSpecies"),
    mode = "directed",
    weighted = TRUE
)

plotDirectedFlow(
    g,
    color_edges_by_weight = TRUE,
    edge_width_range = c(0.5, 2),
    vertex_size = 12,
    vertex_label_cex = 0.7,
    main = "Beta - nSpecies (custom)"
)
```

![Custom directed flow
plot.](visualisation_files/figure-html/plot-directed-flow-1.png)

Custom directed flow plot.

## Summary of plot types

| Object class | `type` argument          | Description                           |
|--------------|--------------------------|---------------------------------------|
| CM           | `"Binary"` (default)     | Directed network, unweighted          |
| CM           | `"nSpecies"`             | Pathways weighted by species count    |
| CM           | `"Consumption"`          | Pathways weighted by consumption flux |
| CM           | `"Production"`           | Pathways weighted by production flux  |
| CM           | `"EffectiveConsumption"` | Effective consuming species           |
| CM           | `"EffectiveProduction"`  | Effective producing species           |
| CMS          | (none)                   | Dendrogram with numbered nodes        |
| CMA          | `"heatmap"`              | Pairwise similarity heatmap           |
| CMA          | `"network"`              | Shared/unique pathway graph           |
| CMA          | `"scores"`               | Bar chart of similarity metrics       |

## Session info

``` r
sessionInfo()
#> R version 4.5.3 (2026-03-11)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] ramen_0.0.0.9001 BiocStyle_2.38.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] SummarizedExperiment_1.40.0     gtable_0.3.6                   
#>  [3] ggplot2_4.0.2                   xfun_0.57                      
#>  [5] bslib_0.10.0                    Biobase_2.70.0                 
#>  [7] lattice_0.22-9                  yulab.utils_0.2.4              
#>  [9] vctrs_0.7.2                     tools_4.5.3                    
#> [11] generics_0.1.4                  stats4_4.5.3                   
#> [13] parallel_4.5.3                  tibble_3.3.1                   
#> [15] pkgconfig_2.0.3                 Matrix_1.7-4                   
#> [17] RColorBrewer_1.1-3              S7_0.2.1                       
#> [19] desc_1.4.3                      S4Vectors_0.48.0               
#> [21] lifecycle_1.0.5                 farver_2.1.2                   
#> [23] compiler_4.5.3                  treeio_1.34.0                  
#> [25] textshaping_1.0.5               Biostrings_2.78.0              
#> [27] Seqinfo_1.0.0                   codetools_0.2-20               
#> [29] htmltools_0.5.9                 sass_0.4.10                    
#> [31] yaml_2.3.12                     lazyeval_0.2.2                 
#> [33] pkgdown_2.2.0                   pillar_1.11.1                  
#> [35] crayon_1.5.3                    jquerylib_0.1.4                
#> [37] tidyr_1.3.2                     BiocParallel_1.44.0            
#> [39] SingleCellExperiment_1.32.0     DelayedArray_0.36.0            
#> [41] cachem_1.1.0                    viridis_0.6.5                  
#> [43] abind_1.4-8                     nlme_3.1-168                   
#> [45] tidyselect_1.2.1                digest_0.6.39                  
#> [47] dplyr_1.2.0                     purrr_1.2.1                    
#> [49] bookdown_0.46                   labeling_0.4.3                 
#> [51] TreeSummarizedExperiment_2.18.0 fastmap_1.2.0                  
#> [53] grid_4.5.3                      cli_3.6.5                      
#> [55] SparseArray_1.10.10             magrittr_2.0.4                 
#> [57] S4Arrays_1.10.1                 ape_5.8-1                      
#> [59] withr_3.0.2                     scales_1.4.0                   
#> [61] rappdirs_0.3.4                  rmarkdown_2.31                 
#> [63] XVector_0.50.0                  matrixStats_1.5.0              
#> [65] igraph_2.2.2                    gridExtra_2.3                  
#> [67] ragg_1.5.2                      evaluate_1.0.5                 
#> [69] knitr_1.51                      GenomicRanges_1.62.1           
#> [71] IRanges_2.44.0                  viridisLite_0.4.3              
#> [73] rlang_1.1.7                     dendextend_1.19.1              
#> [75] Rcpp_1.1.1                      glue_1.8.0                     
#> [77] tidytree_0.4.7                  BiocManager_1.30.27            
#> [79] BiocGenerics_0.56.0             jsonlite_2.0.0                 
#> [81] R6_2.6.1                        MatrixGenerics_1.22.0          
#> [83] systemfonts_1.3.2               fs_2.0.1
```
