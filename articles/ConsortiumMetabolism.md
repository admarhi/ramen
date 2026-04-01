# The ConsortiumMetabolism Class

## Introduction

The `ConsortiumMetabolism` (CM) class is the fundamental building block
of the `ramen` package. It represents a single microbial consortium’s
metabolic exchange network – the directed graph of metabolite
consumption and production by its member species.

A CM object stores:

- **Six assay matrices** (square, metabolite-by-metabolite, sparse):
  Binary, nSpecies, Consumption, Production, EffectiveConsumption,
  EffectiveProduction
- **Pathway data** linking consumed and produced metabolites with
  per-pathway summary statistics
- **An igraph object** of the directed metabolic network

This vignette covers construction, accessors, replacement methods,
subsetting, and inspection of CM objects.

``` r
library(ramen)
```

## Construction

### From an edge list

The
[`ConsortiumMetabolism()`](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md)
constructor accepts a data.frame (or tibble) with columns for species,
metabolite, and flux. Negative fluxes indicate consumption; positive
fluxes indicate production.

``` r
edge_data <- data.frame(
    species = c(
        "Sp_A", "Sp_A", "Sp_A", "Sp_A",
        "Sp_B", "Sp_B", "Sp_B",
        "Sp_C", "Sp_C", "Sp_C"
    ),
    met = c(
        "met1", "met2", "met3", "met4",
        "met1", "met3", "met4",
        "met2", "met3", "met4"
    ),
    flux = c(
        -0.30, 0.53, -2.23, 3.31,
        2.55, 0.34, -1.85,
        -1.30, -0.48, 0.60
    )
)

cm_manual <- ConsortiumMetabolism(
    edge_data,
    name = "manual_example"
)
cm_manual
#> 
#> ── ConsortiumMetabolism
#> Name: "manual_example"
#> Weighted metabolic network with 4 metabolites.
```

Custom column names can be specified via `species_col`,
`metabolite_col`, and `flux_col`:

``` r
df <- data.frame(
    organism = c("Sp_A", "Sp_A", "Sp_B", "Sp_B"),
    compound = c("glc", "ac", "glc", "ac"),
    rate = c(-1.0, 0.5, 0.8, -0.3)
)

cm_custom <- ConsortiumMetabolism(
    df,
    name = "custom_cols",
    species_col = "organism",
    metabolite_col = "compound",
    flux_col = "rate"
)
cm_custom
#> 
#> ── ConsortiumMetabolism
#> Name: "custom_cols"
#> Weighted metabolic network with 2 metabolites.
```

### With `synCM()`

For quick testing and demonstration,
[`synCM()`](https://admarhi.github.io/ramen/reference/synCM.md)
generates synthetic consortia with random species-metabolite
interactions.

``` r
cm_syn <- synCM(
    "SyntheticComm",
    n_species = 5,
    max_met = 8,
    seed = 42
)
cm_syn
#> 
#> ── ConsortiumMetabolism
#> Name: "SyntheticComm"
#> Weighted metabolic network with 8 metabolites.
```

To retrieve the raw edge list without constructing a CM object, set
`cm = FALSE`:

``` r
raw_edges <- synCM(
    "raw",
    n_species = 3,
    max_met = 5,
    seed = 1,
    cm = FALSE
)
head(raw_edges)
#> # A tibble: 6 × 3
#>   species metabolites fluxes
#>   <chr>   <chr>        <dbl>
#> 1 ERH100V met2        -4.79 
#> 2 ERH100V met1         0.989
#> 3 ERU091V met1         1.73 
#> 4 ERU091V met5        -0.916
#> 5 ERU091V met4         4.54 
#> 6 ERU091V met2         1.17
```

### From MiSoSoup data

The package ships with `misosoup24`, a list of 56 metabolic solutions
from [MiSoSoup](https://github.com/sirno/misosoup). Each element is a
tibble with columns `metabolites`, `species`, and `fluxes`.

``` r
data("misosoup24")
length(misosoup24)
#> [1] 56
names(misosoup24)[1:6]
#> [1] "ac_A1R12_1"  "ac_A1R12_10" "ac_A1R12_11" "ac_A1R12_12" "ac_A1R12_13"
#> [6] "ac_A1R12_14"

## Create a CM from the first solution
cm_ms <- ConsortiumMetabolism(
    misosoup24[[1]],
    name = names(misosoup24)[1],
    species_col = "species",
    metabolite_col = "metabolites",
    flux_col = "fluxes"
)
cm_ms
#> 
#> ── ConsortiumMetabolism
#> Name: "ac_A1R12_1"
#> Weighted metabolic network with 14 metabolites.
```

## The `show()` method

Printing a CM object gives a concise summary including the name, whether
the network is weighted, and the number of metabolites.

``` r
cm_syn
#> 
#> ── ConsortiumMetabolism
#> Name: "SyntheticComm"
#> Weighted metabolic network with 8 metabolites.
```

## Accessors

### `name()` and `description()`

``` r
name(cm_syn)
#> [1] "SyntheticComm"
description(cm_syn)
#> character(0)
```

### `species()`

Returns the unique species identifiers.

``` r
species(cm_syn)
#> [1] "MEU9286A" "FMD9122E" "IEV6830P" "DKC1921L" "WEL6294Y"
```

### `metabolites()`

Returns the unique metabolite identifiers.

``` r
metabolites(cm_syn)
#> [1] "met2" "met8" "met1" "met4" "met6" "met5" "met7" "met3"
```

### `pathways()`

Returns the pathway data.frame, which contains per-pathway summary
statistics including the number of species catalysing each pathway,
summed consumption and production fluxes, and effective diversity.

``` r
pw <- pathways(cm_syn)
pw[, c("consumed", "produced", "n_species")]
#> # A tibble: 24 × 3
#>    consumed produced n_species
#>    <chr>    <chr>        <dbl>
#>  1 met8     met2             2
#>  2 met8     met1             2
#>  3 met8     met6             1
#>  4 met4     met2             1
#>  5 met4     met1             1
#>  6 met4     met6             1
#>  7 met6     met5             2
#>  8 met6     met7             2
#>  9 met6     met8             2
#> 10 met6     met4             2
#> # ℹ 14 more rows
```

### `consortia()`

Returns the input data in a tidy format with columns `met`, `species`,
and `flux`.

``` r
head(consortia(cm_syn))
#>    met  species       flux
#> 1 met2 MEU9286A  4.6646866
#> 2 met8 MEU9286A -3.5629242
#> 3 met1 MEU9286A  0.4554387
#> 4 met4 MEU9286A -3.2583977
#> 5 met6 MEU9286A  4.8401184
#> 6 met5 FMD9122E  2.6453737
```

## Replacement methods

### `name<-`

``` r
name(cm_syn) <- "RenamedComm"
name(cm_syn)
#> [1] "RenamedComm"
```

### `description<-`

``` r
description(cm_syn) <- "A synthetic test community"
description(cm_syn)
#> [1] "A synthetic test community"
```

## Assay matrices

A CM object inherits from `TreeSummarizedExperiment` and stores six
square (metabolite-by-metabolite) sparse assay matrices. Rows represent
consumed metabolites and columns represent produced metabolites.

``` r
SummarizedExperiment::assayNames(cm_ms)
#> [1] "Binary"               "nSpecies"             "Consumption"         
#> [4] "Production"           "EffectiveConsumption" "EffectiveProduction"
```

Each matrix encodes a different view of the metabolic network:

| Assay                | Content                                  |
|----------------------|------------------------------------------|
| Binary               | 1 if a pathway exists, 0 otherwise       |
| nSpecies             | Number of species catalysing the pathway |
| Consumption          | Summed consumption flux across species   |
| Production           | Summed production flux across species    |
| EffectiveConsumption | Effective number of consuming species    |
| EffectiveProduction  | Effective number of producing species    |

``` r
## The Binary assay (sparse matrix)
SummarizedExperiment::assay(cm_ms, "Binary")
#> 14 x 14 sparse Matrix of class "dgCMatrix"
#>   [[ suppressing 14 column names 'ac', 'acald', 'ala__D' ... ]]
#>                                   
#> ac     . 1 . . 1 1 1 . 1 . 1 1 1 1
#> acald  1 . 1 1 . 1 . 1 . 1 . . . .
#> ala__D . 1 . . 1 1 1 . 1 . 1 1 1 1
#> ala__L . 1 . . 1 1 1 . 1 . 1 1 1 1
#> asp__L 1 . 1 1 . 1 . 1 . 1 . . . .
#> co2    . . . . . . . . . . . . . .
#> etoh   1 . 1 1 . 1 . 1 . 1 . . . .
#> glu__L . 1 . . 1 1 1 . 1 . 1 1 1 1
#> gly    1 . 1 1 . 1 . 1 . 1 . . . .
#> gthox  . 1 . . 1 1 1 . 1 . 1 1 1 1
#> gthrd  1 . 1 1 . 1 . 1 . 1 . . . .
#> h2o2   1 . 1 1 . 1 . 1 . 1 . . . .
#> h2s    1 . 1 1 . 1 . 1 . 1 . . . .
#> pyr    1 . 1 1 . 1 . 1 . 1 . . . .
```

## Dimensionality

Because CM inherits from `TreeSummarizedExperiment`, standard dimension
accessors work:

``` r
dim(cm_ms)
#> [1] 14 14
nrow(cm_ms)
#> [1] 14
ncol(cm_ms)
#> [1] 14
```

The rows and columns both correspond to the metabolite universe of this
consortium (including a synthetic “media” node when needed).

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
#>  [1] SummarizedExperiment_1.40.0     xfun_0.57                      
#>  [3] bslib_0.10.0                    Biobase_2.70.0                 
#>  [5] lattice_0.22-9                  yulab.utils_0.2.4              
#>  [7] vctrs_0.7.2                     tools_4.5.3                    
#>  [9] generics_0.1.4                  stats4_4.5.3                   
#> [11] parallel_4.5.3                  tibble_3.3.1                   
#> [13] pkgconfig_2.0.3                 Matrix_1.7-4                   
#> [15] RColorBrewer_1.1-3              desc_1.4.3                     
#> [17] S4Vectors_0.48.0                lifecycle_1.0.5                
#> [19] farver_2.1.2                    compiler_4.5.3                 
#> [21] treeio_1.34.0                   textshaping_1.0.5              
#> [23] Biostrings_2.78.0               Seqinfo_1.0.0                  
#> [25] codetools_0.2-20                htmltools_0.5.9                
#> [27] sass_0.4.10                     yaml_2.3.12                    
#> [29] lazyeval_0.2.2                  pkgdown_2.2.0                  
#> [31] pillar_1.11.1                   crayon_1.5.3                   
#> [33] jquerylib_0.1.4                 tidyr_1.3.2                    
#> [35] BiocParallel_1.44.0             SingleCellExperiment_1.32.0    
#> [37] DelayedArray_0.36.0             cachem_1.1.0                   
#> [39] abind_1.4-8                     nlme_3.1-168                   
#> [41] tidyselect_1.2.1                digest_0.6.39                  
#> [43] dplyr_1.2.0                     purrr_1.2.1                    
#> [45] bookdown_0.46                   TreeSummarizedExperiment_2.18.0
#> [47] fastmap_1.2.0                   grid_4.5.3                     
#> [49] cli_3.6.5                       SparseArray_1.10.10            
#> [51] magrittr_2.0.4                  S4Arrays_1.10.1                
#> [53] utf8_1.2.6                      ape_5.8-1                      
#> [55] withr_3.0.2                     scales_1.4.0                   
#> [57] rappdirs_0.3.4                  rmarkdown_2.31                 
#> [59] XVector_0.50.0                  matrixStats_1.5.0              
#> [61] igraph_2.2.2                    ragg_1.5.2                     
#> [63] evaluate_1.0.5                  knitr_1.51                     
#> [65] GenomicRanges_1.62.1            IRanges_2.44.0                 
#> [67] rlang_1.1.7                     Rcpp_1.1.1                     
#> [69] glue_1.8.0                      tidytree_0.4.7                 
#> [71] BiocManager_1.30.27             BiocGenerics_0.56.0            
#> [73] jsonlite_2.0.0                  R6_2.6.1                       
#> [75] MatrixGenerics_1.22.0           systemfonts_1.3.2              
#> [77] fs_2.0.1
```
