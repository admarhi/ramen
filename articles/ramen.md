# Introduction to ramen

## Installation

`ramen` can be installed from Bioconductor:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ramen")
```

``` r
library(ramen)
```

## Overview

`ramen` (**R**econstruction and **A**lignment of **M**icrobial
**E**xchange **N**etworks) analyses the functional similarity of
microbial communities from a metabolic network perspective. Rather than
comparing communities by species composition, `ramen` asks: *which
metabolite-to-metabolite pathways does each community catalyse, and how
similar are the resulting networks?*

### Key concept: pathways

Throughout `ramen`, a **pathway** is a directed edge from one metabolite
to another: “metabolite A is consumed and metabolite B is produced by at
least one species.” This is *not* a biochemical pathway in the KEGG or
MetaCyc sense (e.g. glycolysis). Rather, it captures a single metabolic
exchange coupling within the community. A consortium’s full set of
pathways forms a square metabolite-by-metabolite network that encodes
its collective metabolic capability.

### Package classes

The package provides three S4 classes that build on Bioconductor’s
*[TreeSummarizedExperiment](https://bioconductor.org/packages/3.22/TreeSummarizedExperiment)*:

- **`ConsortiumMetabolism`** (CM) – a single community’s metabolic
  exchange network
- **`ConsortiumMetabolismSet`** (CMS) – a collection of CMs with
  precomputed overlap scores and a dendrogram
- **`ConsortiumMetabolismAlignment`** (CMA) – the result of aligning two
  or more communities

This vignette walks through a complete analysis from data import to
alignment using the bundled `misosoup24` dataset (56 metabolic solutions
from [MiSoSoup](https://github.com/sirno/misosoup)). For details on
alignment metrics, see
[`vignette("alignment", package = "ramen")`](https://admarhi.github.io/ramen/articles/alignment.md);
for a gallery of all plot types, see
[`vignette("visualisation", package = "ramen")`](https://admarhi.github.io/ramen/articles/visualisation.md).

## Data import

### The `misosoup24` dataset

`ramen` ships with `misosoup24`, a list of 56 metabolic solutions. Each
element is a data.frame with columns `metabolites`, `species`, and
`fluxes`, where negative fluxes indicate consumption and positive fluxes
indicate production.

``` r
data("misosoup24")
length(misosoup24)
#> [1] 56
names(misosoup24)[1:8]
#> [1] "ac_A1R12_1"  "ac_A1R12_10" "ac_A1R12_11" "ac_A1R12_12" "ac_A1R12_13"
#> [6] "ac_A1R12_14" "ac_A1R12_15" "ac_A1R12_16"
head(misosoup24[[1]])
#> # A tibble: 6 × 3
#>   metabolite species    flux
#>   <chr>      <chr>     <dbl>
#> 1 ac         A1R12     0.773
#> 2 ac         I2R16   -10.8  
#> 3 acald      A1R12    -1.12 
#> 4 acald      I2R16     1.12 
#> 5 ala__D     A1R12     0.760
#> 6 ala__D     I2R16    -0.760
```

### Importing raw MiSoSoup YAML

For raw MiSoSoup output (nested YAML),
[`importMisosoup()`](https://admarhi.github.io/ramen/reference/importMisosoup.md)
parses the data into structured tibbles of consortia, media, and growth
information.
[`overviewMisosoup()`](https://admarhi.github.io/ramen/reference/overviewMisosoup.md)
provides a quick summary of the data structure before full import.

``` r
## Parse raw MiSoSoup YAML (not run -- requires external data)
raw <- yaml::read_yaml("path/to/misosoup_output.yaml")
overviewMisosoup(raw)
result <- importMisosoup(raw)
names(result) # "consortia", "media", "growth"
```

### Alternative input formats

The
[`ConsortiumMetabolism()`](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md)
constructor accepts any data.frame with columns `species`, `metabolite`,
and `flux`. Custom column names can be specified via `species_col`,
`metabolite_col`, and `flux_col`.

If your data is in wide format with separate columns for consumed and
produced metabolites,
[`pivotCM()`](https://admarhi.github.io/ramen/reference/pivotCM.md)
converts it to the long format that
[`ConsortiumMetabolism()`](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md)
expects:

``` r
wide_data <- data.frame(
    species = c("Sp_A", "Sp_B", "Sp_C"),
    uptake = c("met1", "met2", "met3"),
    secretion = c("met2", "met3", "met1"),
    flux = c(1, 1, 1)
)

long_data <- pivotCM(
    wide_data,
    species = "species",
    from = "uptake",
    to = "secretion",
    flux = "flux"
)
head(long_data)
#> # A tibble: 6 × 3
#>   species met    flux
#>   <chr>   <chr> <dbl>
#> 1 Sp_A    met1     -1
#> 2 Sp_A    met2      1
#> 3 Sp_B    met2     -1
#> 4 Sp_B    met3      1
#> 5 Sp_C    met3     -1
#> 6 Sp_C    met1      1
```

## Building a ConsortiumMetabolism

### From `misosoup24`

``` r
cm1 <- ConsortiumMetabolism(
    misosoup24[[1]],
    name = names(misosoup24)[1]
)
cm1
#> 
#> ── ConsortiumMetabolism
#> Name: "ac_A1R12_1"
#> Weighted metabolic network: 2 species, 14 metabolites, 93 pathways.
```

### Synthetic data with `synCM()`

For quick testing,
[`synCM()`](https://admarhi.github.io/ramen/reference/synCM.md)
generates a random consortium:

``` r
cm_syn <- synCM("Synthetic", n_species = 4, max_met = 8, seed = 42)
cm_syn
#> 
#> ── ConsortiumMetabolism
#> Name: "Synthetic"
#> Weighted metabolic network: 4 species, 8 metabolites, 13 pathways.
```

### Inspecting a CM

A CM inherits from `TreeSummarizedExperiment`, so standard dimension
accessors work. Rows and columns both correspond to metabolites (it is a
square metabolite-by-metabolite matrix):

``` r
dim(cm1)
#> [1] 14 14
```

The six assay matrices encode different views of the network:

``` r
SummarizedExperiment::assayNames(cm1)
#> [1] "Binary"               "nSpecies"             "Consumption"         
#> [4] "Production"           "EffectiveConsumption" "EffectiveProduction"
```

| Assay                | Content                                  |
|----------------------|------------------------------------------|
| Binary               | 1 if a pathway exists, 0 otherwise       |
| nSpecies             | Number of species catalysing the pathway |
| Consumption          | Summed consumption flux                  |
| Production           | Summed production flux                   |
| EffectiveConsumption | Effective number of consuming species    |
| EffectiveProduction  | Effective number of producing species    |

### Accessors

``` r
name(cm1)
#> [1] "ac_A1R12_1"
species(cm1)
#> [1] "A1R12" "I2R16"
metabolites(cm1)
#>  [1] "ac"     "acald"  "ala__D" "ala__L" "asp__L" "co2"    "etoh"   "glu__L"
#>  [9] "gly"    "gthox"  "gthrd"  "h2o2"   "h2s"    "pyr"
```

[`pathways()`](https://admarhi.github.io/ramen/reference/pathways.md)
returns per-pathway summary statistics:

``` r
head(pathways(cm1))
#> # A tibble: 6 × 3
#>   consumed produced n_species
#>   <chr>    <chr>        <dbl>
#> 1 ac       acald            1
#> 2 ac       asp__L           1
#> 3 ac       co2              1
#> 4 ac       etoh             1
#> 5 ac       gly              1
#> 6 ac       gthrd            1
```

[`consortia()`](https://admarhi.github.io/ramen/reference/consortia.md)
returns the tidy input data:

``` r
head(consortia(cm1))
#>      met species        flux
#> 1     ac   A1R12   0.7729250
#> 2     ac   I2R16 -10.7673345
#> 3  acald   A1R12  -1.1209639
#> 4  acald   I2R16   1.1209639
#> 5 ala__D   A1R12   0.7603136
#> 6 ala__D   I2R16  -0.7603136
```

Replacement methods allow renaming:

``` r
name(cm1) <- "first_solution"
name(cm1)
#> [1] "first_solution"
```

## Building a ConsortiumMetabolismSet

A `ConsortiumMetabolismSet` groups multiple CMs. During construction,
`ramen` automatically expands all binary matrices to a shared universal
metabolite space, computes pairwise overlap scores, and builds a
hierarchical clustering dendrogram.

``` r
## Build 20 CMs from misosoup24
cm_list <- lapply(seq_len(20), function(i) {
    ConsortiumMetabolism(
        misosoup24[[i]],
        name = names(misosoup24)[i]
    )
})

cms <- ConsortiumMetabolismSet(cm_list, name = "MiSoSoup_20")
cms
```

### Species across the set

[`species()`](https://admarhi.github.io/ramen/reference/species.md)
returns all species with the number of distinct pathways each catalyses:

``` r
species(cms)
#> # A tibble: 23 × 2
#>    species n_pathways
#>    <chr>        <int>
#>  1 A1R12          302
#>  2 D2M19          132
#>  3 E3M18          108
#>  4 F3R08           99
#>  5 m_3C02          96
#>  6 C2R02           91
#>  7 m_3D05          91
#>  8 C3R12           88
#>  9 I3M07           88
#> 10 B3M02           88
#> # ℹ 13 more rows
```

Filter by metabolic role – generalists (top fraction by pathway count)
or specialists (bottom fraction):

``` r
species(cms, type = "generalists")
#> # A tibble: 4 × 2
#>   species n_pathways
#>   <chr>        <int>
#> 1 A1R12          302
#> 2 D2M19          132
#> 3 E3M18          108
#> 4 F3R08           99
species(cms, type = "specialists")
#> # A tibble: 4 × 2
#>   species n_pathways
#>   <chr>        <int>
#> 1 B3R10           50
#> 2 I2R16           45
#> 3 A3R04           36
#> 4 C1M14           35
```

### Pathway classification

[`pathways()`](https://admarhi.github.io/ramen/reference/pathways.md)
with a `type` argument classifies pathways by prevalence:

``` r
## Pan-consortia: present in many consortia (top fraction)
head(pathways(cms, type = "pan-cons"))
#> # A tibble: 6 × 4
#>   consumed produced n_species n_cons
#>   <chr>    <chr>        <int>  <int>
#> 1 ac       pyr             20     20
#> 2 ala__L   pyr             20     20
#> 3 etoh     co2              5     20
#> 4 h2s      co2             15     20
#> 5 pyr      ac               1     20
#> 6 pyr      ala__L           1     20

## Niche: present in few consortia (bottom fraction)
head(pathways(cms, type = "niche"))
#> # A tibble: 6 × 4
#>   consumed     produced n_species n_cons
#>   <chr>        <chr>        <int>  <int>
#> 1 LalaDgluMdap ac               1      1
#> 2 LalaDgluMdap ala__D           1      1
#> 3 LalaDgluMdap ala__L           1      1
#> 4 LalaDgluMdap co2              1      1
#> 5 LalaDgluMdap etoh             1      1
#> 6 LalaDgluMdap gthox            1      1
```

The `quantileCutoff` parameter controls the threshold (default 0.1).

### Functional groups

[`functionalGroups()`](https://admarhi.github.io/ramen/reference/functionalGroups.md)
clusters species by the Jaccard similarity of their pathway sets.
[`plotFunctionalGroups()`](https://admarhi.github.io/ramen/reference/plotFunctionalGroups.md)
visualizes the resulting dendrogram:

``` r
fg <- functionalGroups(cms)
plotFunctionalGroups(fg, k = 3)
#> Loading required namespace: colorspace
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the dendextend package.
#>   Please report the issue at <https://github.com/talgalili/dendextend/issues>.
#> This warning is displayed once per session.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.
```

![Functional groups
dendrogram.](ramen_files/figure-html/cms-func-groups-1.png)

Functional groups dendrogram.

### Dendrogram and cluster extraction

``` r
plot(cms)
```

![CMS dendrogram with numbered
nodes.](ramen_files/figure-html/cms-dendrogram-1.png)

CMS dendrogram with numbered nodes.

The numbered internal nodes can be used to extract sub-clusters:

``` r
sub_cms <- extractCluster(cms, node_id = 1)
sub_cms
```

## Alignment

### Pairwise alignment

[`align()`](https://admarhi.github.io/ramen/reference/align.md) compares
two CMs and returns a `ConsortiumMetabolismAlignment`:

``` r
cma_pair <- align(cm_list[[1]], cm_list[[2]])
cma_pair
#> 
#> ── ConsortiumMetabolismAlignment
#> Name: NA
#> Type: "pairwise"
#> Metric: "FOS"
#> Score: 0.7634
#> Query: "ac_A1R12_1", Reference: "ac_A1R12_10"
#> Coverage: query 0.763, reference 0.38
```

All similarity metrics are computed automatically:

``` r
scores(cma_pair)
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

Pathway correspondences show shared and unique pathways:

``` r
shared <- pathways(cma_pair, type = "shared")
unique_pw <- pathways(cma_pair, type = "unique")
nrow(shared)
#> [1] 71
nrow(unique_pw$query)
#> [1] 22
nrow(unique_pw$reference)
#> [1] 116
```

### Multiple alignment

[`align()`](https://admarhi.github.io/ramen/reference/align.md) on a CMS
computes all pairwise similarities:

``` r
cma_mult <- align(cms)
cma_mult
```

The similarity matrix and summary scores:

``` r
scores(cma_mult)
#> $mean
#> [1] 0.5858039
#> 
#> $median
#> [1] 0.5651547
#> 
#> $min
#> [1] 0.2362205
#> 
#> $max
#> [1] 1
#> 
#> $sd
#> [1] 0.1764926
#> 
#> $nPairs
#> [1] 190
round(similarityMatrix(cma_mult)[1:5, 1:5], 3)
#>             ac_A1R12_1 ac_A1R12_10 ac_A1R12_11 ac_A1R12_12 ac_A1R12_13
#> ac_A1R12_1       1.000       0.763       0.538       0.785       0.538
#> ac_A1R12_10      0.763       1.000       0.438       0.542       0.549
#> ac_A1R12_11      0.538       0.438       1.000       0.507       1.000
#> ac_A1R12_12      0.785       0.542       0.507       1.000       0.430
#> ac_A1R12_13      0.538       0.549       1.000       0.430       1.000
```

Pathway prevalence across consortia:

``` r
prev <- prevalence(cma_mult)
head(prev[order(-prev$nConsortia), ])
#>     consumed produced nConsortia proportion
#> 59       pyr       ac         20          1
#> 139      pyr   ala__L         20          1
#> 222     etoh      co2         20          1
#> 231      h2s      co2         20          1
#> 238      pyr      co2         20          1
#> 586       ac      pyr         20          1
```

For a detailed treatment of alignment metrics, p-values, and accessor
functions, see
[`vignette("alignment", package = "ramen")`](https://admarhi.github.io/ramen/articles/alignment.md).

## Visualization

Each class has a
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method. Here is
one example per class; for the full gallery, see
[`vignette("visualisation", package = "ramen")`](https://admarhi.github.io/ramen/articles/visualisation.md).

``` r
plot(cm_list[[1]], type = "Binary")
```

![Metabolic network for the first
consortium.](ramen_files/figure-html/plot-cm-1.png)

Metabolic network for the first consortium.

``` r
plot(cms)
```

![Dendrogram of 20 consortia.](ramen_files/figure-html/plot-cms-1.png)

Dendrogram of 20 consortia.

``` r
plot(cma_mult, type = "heatmap")
```

![Pairwise similarity heatmap.](ramen_files/figure-html/plot-cma-1.png)

Pairwise similarity heatmap.

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
#> [1] ramen_0.99.0     BiocStyle_2.38.0
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
#> [19] desc_1.4.3                      S4Vectors_0.48.1               
#> [21] lifecycle_1.0.5                 farver_2.1.2                   
#> [23] compiler_4.5.3                  treeio_1.34.0                  
#> [25] textshaping_1.0.5               Biostrings_2.78.0              
#> [27] Seqinfo_1.0.0                   codetools_0.2-20               
#> [29] htmltools_0.5.9                 sass_0.4.10                    
#> [31] yaml_2.3.12                     lazyeval_0.2.3                 
#> [33] pkgdown_2.2.0                   pillar_1.11.1                  
#> [35] crayon_1.5.3                    jquerylib_0.1.4                
#> [37] tidyr_1.3.2                     BiocParallel_1.44.0            
#> [39] SingleCellExperiment_1.32.0     DelayedArray_0.36.1            
#> [41] cachem_1.1.0                    viridis_0.6.5                  
#> [43] abind_1.4-8                     nlme_3.1-168                   
#> [45] tidyselect_1.2.1                digest_0.6.39                  
#> [47] dplyr_1.2.1                     purrr_1.2.1                    
#> [49] bookdown_0.46                   labeling_0.4.3                 
#> [51] TreeSummarizedExperiment_2.18.0 fastmap_1.2.0                  
#> [53] grid_4.5.3                      cli_3.6.5                      
#> [55] SparseArray_1.10.10             magrittr_2.0.5                 
#> [57] S4Arrays_1.10.1                 utf8_1.2.6                     
#> [59] ape_5.8-1                       withr_3.0.2                    
#> [61] scales_1.4.0                    rappdirs_0.3.4                 
#> [63] rmarkdown_2.31                  XVector_0.50.0                 
#> [65] matrixStats_1.5.0               igraph_2.2.3                   
#> [67] gridExtra_2.3                   ragg_1.5.2                     
#> [69] evaluate_1.0.5                  knitr_1.51                     
#> [71] GenomicRanges_1.62.1            IRanges_2.44.0                 
#> [73] viridisLite_0.4.3               rlang_1.2.0                    
#> [75] dendextend_1.19.1               Rcpp_1.1.1                     
#> [77] glue_1.8.0                      tidytree_0.4.7                 
#> [79] BiocManager_1.30.27             BiocGenerics_0.56.0            
#> [81] jsonlite_2.0.0                  R6_2.6.1                       
#> [83] MatrixGenerics_1.22.0           systemfonts_1.3.2              
#> [85] fs_2.0.1
```
