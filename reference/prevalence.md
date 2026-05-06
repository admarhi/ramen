# Get Pathway Prevalence

Returns pathway prevalence across consortia from a multiple alignment.

## Usage

``` r
prevalence(object)

# S4 method for class 'ConsortiumMetabolismAlignment'
prevalence(object)
```

## Arguments

- object:

  A
  [ConsortiumMetabolismAlignment](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismAlignment.md)
  object of type `"multiple"`.

## Value

A data.frame with columns `consumed`, `produced`, `nConsortia`, and
`proportion`.

## Methods (by class)

- `prevalence(ConsortiumMetabolismAlignment)`: Pathway prevalence from a
  multiple
  [ConsortiumMetabolismAlignment](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismAlignment.md)

## Examples

``` r
cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
cms <- ConsortiumMetabolismSet(cm1, cm2, name = "test")
#> 
#> ── Creating CMS "test" ─────────────────────────────────────────────────────────
#> ℹ Validating 2 <ConsortiumMetabolism> objects
#> ✔ Validating 2 <ConsortiumMetabolism> objects [11ms]
#> 
#> ℹ Collecting metabolites from 2 consortia
#> ✔ Collecting metabolites from 2 consortia [31ms]
#> 
#> ℹ Re-indexing 7 unique metabolites
#> ✔ Re-indexing 7 unique metabolites [27ms]
#> 
#> ℹ Expanding 2 binary matrices to 7-dimensional space
#> ✔ Expanding 2 binary matrices to 7-dimensional space [24ms]
#> 
#> ℹ Computing 7 x 7 levels matrix
#> ✔ Computing 7 x 7 levels matrix [25ms]
#> 
#> ℹ Computing pairwise overlap (1 pairs via crossprod)
#> ✔ Computing pairwise overlap (1 pairs via crossprod) [23ms]
#> 
#> ℹ Assembling pathway data from 2 consortia
#> ✔ Assembling pathway data from 2 consortia [31ms]
#> 
#> ℹ Building dendrogram from 2 x 2 dissimilarity matrix
#> ✔ Building dendrogram from 2 x 2 dissimilarity matrix [22ms]
#> 
#> ℹ Extracting dendrogram node positions
#> ✔ Extracting dendrogram node positions [24ms]
#> 
#> ℹ Collecting 2 consortium graphs
#> CMS "test" created: 2 consortia, 7 metabolites (0.2s)
#> ✔ Collecting 2 consortium graphs [91ms]
#> 
cma <- align(cms)
#> Computing multiple alignment for 2 consortia using "FOS".
prevalence(cma)
#>    consumed produced nConsortia proportion
#> 1      met2     met1          1        0.5
#> 2      met4     met1          1        0.5
#> 3      met5     met1          2        1.0
#> 4      met2     met3          1        0.5
#> 5      met4     met3          1        0.5
#> 6      met5     met3          2        1.0
#> 7      met2     met4          1        0.5
#> 8      met5     met4          1        0.5
#> 9     media     met5          1        0.5
#> 10     met3     met5          1        0.5
#> 11     met6     met5          1        0.5
#> 12     met2     met6          1        0.5
#> 13     met5     met6          1        0.5
```
