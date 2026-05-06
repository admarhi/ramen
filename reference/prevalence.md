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
#> ✔ Collecting metabolites from 2 consortia [29ms]
#> 
#> ℹ Re-indexing 6 unique metabolites
#> ✔ Re-indexing 6 unique metabolites [26ms]
#> 
#> ℹ Expanding 2 binary matrices to 6-dimensional space
#> ✔ Expanding 2 binary matrices to 6-dimensional space [23ms]
#> 
#> ℹ Computing 6 x 6 levels matrix
#> ✔ Computing 6 x 6 levels matrix [24ms]
#> 
#> ℹ Computing pairwise overlap (1 pairs via crossprod)
#> ✔ Computing pairwise overlap (1 pairs via crossprod) [23ms]
#> 
#> ℹ Assembling pathway data from 2 consortia
#> ✔ Assembling pathway data from 2 consortia [30ms]
#> 
#> ℹ Building dendrogram from 2 x 2 dissimilarity matrix
#> ✔ Building dendrogram from 2 x 2 dissimilarity matrix [21ms]
#> 
#> ℹ Extracting dendrogram node positions
#> ✔ Extracting dendrogram node positions [23ms]
#> 
#> ℹ Collecting 2 consortium graphs
#> CMS "test" created: 2 consortia, 6 metabolites (0.2s)
#> ✔ Collecting 2 consortium graphs [86ms]
#> 
cma <- align(cms)
#> Computing multiple alignment for 2 consortia using "FOS".
prevalence(cma)
#>    consumed produced nConsortia proportion
#> 1      met2     met1          2        1.0
#> 2      met3     met1          1        0.5
#> 3      met4     met1          1        0.5
#> 4      met5     met1          1        0.5
#> 5      met6     met1          1        0.5
#> 6      met1     met2          2        1.0
#> 7      met4     met2          2        1.0
#> 8      met5     met2          1        0.5
#> 9      met1     met3          2        1.0
#> 10     met2     met3          1        0.5
#> 11     met4     met3          2        1.0
#> 12     met5     met3          2        1.0
#> 13     met1     met4          1        0.5
#> 14     met2     met4          2        1.0
#> 15     met3     met4          1        0.5
#> 16     met6     met4          1        0.5
#> 17     met1     met5          1        0.5
#> 18     met2     met5          1        0.5
#> 19     met3     met5          1        0.5
#> 20     met4     met5          1        0.5
#> 21     met6     met5          1        0.5
#> 22     met1     met6          1        0.5
#> 23     met4     met6          1        0.5
#> 24     met5     met6          1        0.5
```
