# Get Similarity Matrix

Returns the pairwise similarity matrix from a multiple alignment.

## Usage

``` r
similarityMatrix(object)

# S4 method for class 'ConsortiumMetabolismAlignment'
similarityMatrix(object)
```

## Arguments

- object:

  A
  [ConsortiumMetabolismAlignment](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismAlignment.md)
  object of type `"multiple"`.

## Value

A numeric n x n matrix.

## Methods (by class)

- `similarityMatrix(ConsortiumMetabolismAlignment)`: Similarity matrix
  from a multiple
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
#> ✔ Building dendrogram from 2 x 2 dissimilarity matrix [22ms]
#> 
#> ℹ Extracting dendrogram node positions
#> ✔ Extracting dendrogram node positions [24ms]
#> 
#> ℹ Collecting 2 consortium graphs
#> CMS "test" created: 2 consortia, 6 metabolites (0.2s)
#> ✔ Collecting 2 consortium graphs [79ms]
#> 
cma <- align(cms)
#> Computing multiple alignment for 2 consortia using "FOS".
similarityMatrix(cma)
#>        comm_1 comm_2
#> comm_1    1.0    0.3
#> comm_2    0.3    1.0
```
