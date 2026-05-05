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
  from a
  [ConsortiumMetabolismAlignment](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismAlignment.md).
  For multiple alignments this is the symmetric n x n pairwise matrix;
  for database search alignments it is a 1 x n row vector of the query's
  scores against each database member.

## Examples

``` r
cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
cms <- ConsortiumMetabolismSet(cm1, cm2, name = "test")
#> 
#> ── Creating CMS "test" ─────────────────────────────────────────────────────────
#> ℹ Validating 2 <ConsortiumMetabolism> objects
#> ✔ Validating 2 <ConsortiumMetabolism> objects [12ms]
#> 
#> ℹ Collecting metabolites from 2 consortia
#> ✔ Collecting metabolites from 2 consortia [32ms]
#> 
#> ℹ Re-indexing 7 unique metabolites
#> ✔ Re-indexing 7 unique metabolites [30ms]
#> 
#> ℹ Expanding 2 binary matrices to 7-dimensional space
#> ✔ Expanding 2 binary matrices to 7-dimensional space [25ms]
#> 
#> ℹ Computing 7 x 7 levels matrix
#> ✔ Computing 7 x 7 levels matrix [26ms]
#> 
#> ℹ Computing pairwise overlap (1 pairs via crossprod)
#> ✔ Computing pairwise overlap (1 pairs via crossprod) [26ms]
#> 
#> ℹ Assembling pathway data from 2 consortia
#> ✔ Assembling pathway data from 2 consortia [32ms]
#> 
#> ℹ Building dendrogram from 2 x 2 dissimilarity matrix
#> ✔ Building dendrogram from 2 x 2 dissimilarity matrix [23ms]
#> 
#> ℹ Extracting dendrogram node positions
#> ✔ Extracting dendrogram node positions [26ms]
#> 
#> ℹ Collecting 2 consortium graphs
#> CMS "test" created: 2 consortia, 7 metabolites (0.2s)
#> ✔ Collecting 2 consortium graphs [89ms]
#> 
cma <- align(cms)
#> Computing multiple alignment for 2 consortia using "FOS".
similarityMatrix(cma)
#>        comm_1 comm_2
#> comm_1    1.0    0.5
#> comm_2    0.5    1.0
```
