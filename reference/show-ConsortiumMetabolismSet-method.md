# Show Method for `ConsortiumMetabolismSet` Object

Show Method for `ConsortiumMetabolismSet` Object

## Usage

``` r
# S4 method for class 'ConsortiumMetabolismSet'
show(object)
```

## Arguments

- object:

  An object of class `ConsortiumMetabolismSet`

## Value

The object, invisibly.

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
#> ✔ Collecting metabolites from 2 consortia [31ms]
#> 
#> ℹ Re-indexing 7 unique metabolites
#> ✔ Re-indexing 7 unique metabolites [28ms]
#> 
#> ℹ Expanding 2 binary matrices to 7-dimensional space
#> ✔ Expanding 2 binary matrices to 7-dimensional space [24ms]
#> 
#> ℹ Computing 7 x 7 levels matrix
#> ✔ Computing 7 x 7 levels matrix [25ms]
#> 
#> ℹ Computing pairwise overlap (1 pairs via crossprod)
#> ✔ Computing pairwise overlap (1 pairs via crossprod) [25ms]
#> 
#> ℹ Assembling pathway data from 2 consortia
#> ✔ Assembling pathway data from 2 consortia [32ms]
#> 
#> ℹ Building dendrogram from 2 x 2 dissimilarity matrix
#> ✔ Building dendrogram from 2 x 2 dissimilarity matrix [24ms]
#> 
#> ℹ Extracting dendrogram node positions
#> ✔ Extracting dendrogram node positions [25ms]
#> 
#> ℹ Collecting 2 consortium graphs
#> CMS "test" created: 2 consortia, 7 metabolites (0.2s)
#> ✔ Collecting 2 consortium graphs [99ms]
#> 
show(cms)
#> 
#> ── ConsortiumMetabolismSet 
#> Name: "test"
#> 2 consortia, 7 species, 7 metabolites.
#> Community size (species): min 3, mean 3.5, max 4.
#> Community size (metabolites): min 3, mean 4.5, max 6.
#> Pathways: 0 pan-cons, 17 niche, 0 core, 10 aux (quantile = 0.1).
#> Species: 1 generalists, 2 specialists (quantile = 0.15).
```
