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
#> ✔ Collecting metabolites from 2 consortia [30ms]
#> 
#> ℹ Re-indexing 5 unique metabolites
#> ✔ Re-indexing 5 unique metabolites [36ms]
#> 
#> ℹ Expanding 2 binary matrices to 5-dimensional space
#> ✔ Expanding 2 binary matrices to 5-dimensional space [25ms]
#> 
#> ℹ Computing 5 x 5 levels matrix
#> ✔ Computing 5 x 5 levels matrix [25ms]
#> 
#> ℹ Computing pairwise overlap (1 pairs via crossprod)
#> ✔ Computing pairwise overlap (1 pairs via crossprod) [23ms]
#> 
#> ℹ Assembling pathway data from 2 consortia
#> ✔ Assembling pathway data from 2 consortia [30ms]
#> 
#> ℹ Building dendrogram from 2 x 2 overlap matrix
#> ✔ Building dendrogram from 2 x 2 overlap matrix [21ms]
#> 
#> ℹ Extracting dendrogram node positions
#> ✔ Extracting dendrogram node positions [23ms]
#> 
#> ℹ Collecting 2 consortium graphs
#> CMS "test" created: 2 consortia, 5 metabolites (0.2s)
#> ✔ Collecting 2 consortium graphs [84ms]
#> 
show(cms)
#> 
#> ── ConsortiumMetabolismSet 
#> Name: "test"
#> 2 consortia, 7 species, 5 metabolites.
#> Community size (species): min 3, mean 3.5, max 4.
#> Community size (metabolites): min 4, mean 4.5, max 5.
#> Pathways: 1 pan-cons, 10 niche, 0 core, 8 aux (quantileCutoff = 0.1; buckets
#> may overlap, see ?pathways).
#> Species: 2 generalists, 4 specialists (quantileCutoff = 0.15).
```
