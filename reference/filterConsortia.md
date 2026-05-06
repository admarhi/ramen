# Filter Consortia from a Set

Selects a subset of consortia from a `ConsortiumMetabolismSet`,
returning a fully recomputed `ConsortiumMetabolismSet` containing only
the selected consortia.

## Usage

``` r
filterConsortia(object, i)

# S4 method for class 'ConsortiumMetabolismSet'
filterConsortia(object, i)
```

## Arguments

- object:

  A `ConsortiumMetabolismSet` object.

- i:

  Integer vector, character vector of consortium names, or logical
  vector.

## Value

A `ConsortiumMetabolismSet` containing only the selected consortia.

## Methods (by class)

- `filterConsortia(ConsortiumMetabolismSet)`: Filter consortia from a
  `ConsortiumMetabolismSet`

## Examples

``` r
cm1 <- synCM("a", n_species = 3, max_met = 5)
cm2 <- synCM("b", n_species = 3, max_met = 5)
cms <- ConsortiumMetabolismSet(cm1, cm2, name = "test")
#> 
#> ── Creating CMS "test" ─────────────────────────────────────────────────────────
#> ℹ Validating 2 <ConsortiumMetabolism> objects
#> ✔ Validating 2 <ConsortiumMetabolism> objects [10ms]
#> 
#> ℹ Collecting metabolites from 2 consortia
#> ✔ Collecting metabolites from 2 consortia [28ms]
#> 
#> ℹ Re-indexing 5 unique metabolites
#> ✔ Re-indexing 5 unique metabolites [25ms]
#> 
#> ℹ Expanding 2 binary matrices to 5-dimensional space
#> ✔ Expanding 2 binary matrices to 5-dimensional space [22ms]
#> 
#> ℹ Computing 5 x 5 levels matrix
#> ✔ Computing 5 x 5 levels matrix [23ms]
#> 
#> ℹ Computing pairwise overlap (1 pairs via crossprod)
#> ✔ Computing pairwise overlap (1 pairs via crossprod) [22ms]
#> 
#> ℹ Assembling pathway data from 2 consortia
#> ✔ Assembling pathway data from 2 consortia [29ms]
#> 
#> ℹ Building dendrogram from 2 x 2 dissimilarity matrix
#> ✔ Building dendrogram from 2 x 2 dissimilarity matrix [21ms]
#> 
#> ℹ Extracting dendrogram node positions
#> ✔ Extracting dendrogram node positions [23ms]
#> 
#> ℹ Collecting 2 consortium graphs
#> CMS "test" created: 2 consortia, 5 metabolites (0.2s)
#> ✔ Collecting 2 consortium graphs [78ms]
#> 
filterConsortia(cms, 1L)
#> 
#> ── ConsortiumMetabolismSet 
#> Name: "test"
#> 1 consortia, 3 species, 5 metabolites.
#> Community size (species): min 3, mean 3, max 3.
#> Community size (metabolites): min 5, mean 5, max 5.
```
