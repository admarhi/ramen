# Get the Constituent Consortia

Returns the `ConsortiumMetabolism` objects that a container is built
from. Defined for `ConsortiumMetabolismSet` (the consortia in the set)
and `ConsortiumMetabolismAlignment` (the consortia that produced the
alignment). The plural noun reflects that the result is always a
collection; for a single `ConsortiumMetabolism`, use
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) to obtain
the underlying edge list.

## Usage

``` r
consortia(object)

# S4 method for class 'ConsortiumMetabolismAlignment'
consortia(object)

# S4 method for class 'ConsortiumMetabolismSet'
consortia(object)
```

## Arguments

- object:

  A `ConsortiumMetabolismSet` object.

## Value

A named list of `ConsortiumMetabolism` objects.

## Methods (by class)

- `consortia(ConsortiumMetabolismAlignment)`: Not applicable to
  `ConsortiumMetabolismAlignment`. By design the class is a lightweight
  result object: it records the names of its inputs (`QueryName` /
  `ReferenceName`) but does not retain copies of the
  `ConsortiumMetabolism` objects themselves. Look up the originating CMs
  by name from the source `ConsortiumMetabolismSet`.

- `consortia(ConsortiumMetabolismSet)`: Get the list of
  `ConsortiumMetabolism` objects from a `ConsortiumMetabolismSet`.

## Examples

``` r
cm1 <- synCM("a", n_species = 3, max_met = 5)
cm2 <- synCM("b", n_species = 3, max_met = 5)
cms <- ConsortiumMetabolismSet(list(cm1, cm2), name = "demo")
#> 
#> ── Creating CMS "demo" ─────────────────────────────────────────────────────────
#> ℹ Validating 2 <ConsortiumMetabolism> objects
#> ✔ Validating 2 <ConsortiumMetabolism> objects [11ms]
#> 
#> ℹ Collecting metabolites from 2 consortia
#> ✔ Collecting metabolites from 2 consortia [30ms]
#> 
#> ℹ Re-indexing 6 unique metabolites
#> ✔ Re-indexing 6 unique metabolites [26ms]
#> 
#> ℹ Expanding 2 binary matrices to 6-dimensional space
#> ✔ Expanding 2 binary matrices to 6-dimensional space [24ms]
#> 
#> ℹ Computing 6 x 6 levels matrix
#> ✔ Computing 6 x 6 levels matrix [25ms]
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
#> ✔ Extracting dendrogram node positions [24ms]
#> 
#> ℹ Collecting 2 consortium graphs
#> CMS "demo" created: 2 consortia, 6 metabolites (0.2s)
#> ✔ Collecting 2 consortium graphs [90ms]
#> 
consortia(cms)
#> [[1]]
#> 
#> ── ConsortiumMetabolism 
#> Name: "a"
#> Weighted metabolic network: 3 species, 3 metabolites, 3 pathways.
#> Pathways per species: min 1, mean 1.3, max 2.
#> 
#> [[2]]
#> 
#> ── ConsortiumMetabolism 
#> Name: "b"
#> Weighted metabolic network: 3 species, 5 metabolites, 10 pathways.
#> Pathways per species: min 1, mean 3.7, max 6.
#> 
```
