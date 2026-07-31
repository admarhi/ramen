# Coerce a ConsortiumMetabolismSet to a data.frame

Row-binds the per-consortium edge lists of every `ConsortiumMetabolism`
in the set, prefixing each row with a `consortium` column.

## Usage

``` r
# S4 method for class 'ConsortiumMetabolismSet'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A
  [`ConsortiumMetabolismSet`](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismSet.md)
  object.

- row.names:

  Ignored.

- optional:

  Ignored.

- ...:

  Additional arguments (currently unused).

## Value

A `data.frame` with columns `consortium`, `met`, `species`, and `flux`.
Empty sets return a 0-row `data.frame` with the same column names.

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
#> ✔ Collecting metabolites from 2 consortia [29ms]
#> 
#> ℹ Re-indexing 6 unique metabolites
#> ✔ Re-indexing 6 unique metabolites [26ms]
#> 
#> ℹ Expanding 2 binary matrices to 6-dimensional space
#> ✔ Expanding 2 binary matrices to 6-dimensional space [23ms]
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
#> ℹ Building dendrogram from 2 x 2 overlap matrix
#> ✔ Building dendrogram from 2 x 2 overlap matrix [22ms]
#> 
#> ℹ Extracting dendrogram node positions
#> ✔ Extracting dendrogram node positions [24ms]
#> 
#> ℹ Collecting 2 consortium graphs
#> CMS "demo" created: 2 consortia, 6 metabolites (0.2s)
#> ✔ Collecting 2 consortium graphs [88ms]
#> 
head(as.data.frame(cms))
#>   consortium  met  species      flux
#> 1          a met4 YSS7968P  1.641772
#> 2          a met4 FHC5781B -2.708004
#> 3          a met3 FHC5781B  2.814971
#> 4          a met3 JOR2426N -2.814971
#> 5          a met4 JOR2426N  1.066231
#> 6          a met1 FHC5781B -2.575173
```
