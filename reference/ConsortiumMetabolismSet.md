# Set of `ConsortiumMetabolism` Objects

Creates a `ConsortiumMetabolismSet` combining multiple
`ConsortiumMetabolism` objects into a unified metabolite space. Computes
pairwise overlap scores and builds a dendrogram for clustering.

## Usage

``` r
ConsortiumMetabolismSet(
  ...,
  name = NA_character_,
  desc = NA_character_,
  linkage = "complete",
  verbose = TRUE
)
```

## Arguments

- ...:

  Lists or individual `ConsortiumMetabolism` objects.

- name:

  Character scalar giving the name of the set.

- desc:

  Optional short description of the set.

- linkage:

  Character scalar specifying the agglomeration method for hierarchical
  clustering of the overlap matrix. Passed to
  [`hclust`](https://rdrr.io/r/stats/hclust.html) as the `method`
  argument. Defaults to `"complete"`, which produces compact clusters
  where every pair of consortia within a cluster has dissimilarity below
  the merge threshold.

- verbose:

  Logical scalar. If `TRUE` (default), prints progress messages during
  construction.

## Value

A `ConsortiumMetabolismSet` object.

## Slots

- `Name`:

  character. Display name for the set.

- `Consortia`:

  list. List of `ConsortiumMetabolism` objects.

- `Description`:

  character. Optional short description.

- `OverlapMatrix`:

  matrix. Pairwise dissimilarity matrix (1 - overlap) between consortia.

- `Dendrogram`:

  list. Hierarchical clustering dendrogram.

- `NodeData`:

  data.frame. Internal node positions from the dendrogram.

- `Graphs`:

  list. Named list of igraph objects, one per consortium.

- `BinaryMatrices`:

  list. Named list of binary matrices expanded to universal metabolite
  space.

- `Pathways`:

  data.frame. Combined pathway list from all consortia with re-indexed
  metabolite positions.

- `Metabolites`:

  data.frame. Metabolite mapping between per-consortium and universal
  indices.

## See also

[ConsortiumMetabolism](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md),
[TreeSummarizedExperiment-class](https://rdrr.io/pkg/TreeSummarizedExperiment/man/TreeSummarizedExperiment-class.html)

## Examples

``` r
cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
cms <- ConsortiumMetabolismSet(cm1, cm2, name = "example")
#> 
#> ── Creating CMS "example" ──────────────────────────────────────────────────────
#> ℹ Validating 2 <ConsortiumMetabolism> objects
#> ✔ Validating 2 <ConsortiumMetabolism> objects [16ms]
#> 
#> ℹ Collecting metabolites from 2 consortia
#> ✔ Collecting metabolites from 2 consortia [50ms]
#> 
#> ℹ Re-indexing 7 unique metabolites
#> ✔ Re-indexing 7 unique metabolites [43ms]
#> 
#> ℹ Expanding 2 binary matrices to 7-dimensional space
#> ✔ Expanding 2 binary matrices to 7-dimensional space [29ms]
#> 
#> ℹ Computing 7 x 7 levels matrix
#> ✔ Computing 7 x 7 levels matrix [32ms]
#> 
#> ℹ Computing pairwise overlap (1 pairs via crossprod)
#> ✔ Computing pairwise overlap (1 pairs via crossprod) [31ms]
#> 
#> ℹ Assembling pathway data from 2 consortia
#> ✔ Assembling pathway data from 2 consortia [39ms]
#> 
#> ℹ Building dendrogram from 2 x 2 dissimilarity matrix
#> ✔ Building dendrogram from 2 x 2 dissimilarity matrix [27ms]
#> 
#> ℹ Extracting dendrogram node positions
#> ✔ Extracting dendrogram node positions [30ms]
#> 
#> ℹ Collecting 2 consortium graphs
#> CMS "example" created: 2 consortia, 7 metabolites (0.3s)
#> ✔ Collecting 2 consortium graphs [107ms]
#> 
cms
#> 
#> ── ConsortiumMetabolismSet 
#> Name: "example"
#> 2 consortia, 7 species, 7 metabolites.
#> Community size (species): min 3, mean 3.5, max 4.
#> Community size (metabolites): min 2, mean 4, max 6.
```
