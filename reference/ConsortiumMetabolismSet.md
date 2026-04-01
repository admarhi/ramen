# Set of `ConsortiumMetabolism` Objects

Creates a `ConsortiumMetabolismSet` combining multiple
`ConsortiumMetabolism` objects into a unified metabolite space. Computes
pairwise overlap scores and builds a dendrogram for clustering.

## Usage

``` r
ConsortiumMetabolismSet(..., name = NA_character_, desc = NA_character_)
```

## Arguments

- ...:

  Lists or individual `ConsortiumMetabolism` objects.

- name:

  Character scalar giving the name of the set.

- desc:

  Optional short description of the set.

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
[TreeSummarizedExperiment](https://rdrr.io/pkg/TreeSummarizedExperiment/man/TreeSummarizedExperiment-class.html)

## Examples

``` r
cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
cms <- ConsortiumMetabolismSet(cm1, cm2, name = "example")
#> 
#> ── Creating CMS "example" ──────────────────────────────────────────────────────
#> ℹ Validating 2 <ConsortiumMetabolism> objects
#> ✔ Validating 2 <ConsortiumMetabolism> objects [14ms]
#> 
#> ℹ Collecting metabolites from 2 consortia
#> ✔ Collecting metabolites from 2 consortia [42ms]
#> 
#> ℹ Re-indexing 6 unique metabolites
#> ✔ Re-indexing 6 unique metabolites [35ms]
#> 
#> ℹ Expanding 2 binary matrices to 6-dimensional space
#> ✔ Expanding 2 binary matrices to 6-dimensional space [23ms]
#> 
#> ℹ Computing 6 x 6 levels matrix
#> ✔ Computing 6 x 6 levels matrix [26ms]
#> 
#> ℹ Computing pairwise overlap (1 pairs via crossprod)
#> ✔ Computing pairwise overlap (1 pairs via crossprod) [27ms]
#> 
#> ℹ Assembling pathway data from 2 consortia
#> ✔ Assembling pathway data from 2 consortia [42ms]
#> 
#> ℹ Building dendrogram from 2 x 2 dissimilarity matrix
#> ✔ Building dendrogram from 2 x 2 dissimilarity matrix [24ms]
#> 
#> ℹ Extracting dendrogram node positions
#> ✔ Extracting dendrogram node positions [168ms]
#> 
#> ℹ Collecting 2 consortium graphs
#> ✔ CMS "example" created: 2 consortia, 6 metabolites (0.4s)
#> ℹ Collecting 2 consortium graphs
#> ✔ Collecting 2 consortium graphs [87ms]
#> 
cms
#> 
#> ── ConsortiumMetabolismSet 
#> Name: "example"
#> Containing 2 consortia.
#> Description: NA
```
