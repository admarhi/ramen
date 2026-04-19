# Subset ramen Objects by Metabolite

Subsetting methods for ramen S4 classes. For `ConsortiumMetabolismSet`,
`[` follows the `TreeSummarizedExperiment` contract: the assays are
*metabolite x metabolite* (m x m), so `i` and `j` index the metabolite
space, not the consortium list. After subsetting, all custom slots are
updated eagerly: `BinaryMatrices` is subsetted to the remaining
metabolite space and `OverlapMatrix`, `Dendrogram`, `NodeData`,
`Pathways`, `Metabolites`, and `Graphs` are recomputed accordingly.
`@Consortia` is not modified — each `ConsortiumMetabolism` object
retains its own metabolite space.

**Note on consortium-level selection:** to extract a subset of consortia
rather than a subset of the metabolite space, use
[`extractCluster`](https://admarhi.github.io/ramen/reference/extractCluster.md)
(dendrogram-based) or the planned
[`filterConsortia()`](https://admarhi.github.io/ramen/reference/filterConsortia.md)
method. Metabolite subsetting is appropriate when the analysis should
focus on a specific part of the metabolic network — for example,
subsetting to the core-pathway metabolites
(`pathways(cms, type = "core")`) to compare consortia on their shared
metabolic backbone only.

`ConsortiumMetabolism` and `ConsortiumMetabolismAlignment` do not
support `[` and return an informative error.

## Usage

``` r
# S4 method for class 'ConsortiumMetabolism,ANY,ANY,ANY'
x[i, j, ..., drop = FALSE]

# S4 method for class 'ConsortiumMetabolismSet,ANY,ANY,ANY'
x[i, j, ..., drop = FALSE]

# S4 method for class 'ConsortiumMetabolismAlignment,ANY,ANY,ANY'
x[i, j, ..., drop = FALSE]
```

## Arguments

- x:

  A ramen S4 object.

- i:

  Row (metabolite) indices — integer, logical, or character names from
  `metabolites(cms)`.

- j:

  Column indices. Must be identical to `i` (or omitted): since assays
  are m x m, asymmetric subsetting is not meaningful.

- ...:

  Additional arguments passed to the TSE parent method.

- drop:

  Ignored; retained for S4 signature compatibility.

## Value

For `ConsortiumMetabolismSet`: a subsetted object with all custom slots
updated to reflect the remaining metabolite space. For
`ConsortiumMetabolism` and `ConsortiumMetabolismAlignment`: an error.

## Functions

- `x[i`: Subsetting a `ConsortiumMetabolism` is not supported.

- `x[i`: Subset a `ConsortiumMetabolismSet` by metabolite index with
  full slot synchronisation.

- `x[i`: Subsetting a `ConsortiumMetabolismAlignment` is not supported.

## See also

[`extractCluster`](https://admarhi.github.io/ramen/reference/extractCluster.md)
for dendrogram-based consortium selection.

## Examples

``` r
cm1 <- synCM("a", n_species = 3, max_met = 6, seed = 1)
cm2 <- synCM("b", n_species = 3, max_met = 6, seed = 2)
cms <- ConsortiumMetabolismSet(
    cm1, cm2, name = "test", verbose = FALSE
)
## Subset to first 3 metabolites (m x m assay semantics)
sub <- cms[seq_len(3), seq_len(3)]
nrow(sub@Metabolites)
#> [1] 3
```
