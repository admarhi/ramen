# Constructor for `ConsortiumMetabolismAlignment` Objects

Builds a valid `ConsortiumMetabolismAlignment` by constructing a TSE
base and setting CMA-specific slots with sensible defaults.

## Usage

``` r
ConsortiumMetabolismAlignment(..., tse = NULL)
```

## Arguments

- ...:

  Named CMA slot values to set (e.g., `Type = "pairwise"`,
  `PrimaryScore = 0.8`).

- tse:

  An optional
  [`TreeSummarizedExperiment`](https://rdrr.io/pkg/TreeSummarizedExperiment/man/TreeSummarizedExperiment-constructor.html)
  to use as the base. Default creates an empty TSE.

## Value

A validated `ConsortiumMetabolismAlignment` object.

## Slots

- `Name`:

  character. Display name for the alignment.

- `Description`:

  character. Optional short description.

- `Type`:

  character. One of `"pairwise"`, `"multiple"`, or `"search"`.

- `Metric`:

  character. Similarity metric used (e.g., `"FOS"`).

- `Params`:

  list. Additional parameters passed to the alignment method.

- `QueryName`:

  character. Name of the query consortium (pairwise).

- `ReferenceName`:

  character. Name of the reference consortium (pairwise).

- `Scores`:

  list. Named list of all computed score components.

- `PrimaryScore`:

  numeric. Primary similarity score, between 0 and 1.

- `Pvalue`:

  numeric. Permutation p-value, between 0 and 1.

- `SharedPathways`:

  data.frame. Pathways present in both query and reference.

- `UniqueQuery`:

  data.frame. Pathways unique to the query.

- `UniqueReference`:

  data.frame. Pathways unique to the reference.

- `SimilarityMatrix`:

  matrix. Pairwise similarity matrix (multiple alignment).

- `ConsensusPathways`:

  data.frame. Consensus network pathways (multiple alignment).

- `Prevalence`:

  data.frame. Pathway prevalence across consortia (multiple alignment).

- `Dendrogram`:

  list. Hierarchical clustering dendrogram (multiple alignment).

- `Pathways`:

  data.frame. Combined pathway list.

- `Graphs`:

  list. List of igraph objects.

- `Metabolites`:

  data.frame. Metabolite mapping table.

## See also

[TreeSummarizedExperiment](https://rdrr.io/pkg/TreeSummarizedExperiment/man/TreeSummarizedExperiment-class.html),
[`align`](https://admarhi.github.io/ramen/reference/align.md)

## Examples

``` r
# Empty alignment
cma <- ConsortiumMetabolismAlignment()

# Pairwise alignment via align()
cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
cma <- align(cm1, cm2)
cma
#> 
#> ── ConsortiumMetabolismAlignment 
#> Name: NA
#> Type: "pairwise"
#> Metric: "FOS"
#> Score: 0
#> Query: "comm_1", Reference: "comm_2"
```
