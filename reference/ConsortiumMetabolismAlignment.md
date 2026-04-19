# Constructor for `ConsortiumMetabolismAlignment` Objects

Builds a valid `ConsortiumMetabolismAlignment` by setting CMA-specific
slots with sensible defaults.

Exports alignment results to a plain `data.frame` suitable for
downstream analysis. For `"pairwise"` alignments the three pathway sets
(`SharedPathways`, `UniqueQuery`, `UniqueReference`) are row-bound with
a `pathway_type` column added. For `"multiple"` alignments
`ConsensusPathways` is returned. For all other types (or empty objects)
`Pathways` is returned, falling back to an empty
[`data.frame()`](https://rdrr.io/r/base/data.frame.html).

## Usage

``` r
ConsortiumMetabolismAlignment(...)

# S4 method for class 'ConsortiumMetabolismAlignment'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- ...:

  Additional arguments (currently unused).

- x:

  A `ConsortiumMetabolismAlignment` object.

- row.names:

  `NULL` or a character vector of row names; ignored.

- optional:

  Logical; ignored.

## Value

A validated `ConsortiumMetabolismAlignment` object.

A `data.frame`. For pairwise alignments the result contains at minimum
`consumed`, `produced`, and `pathway_type` columns.

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
#> Name: "comm_1 vs comm_2"
#> Type: "pairwise"
#> Metric: "FOS"
#> Score: 0
#> Query: "comm_1", Reference: "comm_2"
#> Coverage: query 0, reference 0

cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
cma <- align(cm1, cm2)
df <- as.data.frame(cma)
head(df)
#>   consumed produced querySpecies referenceSpecies     pathway_type
#> 1     met5     met4 TFL8958R....     JGE334P,....           shared
#> 2     met1     met5 DOR2356X....     JYF8221A....           shared
#> 3    media     met5                                   unique_query
#> 4     met4     met5                                   unique_query
#> 5     met5     met1                               unique_reference
#> 6     met6     met1                               unique_reference
```
