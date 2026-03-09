# Get Functional Groups

Calculates and returns functional groups based on metabolic reactions.
For `ConsortiumMetabolismSet` objects, this involves analyzing shared
reactions across species to identify clusters of species with similar
metabolic capabilities.

## Usage

``` r
getFunctionalGroups(object, k = 4, ...)

# S4 method for class 'ConsortiumMetabolismSet'
getFunctionalGroups(object, k = 4, label_size = 6, label_colours = NULL)
```

## Arguments

- object:

  A `ConsortiumMetabolismSet` object.

- k:

  An integer scalar specifying the number of clusters to color in the
  dendrogram.

- ...:

  Additional arguments to be passed to specific methods.

- label_size:

  Numeric scalar specifying label size in the output plot.

- label_colours:

  If not `NULL`, tb with two columns, label and colour.

## Value

A dendrogram object representing the hierarchical clustering of species
into functional groups. The plot of the dendrogram is also displayed.

## Details

This method is currently implemented for `ConsortiumMetabolismSet`
objects. Future versions will extend functionality to
`ConsortiumMetabolismAlignment` objects to allow for comparative
functional group analysis across different alignments.
