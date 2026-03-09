# Get Edges From a `ConsortiumMetabolism` Object

## Usage

``` r
getEdges(object, ...)

# S4 method for class 'ConsortiumMetabolism'
getEdges(object)

# S4 method for class 'ConsortiumMetabolismSet'
getEdges(
  object,
  type = c("all", "pan-cons", "niche", "core", "aux"),
  quantileCutoff = 0.1
)
```

## Arguments

- object:

  A `ConsortiumMetabolism` object

- ...:

  Object specific arguments.

- type:

  Character scalar giving the type of edges to output.

- quantileCutoff:

  Numeric scalar between 0 and 1 giving the quantile threshold to use
  for filtering edges. For "pan-cons" and "core" types, edges above
  `1 - quantileCutoff` are returned. For "niche" and "aux" types, edges
  below `quantileCutoff` are returned. Defaults to 0.1 (i.e., top/bottom
  10\\

A tibble containing edge information including:

- consumed/produced metabolites

- number of species involved

- consumption/production sums

- effective consumption/production metrics

Retrieves the edges representing metabolic interactions between species.
The argument `type` can be used to return only specfic types of Edges
from a `ConsortiumMetabolismSet` type object.

- `all` will return all edges

- `pan-cons` will return only edges that exist in all consortia

- `niche` will return niche consortia. A niche is defined a

Methods (by class)

- `getEdges(ConsortiumMetabolism)`: Get Edges From a
  `ConsortiumMetabolism` Object

- `getEdges(ConsortiumMetabolismSet)`: Get Edges From a
  `ConsortiumMetabolismSet` Object
