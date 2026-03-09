# Compare Species

Compares the metabolisms of two species and outputs a list of tibbles,
containing the tibbles intersection, unique, and consistency. As the
name suggests, intersection, unique contain the intersection and unique
pathways per species compared, while the consistency tibble contains
information on whether or not a specie's set of pathways is consistent
in all of the consortia in which it is present or not. If all species
contain the same edges in all consortia in which they appear, this
tibble will be returned with 0 rows. For `ConsortiumMetabolismSet`
objects.

## Usage

``` r
compareSpecies(object, species)

# S4 method for class 'ConsortiumMetabolismSet'
compareSpecies(object, species)
```

## Arguments

- object:

  A `ConsortiumMetabolismSet` object.

- species:

  A character vector of species names to compare

## Value

A list of tibbles.

## Details

This method is currently implemented for `ConsortiumMetabolismSet`
objects. Future versions will extend functionality to
`ConsortiumMetabolism` objects to allow for the analysis of species
within a single consortium different alignments.

## Methods (by class)

- `compareSpecies(ConsortiumMetabolismSet)`: Compare Species in a
  ConsortiumMetabolismSet
