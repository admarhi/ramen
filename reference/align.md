# Align a `ConsortiaMetabolismSet` Object

Creates an alignment of multiple consortium metabolisms to identify
common metabolic patterns and interactions across different communities.

## Usage

``` r
align(object, name)

# S4 method for class 'ConsortiumMetabolismSet'
align(object, name)
```

## Arguments

- object:

  A `ConsortiumMetabolismSet` object containing multiple consortium
  metabolisms to align

- name:

  Character scalar giving name of the alignment, if `NULL` inherits from
  the `ConsortiumMetabolismSet` object.

## Value

A `ConsortiumMetabolismAlignment` object containing the aligned
metabolic networks and associated metrics

## Methods (by class)

- `align(ConsortiumMetabolismSet)`: Align a `ConsortiumMetabolismSet`
  Object
