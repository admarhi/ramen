# Compare Two Species by Pathway Set

Computes similarity metrics between the pathway sets of two species. Two
dispatch modes are supported:

- `compareSpecies(cm, sp1, sp2)`: compare two species within the same
  `ConsortiumMetabolism` object.

- `compareSpecies(cm1, cm2, sp1, sp2)`: compare one species from each of
  two `ConsortiumMetabolism` objects (e.g. the same species under
  different growth conditions).

The pathway set for a species is the set of (consumed, produced) pairs
in which that species participates.

## Usage

``` r
compareSpecies(x, y, ...)

# S4 method for class 'ConsortiumMetabolism,character'
compareSpecies(x, y, sp2, ...)

# S4 method for class 'ConsortiumMetabolism,ConsortiumMetabolism'
compareSpecies(x, y, sp1, sp2, ...)
```

## Arguments

- x:

  A
  [ConsortiumMetabolism](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md)
  object (first).

- y:

  A
  [ConsortiumMetabolism](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md)
  object (second).

- ...:

  For same-CM: `sp2` (character scalar naming the second species). For
  cross-CM: `sp1` and `sp2` (character scalars naming one species from
  each consortium).

- sp2:

  Character scalar; species from `y`.

- sp1:

  Character scalar; species from `x`.

## Value

A named list with elements:

- `fos`:

  Szymkiewicz-Simpson overlap score (intersection over min set size).

- `jaccard`:

  Jaccard similarity (intersection over union).

- `n_shared`:

  Number of shared pathways.

- `n_unique_sp1`:

  Pathways only in sp1.

- `n_unique_sp2`:

  Pathways only in sp2.

## Methods (by class)

- `compareSpecies(x = ConsortiumMetabolism, y = character)`: Compare two
  species within the same
  [ConsortiumMetabolism](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md)

- `compareSpecies(x = ConsortiumMetabolism, y = ConsortiumMetabolism)`:
  Compare one species from each of two
  [ConsortiumMetabolism](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md)
  objects

## Examples

``` r
cm <- synCM("test", n_species = 3, max_met = 5)
sp <- species(cm)
compareSpecies(cm, sp[1], sp[2])
#> $fos
#> [1] 0.3333333
#> 
#> $jaccard
#> [1] 0.125
#> 
#> $n_shared
#> [1] 1
#> 
#> $n_unique_sp1
#> [1] 2
#> 
#> $n_unique_sp2
#> [1] 5
#> 
```
