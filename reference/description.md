# Get or Set Object Description

Returns or sets the description of a ramen object.

## Usage

``` r
description(object)

description(object) <- value

# S4 method for class 'ConsortiumMetabolism'
description(object)

# S4 method for class 'ConsortiumMetabolismSet'
description(object)

# S4 method for class 'ConsortiumMetabolismAlignment'
description(object)

# S4 method for class 'ConsortiumMetabolism'
description(object) <- value

# S4 method for class 'ConsortiumMetabolismSet'
description(object) <- value

# S4 method for class 'ConsortiumMetabolismAlignment'
description(object) <- value
```

## Arguments

- object:

  A
  [ConsortiumMetabolism](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md),
  [ConsortiumMetabolismSet](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismSet.md),
  or
  [ConsortiumMetabolismAlignment](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismAlignment.md)
  object.

- value:

  Character scalar specifying the new description.

## Value

A character scalar (getter), or the modified object (setter).

## Examples

``` r
cm <- synCM("test", n_species = 3, max_met = 5)
description(cm)
#> character(0)
```
