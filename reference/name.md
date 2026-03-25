# Get or Set Object Name

Returns or sets the name of a ramen object.

## Usage

``` r
name(object)

name(object) <- value

# S4 method for class 'ConsortiumMetabolism'
name(object)

# S4 method for class 'ConsortiumMetabolismSet'
name(object)

# S4 method for class 'ConsortiumMetabolismAlignment'
name(object)

# S4 method for class 'ConsortiumMetabolism'
name(object) <- value

# S4 method for class 'ConsortiumMetabolismSet'
name(object) <- value

# S4 method for class 'ConsortiumMetabolismAlignment'
name(object) <- value
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

  Character scalar specifying the new name.

## Value

A character scalar (getter), or the modified object (setter).

## Examples

``` r
cm <- synCM("test", n_species = 3, max_met = 5)
name(cm)
#> [1] "test"
```
