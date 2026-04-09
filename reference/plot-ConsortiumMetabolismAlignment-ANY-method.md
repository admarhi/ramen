# Plot a ConsortiumMetabolismAlignment object

Plot a ConsortiumMetabolismAlignment object

## Usage

``` r
# S4 method for class 'ConsortiumMetabolismAlignment,ANY'
plot(x, type = NULL)
```

## Arguments

- x:

  A `ConsortiumMetabolismAlignment` object.

- type:

  Character specifying the plot type: `"heatmap"`, `"network"`, or
  `"scores"`.

## Value

For `"heatmap"` and `"scores"`, a `ggplot` object. For `"network"`,
invisibly returns `NULL` (base igraph plot).

## Examples

``` r
cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
cma <- align(cm1, cm2)
plot(cma)
```
