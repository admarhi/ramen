# Plot a ConsortiumMetabolism object

Plot a ConsortiumMetabolism object

## Usage

``` r
# S4 method for class 'ConsortiumMetabolism,ANY'
plot(
  x,
  type = c("Binary", "nSpecies", "Consumption", "Production", "EffectiveConsumption",
    "EffectiveProduction", "nEffectiveSpeciesConsumption", "nEffectiveSpeciesProduction")
)
```

## Arguments

- x:

  A `ConsortiumMetabolism` object.

- type:

  Character specifying the assay to plot.

## Value

A `ggplot` object rendered with ggraph.

## Examples

``` r
cm <- synCM("test", n_species = 3, max_met = 5)
plot(cm)
```
