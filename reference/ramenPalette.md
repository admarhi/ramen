# Ramen package palette

Named list of colour vectors used by every plotting function in ramen.
Centralising the palette here ensures that every figure produced by the
package shares the same colour vocabulary. Where possible, colours are
drawn from the Okabe-Ito palette (colour-blind safe) and from
ColorBrewer sequential ramps.

## Usage

``` r
ramenPalette
```

## Format

A named `list` with elements:

- `nodeRole`:

  Three colours mapping `source`/`intermediate`/`sink`.

- `edgeCategorical`:

  Three colours mapping `shared`/`query`/`reference` edge categories in
  pairwise alignment networks.

- `edgeWeight`:

  Two-colour ramp (`low`/`high`) for continuous edge-weight gradients.

- `heatmapFill`:

  Three-stop ramp (`low`/`mid`/`high`) for similarity heatmaps.

- `bar`:

  Single colour for bar fills.

- `cluster`:

  Eight-colour Okabe-Ito-derived palette, recycled for k-cluster
  dendrogram colouring.

## Value

A named `list` of character vectors.

## Examples

``` r
ramenPalette$nodeRole
#>       source intermediate         sink 
#>    "#0072B2"    "#F0E442"    "#D55E00" 
ramenPalette$edgeCategorical
#>    shared     query reference 
#> "#000000" "#0072B2" "#D55E00" 
```
