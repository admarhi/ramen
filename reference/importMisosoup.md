# Import and Process Misosoup Data

Processes raw misosoup data into a structured format containing
consortia, media, and growth information. The function handles data
cleaning, transformation, and organization of metabolic flux data.

## Usage

``` r
importMisosoup(data)
```

## Arguments

- data:

  A nested list containing misosoup simulation results. The structure
  should be data\[\[substrate\]\]\[\[focal_strain\]\] where each element
  contains solution data.

## Value

A list containing three tibbles:

- consortia:

  Metabolic flux data for each species in the consortium

- media:

  Media composition data

- growth:

  Growth information for each solution

## See also

[`overviewMisosoup`](https://admarhi.github.io/ramen/reference/overviewMisosoup.md)
for a summary of the input data

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires raw MiSoSoup YAML data
raw <- yaml::read_yaml("misosoup_output.yaml")
result <- importMisosoup(raw)
str(result, max.level = 1)
} # }
```
