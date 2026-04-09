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
# \donttest{
# Requires raw MiSoSoup YAML data
raw <- yaml::read_yaml("misosoup_output.yaml")
#> Warning: cannot open file 'misosoup_output.yaml': No such file or directory
#> Error in file(file, "rt", encoding = fileEncoding): cannot open the connection
result <- importMisosoup(raw)
#> Error in dplyr::relocate(dplyr::bind_rows(Map(function(x, y) {    entries <- data[[x]][[y]]    dplyr::bind_rows(Map(function(z, idx) {        dplyr::mutate(tidyr::pivot_longer(tibble::as_tibble(z$solution),             cols = dplyr::everything(), names_to = "rxn", values_to = "flux"),             substrate = x, focal_strain = y, solution = idx)    }, entries, names(entries)))}, tb_import$substrate, tb_import$focal_strain)), "rxn", "flux",     .after = "solution"): Can't select columns that don't exist.
#> ✖ Column `rxn` doesn't exist.
str(result, max.level = 1)
#> Error: object 'result' not found
# }
```
