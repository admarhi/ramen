# Overview of Misosoup Data Structure

Provides a summary of the misosoup data structure, including the number
of consortia and zero-growth solutions for each substrate and focal
strain combination.

## Usage

``` r
overviewMisosoup(data)
```

## Arguments

- data:

  A nested list containing misosoup simulation results. The structure
  should be data\[\[substrate\]\]\[\[focal_strain\]\] where each element
  contains solution data.

## Value

A tibble with columns:

- substrate:

  Substrate identifier

- focal_strain:

  Focal strain identifier

- n_cons:

  Number of consortia solutions

- n_zero_growth:

  Number of zero-growth solutions

## See also

[`importMisosoup`](https://admarhi.github.io/ramen/reference/importMisosoup.md)
for processing the full data

## Examples

``` r
# \donttest{
# Requires raw MiSoSoup YAML data
raw <- yaml::read_yaml("misosoup_output.yaml")
#> Warning: cannot open file 'misosoup_output.yaml': No such file or directory
#> Error in file(file, "rt", encoding = fileEncoding): cannot open the connection
overviewMisosoup(raw)
#> # A tibble: 0 × 4
#> # ℹ 4 variables: substrate <chr>, focal_strain <lgl>, n_cons <int>,
#> #   n_zero_growth <dbl>
# }
```
