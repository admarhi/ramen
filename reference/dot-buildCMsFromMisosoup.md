# Build CMs from a raw MiSoSoup list

Iterates over `raw[[substrate]][[sec]][[sol]]`, parses each viable
solution, and constructs a
[ConsortiumMetabolism](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md)
object. Zero-growth solutions (empty community) are silently skipped.

## Usage

``` r
.buildCMsFromMisosoup(
  raw,
  normalize_ids = TRUE,
  biomassPattern = "(?i)^.*?(?:^|_)(?:biomass|growth)_"
)
```

## Arguments

- raw:

  Nested list from
  [`yaml::read_yaml()`](https://yaml.r-lib.org/reference/read_yaml.html).

- normalize_ids:

  Whether to normalize metabolite IDs.

- biomassPattern:

  Perl regex identifying biomass reaction rows.

## Value

Named list of CM objects.
