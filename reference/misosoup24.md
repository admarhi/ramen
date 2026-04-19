# Example Data from MiSoSoup 2024

List of 56 solutions from MiSoSoup for different focal strains in
different media.

## Format

A list of 56 tibbles, each with columns:

- species:

  Character, species identifier.

- metabolite:

  Character, metabolite identifier following the BiGG namespace (King et
  al. 2016). Common examples: `ac` (acetate), `co2` (CO2), `etoh`
  (ethanol), `pyr` (pyruvate), `ala__L` (L-alanine), `ala__D`
  (D-alanine). The double-underscore convention encodes stereochemistry
  (e.g. `__L` = L-isomer, `__D` = D-isomer). See <https://bigg.ucsd.edu>
  for the full namespace.

- flux:

  Numeric, metabolic flux (negative for consumption, positive for
  production).

## Value

A list of 56 tibbles.

## Details

Metabolite identifiers follow the BiGG (Biochemically, Genetically and
Genomically structured) knowledge base namespace
(<https://bigg.ucsd.edu>; King et al. 2016). Raw exchange reaction IDs
from flux-balance simulators (e.g. `R_EX_ac_e`, `EX_ac_e`) are
normalized to bare metabolite names by
[`importMisosoup()`](https://admarhi.github.io/ramen/reference/importMisosoup.md)
via `normalize_ids = TRUE` (default).

## References

King ZA et al. (2016) BiGG Models: A platform for integrating,
standardizing and sharing genome-scale models. *Nucleic Acids Research*,
44(D1), D515–D522.
[doi:10.1093/nar/gkv1049](https://doi.org/10.1093/nar/gkv1049)

## Examples

``` r
data("misosoup24")
head(misosoup24[[1]])
#> # A tibble: 6 × 3
#>   metabolite species    flux
#>   <chr>      <chr>     <dbl>
#> 1 ac         A1R12     0.773
#> 2 ac         I2R16   -10.8  
#> 3 acald      A1R12    -1.12 
#> 4 acald      I2R16     1.12 
#> 5 ala__D     A1R12     0.760
#> 6 ala__D     I2R16    -0.760
```
