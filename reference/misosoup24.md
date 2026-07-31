# Example Data from MiSoSoup 2024

Pre-imported exchange-flux tables for 56 microbial-community solutions
enumerated by MiSoSoup. Each list element captures a single alternative
optimal community around a focal strain on a single carbon source, in
the flat `metabolite` / `species` / `flux` long form that
[`ConsortiumMetabolism()`](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md)
consumes directly.

## Format

A named list of 56 data.frames, each with columns:

- species:

  Character, species identifier. The species flagged as *focal* appears
  in the list element name (see naming convention in *Details*); the
  remaining species are the *support* community selected by the MILP
  solver to keep the focal strain viable.

- metabolite:

  Character, metabolite identifier following the BiGG namespace (King et
  al. 2016). Common examples: `ac` (acetate), `co2` (CO2), `etoh`
  (ethanol), `pyr` (pyruvate), `ala__L` (L-alanine), `ala__D`
  (D-alanine). The double-underscore convention encodes stereochemistry
  (e.g. `__L` = L-isomer, `__D` = D-isomer). See <https://bigg.ucsd.edu>
  for the full namespace.

- flux:

  Numeric, metabolic flux. Negative values denote consumption (uptake by
  the species), positive values denote production (secretion).

## Value

A named list of 56 data.frames.

## Details

**What MiSoSoup is.** MiSoSoup is a Mixed-Integer Linear Programming
(MILP) enumerator that, given a community-level metabolic-modelling
objective, returns *multiple alternative optimal community compositions*
satisfying that objective. The output ramen consumes is the per-solution
exchange-flux table: which species exchange which metabolites, and at
what flux, in each viable community.

**Focal strain and support community.** Each MiSoSoup run is anchored on
a *focal strain* – the species the solver is instructed to keep alive
while exploring viable communities around it. The other species in each
solution are the *support*: the MILP picks them so that the focal strain
achieves the requested objective. A single focal-strain / substrate
combination typically yields many alternative optima with distinct
support compositions; these become the separate list elements that share
a substrate / strain prefix.

**Biological scope.** The 56 solutions span three focal strains
(`A1R12`, `A3R04`, `B3M02`) on three carbon sources: acetate (`ac`),
citrate (`cit`), and fructose-6-phosphate (`f6p`). The growth medium is
a minimal medium with the named substrate as the sole carbon source (see
source paper for the full medium composition).

**Naming convention.** List element names follow the pattern
`{substrate}_{focal-strain}_{solution-id}`. For example, `ac_A1R12_1` is
alternative optimum 1 for focal strain `A1R12` on acetate. The substrate
/ strain grouping can be recovered by splitting the names on
underscores.

**Caution – alternative optima.** Two solutions with the same substrate
/ strain prefix are *alternative optima of the same MILP*, not
biologically independent communities. Their pairwise FOS overlap is
expected to be high by construction; treat such overlaps as
solver-consistency diagnostics rather than ecological signal.

**Relationship to
[`importMisosoup()`](https://admarhi.github.io/ramen/reference/importMisosoup.md).**
For raw MiSoSoup YAML output, use
[`importMisosoup()`](https://admarhi.github.io/ramen/reference/importMisosoup.md),
which parses the nested format and attaches per-solution biomass and
medium metadata to the resulting CMs. The `misosoup24` list is already
flattened to the long form that
[`ConsortiumMetabolism()`](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md)
accepts directly, so medium composition is implicit rather than stored
on each element.

**Metabolite identifiers.** BiGG (Biochemically, Genetically and
Genomically structured) knowledge base identifiers
(<https://bigg.ucsd.edu>; King et al. 2016). Raw exchange reaction IDs
from flux-balance simulators (e.g. `R_EX_ac_e`, `EX_ac_e`) are
normalised to bare metabolite names by
[`importMisosoup()`](https://admarhi.github.io/ramen/reference/importMisosoup.md)
via `normalize_ids = TRUE` (default).

## References

King ZA et al. (2016) BiGG Models: A platform for integrating,
standardizing and sharing genome-scale models. *Nucleic Acids Research*,
44(D1), D515–D522.
[doi:10.1093/nar/gkv1049](https://doi.org/10.1093/nar/gkv1049)

## See also

[`importMisosoup`](https://admarhi.github.io/ramen/reference/importMisosoup.md),
[`ConsortiumMetabolism`](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md)

## Examples

``` r
data("misosoup24")
length(misosoup24)
#> [1] 56
head(names(misosoup24))
#> [1] "ac_A1R12_1"  "ac_A1R12_10" "ac_A1R12_11" "ac_A1R12_12" "ac_A1R12_13"
#> [6] "ac_A1R12_14"
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
