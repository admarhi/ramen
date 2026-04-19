# Parse a single MiSoSoup solution entry

Extracts the species list from the `community` block, classifies each
reaction in the `solution` block as growth / species-scoped exchange /
media-level exchange, and returns the three components as a list.

## Usage

``` r
.parseMisosoupSolution(sol, normalize_ids = TRUE)
```

## Arguments

- sol:

  A single solution entry with `community` and `solution` fields.

- normalize_ids:

  Whether to decode `__NN__` escapes and strip `R_EX_`/`_e` from
  metabolite IDs.

## Value

A list with `consortia` (tibble of species, metabolite, flux for
species-scoped exchanges), `media` (tibble of metabolite, flux for
media-level bounds), and `growth` (named numeric vector).

## Details

Species-scoped exchanges are matched by stripping `_<species>_i`
suffixes against the known species list (from `community`), which is
robust to both escaped metabolite names and species IDs containing
underscores (the old `_e_` delimiter split failed on both).

The growth row key is matched against both `Growth_<species>` (older
MiSoSoup runs) and `R_Biomass_<species>` (newer runs).
