# Parse a single MiSoSoup solution entry

Extracts the species list from the `community` block, classifies each
reaction in the `solution` block as growth / species-scoped exchange /
media-level exchange, and returns the three components as a list.

## Usage

``` r
.parseMisosoupSolution(
  sol,
  normalize_ids = TRUE,
  biomassPattern = "(?i)^.*?(?:^|_)(?:biomass|growth)_"
)
```

## Arguments

- sol:

  A single solution entry with `community` and `solution` fields.

- normalize_ids:

  Whether to decode `__NN__` escapes and strip `R_EX_`/`_e` from
  metabolite IDs.

- biomassPattern:

  Perl regex identifying biomass rows; must consume the full prefix
  through the species-ID separator (used for both detection and species
  extraction via [`sub()`](https://rdrr.io/r/base/grep.html)).

## Value

A list with `consortia` (tibble of species, metabolite, flux for
species-scoped exchanges), `media` (tibble of metabolite, flux for
media-level bounds), `growth` (named numeric vector), and
`communityGrowth` (scalar `community_growth` value or `NA_real_`).

## Details

Species-scoped exchanges are matched by stripping `_<species>_i`
suffixes against the known species list (from `community`), which is
robust to both escaped metabolite names and species IDs containing
underscores (the old `_e_` delimiter split failed on both).

The biomass-reaction name is model-dependent (MiSoSoup carries it
through from the underlying GEM), so matching is done with a
configurable Perl regex that consumes the full prefix through the
species-ID separator. The default handles the three observed variants:
`Growth_<sp>`, `R_Biomass_<sp>`, and `R_R_BIOMASS_<sp>`.

A `community_growth` summary row (MiSoSoup's community-level growth
rate, distinct from per-species biomass fluxes) is captured as a scalar
in the returned list rather than treated as a media flux. Non-suffixed
`R_EX_<metabolite>` rows remain in media — they are legitimate
community-level exchange totals.
