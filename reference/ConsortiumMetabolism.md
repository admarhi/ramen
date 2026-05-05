# Functional Microbiome Representation

Creates a `ConsortiumMetabolism` object representing metabolic
interactions in a microbial community. The object stores metabolite
consumption and production by different species, along with flux sums
and effective fluxes.

## Usage

``` r
ConsortiumMetabolism(
  data,
  name = NA_character_,
  growth = NULL,
  species_col = "species",
  metabolite_col = "metabolite",
  flux_col = "flux",
  verbose = FALSE
)
```

## Arguments

- data:

  a data.frame with columns for species, metabolites and fluxes. Fluxes
  can be weighted or unweighted (magnitude 1).

- name:

  Character scalar giving the consortium name.

- growth:

  Optional named numeric vector of per-species growth rates (e.g. FBA
  objective values). Names must match species in `data`. Stored in
  `metadata(cm)$growth` and retrievable via the
  [`growth`](https://admarhi.github.io/ramen/reference/growth.md)
  accessor.

- species_col:

  Character scalar for the species column name, defaults to `"species"`.

- metabolite_col:

  Character scalar for the metabolite column name, defaults to
  `"metabolite"`.

- flux_col:

  Character scalar for the flux column name, defaults to `"flux"`.

- verbose:

  Logical scalar. If `TRUE`, prints progress messages during
  construction. Defaults to `FALSE`.

## Value

A `ConsortiumMetabolism` object.

## Slots

- `Name`:

  character. Display name for the consortium.

- `Description`:

  character. Optional short description.

- `Pathways`:

  data.frame. Pathway list of metabolic interactions with per-pathway
  metrics (species, flux sums, effective diversity).

- `Weighted`:

  logical. Whether flux magnitudes are used.

- `InputData`:

  data.frame. Original input data (species, metabolite, flux columns).

- `Metabolites`:

  character. Unique metabolite identifiers.

- `Graphs`:

  list. List containing an igraph object of the metabolic network.

## See also

[TreeSummarizedExperiment-class](https://rdrr.io/pkg/TreeSummarizedExperiment/man/TreeSummarizedExperiment-class.html)

## Examples

``` r
cm <- synCM("example", n_species = 3, max_met = 5)
cm
#> 
#> ── ConsortiumMetabolism 
#> Name: "example"
#> Weighted metabolic network: 3 species, 3 metabolites, 4 pathways.
```
