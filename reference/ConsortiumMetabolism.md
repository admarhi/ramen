# Functional Microbiome Representation based on TreeSummarizedExperiment

Creates a ConsortiumMetabolism object representing metabolic
interactions in a microbial community. The object contains information
about metabolite consumption and production by different species, along
with various metrics like flux sums and effective fluxes.

## Usage

``` r
ConsortiumMetabolism(
  data,
  name = NA_character_,
  species_col = "species",
  metabolite_col = "met",
  flux_col = "flux",
  ...
)
```

## Arguments

- data:

  a DataFrame-like object that includes columns specifying the species,
  metabolites and fluxes in the microbiome. The fluxes can either be
  weighted or unweighted (all of magnitude 1).

- name:

  a character scalar specifying the name of the Microbiome

- species_col:

  Character scalar specifying the name of the species column, defaults
  to 'species'.

- metabolite_col:

  Character scalar specifying the name of the metabolite column,
  defaults to 'met'.

- flux_col:

  Character scalar specifying the name of the flux column, defaults to
  'flux'.

- ...:

  Additional arguments to be passed to the constructor.

## Value

A ConsortiumMetabolism object containing:

- Assays for binary interactions, edge counts, consumption/production
  metrics

- Row and column data about metabolites

- Graph representation of the metabolic network

- Original input data and computed edge information
