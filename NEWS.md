# ramen 0.99.0

Development release prepared for Bioconductor submission. This
release consolidates all pre-submission development since the initial
development snapshot (0.0.0.9001). As ramen has not had a released
version, none of the breaking changes below carries a deprecation
cycle.

## Breaking changes

* `importMisosoup()` now returns a single `ConsortiumMetabolismSet`.
  Earlier development versions returned a list of three tibbles; the
  importer now parses MiSoSoup YAML directly into a container of
  `ConsortiumMetabolism` objects, with media and growth carried as
  per-CM `metadata()`.
* `species()` on a `ConsortiumMetabolismSet` now returns a character
  vector of species identifiers. It previously returned a tibble.
  Use the new `speciesSummary()` for a per-species tabular
  breakdown.
* `ConsortiumMetabolismAlignment` is now a free-standing S4 class and
  no longer inherits from `TreeSummarizedExperiment`. Its slots are
  bespoke to alignment results. `ConsortiumMetabolism` and
  `ConsortiumMetabolismSet` continue to inherit from
  `TreeSummarizedExperiment`.
* `EffectiveConsumption` and `EffectiveProduction` assays now store
  the flux-corrected effective flux F * 2^H(p) (matching equation 2
  of the underlying thesis), not the bare Hill-1 perplexity 2^H(p).
  Same units as the `Consumption` / `Production` assays.
* Two new assays expose the previous quantity under explicit names:
  `nEffectiveSpeciesConsumption` and `nEffectiveSpeciesProduction`
  carry the Hill-1 effective number of contributing species
  (unitless, in [1, S]), mirroring the existing `nSpecies` count.
  Algebraic identity: `EffectiveConsumption = Consumption *
  nEffectiveSpeciesConsumption` (modulo two-decimal rounding).
* `consortia()` no longer has a `ConsortiumMetabolism` method. The
  plural noun applies only to containers of consortia, so the
  accessor is now scoped to `ConsortiumMetabolismSet` (returns the
  list of constituent CMs). For the underlying edge list of a single
  `ConsortiumMetabolism`, use `as.data.frame(cm)`. `consortia(cma)`
  raises an error by design -- a `ConsortiumMetabolismAlignment`
  records only its inputs' names, not copies of them.
* `plotDirectedFlow()` parameters renamed to lowerCamelCase
  throughout, with British-English spelling for colour-related
  arguments. The size-bearing arguments (`nodeSize`,
  `nodeLabelSize`, `edgeArrowSize`) now use ggraph millimetre units
  rather than igraph `cex` factors, with defaults adjusted for
  legibility. This is a hard rename -- no aliases for the old names.

## New features

* New importer `importSmetana()` reads raw SMETANA detailed output
  into a `ConsortiumMetabolismSet`.
* New `compareSpecies()` compares the metabolic role of species
  within a single CM and across CMs.
* New `speciesSummary()` returns a per-species tabular summary for a
  `ConsortiumMetabolismSet`.
* New `filterConsortia()` selects a subset of consortia from a
  `ConsortiumMetabolismSet`.
* `align()` gains a database-search mode: `align(CM, CMS)` searches a
  query consortium against a set and returns ranked matches with
  permutation p-values.
* `[` subsetting is now defined for `ConsortiumMetabolismSet`.
* New `growth()` accessor; `ConsortiumMetabolism()` gains a `growth`
  argument for attaching per-species growth rates.
* `functionalGroups()` now has a `ConsortiumMetabolism` method,
  complementing the existing `ConsortiumMetabolismSet` method.
* `metabolites()` on a `ConsortiumMetabolism` gains `species` and
  `direction` arguments for direction-aware metabolite queries.
* `as.data.frame()` is now defined for `ConsortiumMetabolism`
  (returns the edge list with columns `met`, `species`, `flux`) and
  `ConsortiumMetabolismSet` (row-binds per-CM edges and prefixes a
  `consortium` column).
* `importMisosoup()` gains a `biomassPattern` argument (Perl regex)
  so users can point it at arbitrary biomass-reaction names. It also
  records CMSC/MSC provenance on each returned CM as
  `metadata(cm)$misosoupMode` and `metadata(cm)$focalStrain`, and
  captures the `community_growth` summary row as
  `metadata(cm)$communityGrowth`.

## Visualisation

* Network plots migrated to ggraph and unified under a shared
  `theme_ramen()` and the exported `ramenPalette`. Plots use a white
  background, bottom legends, a colour-blind-safe node palette, and a
  blue edge-weight scale; non-data axes and grid lines are stripped.
* `show()` output enriched with per-species, per-class, and
  alignment-level detail across all three classes.

## Bug fixes

* `ConsortiumMetabolism()` and `ConsortiumMetabolismSet()` reject
  malformed input (non-finite flux, empty-string species or
  metabolite, NULL / empty / wrong-type arguments) with actionable
  `cli` errors instead of leaking R, S4, or dplyr internals.
* `new("ConsortiumMetabolism")`, `new("ConsortiumMetabolismSet")`,
  and other misuse paths now produce informative errors with hints.
* CM and CMS validity now accept an `NA` name.
* `species(CMS)` honours the `type` argument; `pathways(CMS)`
  niche / aux buckets now include ties at the quantile cutoff.
* `metabolites()` on a subset `ConsortiumMetabolismSet` no longer
  returns the parent union.
* `ConsortiumMetabolismAlignment` accessors no longer error on
  uninitialised objects.
* `plotFunctionalGroups()` guards against degenerate small
  dendrograms.
* Weighted-data detection in the `ConsortiumMetabolism()`
  constructor is more robust.
* Several fixes to `importMisosoup()` / `overviewMisosoup()` and a
  batch of fixes from the pre-submission usability audit.

## Documentation

* All user-facing messages, documentation, and comments standardised
  on British English spelling; `Language: en-GB` declared in
  DESCRIPTION.
* Vignettes synced with the 8-assay schema and expanded with a
  mathematical-formulation section for alignment, a TSE-interop
  section, a glossary, and species-name-normalisation guidance.
* Added `inst/CITATION` and `inst/extdata/misosoup/README.md`.

## Internal

* Minimum R version lowered to 4.4.0.
* Alignment method bodies split into named helpers; source formatted
  with air and linted clean.
* Added a test-coverage workflow (Codecov).

# ramen 0.0.0.9001

## New Features

* Initial development release.
* Three core S4 classes: `ConsortiumMetabolism`, `ConsortiumMetabolismSet`,
  and `ConsortiumMetabolismAlignment`, all inheriting from
  `TreeSummarizedExperiment`.
* Import microbial metabolic network data from MiSoSoup YAML format via
  `importMisosoup()`.
* Pairwise and multiple alignment of consortium metabolisms using functional
  overlap scores (FOS), Jaccard, Bray-Curtis, and redundancy overlap metrics.
* Visualisation methods including heatmaps, network plots, and alignment
  score plots.
* Functional group analysis with hierarchical clustering of species by
  shared metabolic pathways.
