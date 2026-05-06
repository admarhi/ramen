# ramen (development version)

## Breaking changes

* `EffectiveConsumption` and `EffectiveProduction` assays now store the
  flux-corrected effective flux $F \cdot 2^{H(p)}$ (matching equation 2
  of the underlying thesis), not the bare Hill-1 perplexity $2^{H(p)}$.
  Same units as the `Consumption` / `Production` assays. Pre-Bioconductor
  release, no users yet, so no deprecation cycle.
* Two new assays expose the previous quantity under explicit names:
  `nEffectiveSpeciesConsumption` and `nEffectiveSpeciesProduction` carry
  the Hill-1 effective number of contributing species (unitless,
  $\in [1, S]$), mirroring the existing `nSpecies` count.
  Algebraic identity: `EffectiveConsumption = Consumption *
  nEffectiveSpeciesConsumption` (modulo two-decimal rounding).
* `consortia()` no longer has a `ConsortiumMetabolism` method. The
  plural noun applies only to containers of consortia, so the accessor
  is now scoped to `ConsortiumMetabolismSet` (returns the list of
  constituent CMs). For the underlying edge list of a single
  `ConsortiumMetabolism`, use `as.data.frame(cm)`.
  `consortia(cma)` continues to raise an error -- by design,
  `ConsortiumMetabolismAlignment` is a result object that records
  only its inputs' names rather than retaining copies of them.
* `plotDirectedFlow()` parameters renamed to lowerCamelCase
  throughout, with British-English spelling for colour-related
  args. The size-bearing arguments (`nodeSize`, `nodeLabelSize`,
  `edgeArrowSize`) now use ggraph millimetre units rather than
  igraph `cex` factors, with defaults adjusted for legibility.
  This is a hard rename -- no aliases for the old names.

## New methods

* `as.data.frame()` is now defined for `ConsortiumMetabolism` (returns
  the edge list with columns `met`, `species`, `flux`) and
  `ConsortiumMetabolismSet` (row-binds per-CM edges and prefixes a
  `consortium` column). The pre-existing
  `ConsortiumMetabolismAlignment` method is unchanged.

## Other changes

* `importMisosoup()` gains a `biomassPattern` argument (Perl regex) so
  users can point it at arbitrary biomass-reaction names. The default
  matches `Growth_<sp>`, `R_Biomass_<sp>`, and `R_R_BIOMASS_<sp>`
  case-insensitively (MiSoSoup's biomass reaction name is
  model-dependent).
* `importMisosoup()` now records CMSC/MSC provenance on each returned
  CM as `metadata(cm)$misosoupMode` (`"CMSC"` or `"MSC"`) and
  `metadata(cm)$focalStrain` (character or `NA`). The `community_growth`
  summary row is captured as `metadata(cm)$communityGrowth` instead of
  leaking into the media bucket.
* Internal documentation, comments, and user-facing CLI messages
  standardised on British English spelling (`colour`, `centre`,
  `organise`, etc.). External-API identifiers from upstream packages
  (igraph, ggplot2, ColorBrewer) keep their original spelling. The
  `Language: en-GB` DESCRIPTION field will follow in a separate commit.

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
