# ramen (development version)

## Breaking changes

* `consortia()` no longer has a `ConsortiumMetabolism` method. The
  plural noun applies only to containers of consortia, so the accessor
  is now scoped to `ConsortiumMetabolismSet` (returns the list of
  constituent CMs). For the underlying edge list of a single
  `ConsortiumMetabolism`, use `as.data.frame(cm)`.
  `consortia(cma)` continues to raise an error -- by design,
  `ConsortiumMetabolismAlignment` is a result object that records
  only its inputs' names rather than retaining copies of them.

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
* Visualization methods including heatmaps, network plots, and alignment
  score plots.
* Functional group analysis with hierarchical clustering of species by
  shared metabolic pathways.
