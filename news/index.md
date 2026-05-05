# Changelog

## ramen (development version)

- [`importMisosoup()`](https://admarhi.github.io/ramen/reference/importMisosoup.md)
  gains a `biomassPattern` argument (Perl regex) so users can point it
  at arbitrary biomass-reaction names. The default matches
  `Growth_<sp>`, `R_Biomass_<sp>`, and `R_R_BIOMASS_<sp>`
  case-insensitively (MiSoSoup’s biomass reaction name is
  model-dependent).
- [`importMisosoup()`](https://admarhi.github.io/ramen/reference/importMisosoup.md)
  now records CMSC/MSC provenance on each returned CM as
  `metadata(cm)$misosoupMode` (`"CMSC"` or `"MSC"`) and
  `metadata(cm)$focalStrain` (character or `NA`). The `community_growth`
  summary row is captured as `metadata(cm)$communityGrowth` instead of
  leaking into the media bucket.

## ramen 0.0.0.9001

### New Features

- Initial development release.
- Three core S4 classes: `ConsortiumMetabolism`,
  `ConsortiumMetabolismSet`, and `ConsortiumMetabolismAlignment`, all
  inheriting from `TreeSummarizedExperiment`.
- Import microbial metabolic network data from MiSoSoup YAML format via
  [`importMisosoup()`](https://admarhi.github.io/ramen/reference/importMisosoup.md).
- Pairwise and multiple alignment of consortium metabolisms using
  functional overlap scores (FOS), Jaccard, Bray-Curtis, and redundancy
  overlap metrics.
- Visualization methods including heatmaps, network plots, and alignment
  score plots.
- Functional group analysis with hierarchical clustering of species by
  shared metabolic pathways.
