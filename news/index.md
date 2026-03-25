# Changelog

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
