# MiSoSoup example data

This directory holds the raw exchange-flux tables behind the
`misosoup24` example dataset shipped with the ramen package
(`data("misosoup24")`).

## Contents

56 CSV files, one per microbial-community solution. Each file is a
long-format exchange-flux table with three columns:

| Column       | Type      | Description |
|--------------|-----------|-------------|
| `metabolite` | character | Metabolite identifier in the BiGG namespace. |
| `species`    | character | Species / strain identifier. |
| `flux`       | numeric   | Exchange flux. Negative = consumption (uptake), positive = production (secretion). |

File names follow the pattern
`{substrate}_{focal-strain}_{solution-id}.csv` -- for example,
`ac_A1R12_1.csv` is alternative optimum 1 for focal strain `A1R12`
on acetate.

## Provenance

These tables are a fixed, curated subset of a December 2024 MiSoSoup
enumeration run. MiSoSoup is a Mixed-Integer Linear Programming
enumerator that returns multiple alternative optimal community
compositions for a community-level metabolic-modelling objective.

- **Upstream source.** The full run is published as
  `misosoup/misosoup_202412.yaml` in the companion data repository
  <https://github.com/admarhi/ramen-data>. The YAML was provided by
  Alberto Pascual-Garcia.
- **Subset scope.** The 56 solutions span three focal strains
  (`A1R12`, `A3R04`, `B3M02`) on three carbon sources: acetate
  (`ac`, 26 solutions), citrate (`cit`, 26), and
  fructose-6-phosphate (`f6p`, 4). Each run uses a minimal medium
  with the named substrate as the sole carbon source.
- **Upstream method parameters.** The MiSoSoup tool version,
  genome-scale models, exact medium composition, MILP solver, and
  random seed are set upstream and are not reproduced here. See the
  MiSoSoup publication and the `ramen-data` repository for the
  authoritative description of the enumeration run.

Because the solutions are alternative optima, two files sharing a
substrate / focal-strain prefix are not biologically independent
communities -- see `?misosoup24` for the full caveat.

## Regenerating the package data

`data/misosoup24.rda` is rebuilt from the CSV files in this
directory by `data-raw/example-communities.R`:

```r
source("data-raw/example-communities.R")
```

The CSV files themselves are a fixed snapshot of the upstream
`misosoup_202412.yaml`; that YAML can be parsed into a
`ConsortiumMetabolismSet` with `importMisosoup()`.

## See also

- `?misosoup24` -- full dataset documentation.
- `?importMisosoup` -- importer for raw MiSoSoup YAML output.
