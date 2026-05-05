# Import MiSoSoup YAML Output

Imports raw MiSoSoup YAML output into a
[ConsortiumMetabolismSet](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismSet.md)
object. MiSoSoup YAMLs describe many consortia per file (one per
`substrate` / second-level-key / `solution` triple), so a single call
always produces a CMS regardless of whether the input is a single file,
a directory of files, or a pre-loaded nested list.

## Usage

``` r
importMisosoup(
  data,
  name = NULL,
  normalize_ids = TRUE,
  biomassPattern = "(?i)^.*?(?:^|_)(?:biomass|growth)_",
  verbose = TRUE
)
```

## Arguments

- data:

  Either a file path to a single MiSoSoup YAML, a directory path
  containing multiple YAML files, or a pre-loaded nested list from
  [`yaml::read_yaml()`](https://yaml.r-lib.org/reference/read_yaml.html).

- name:

  Name for the returned `ConsortiumMetabolismSet`. If `NULL` (default),
  the name is derived from the input filename (for a file) or directory
  basename (for a directory). **Required** when `data` is a pre-loaded
  list.

- normalize_ids:

  If `TRUE` (default), normalize metabolite IDs via
  [`.normalizeBiggIds()`](https://admarhi.github.io/ramen/reference/dot-normalizeBiggIds.md)
  (strip `R_EX_` prefix, `_e` suffix) and decode COBRApy `__NN__` escape
  sequences via
  [`.decodeBiggEscapes()`](https://admarhi.github.io/ramen/reference/dot-decodeBiggEscapes.md).
  If `FALSE`, keep metabolite IDs in their raw form.

- biomassPattern:

  Perl-compatible regex identifying biomass reaction rows. The pattern
  must match *and consume* the entire prefix from the start of the
  reaction ID up through the separator before the species ID, because it
  is used for both detection and species extraction (via
  [`sub()`](https://rdrr.io/r/base/grep.html)). Defaults to
  `"(?i)^.*?(?:^|_)(?:biomass|growth)_"`, which covers the three
  observed MiSoSoup variants (`Growth_<sp>`, `R_Biomass_<sp>`,
  `R_R_BIOMASS_<sp>`) case-insensitively and correctly rejects the
  `community_growth` summary row. Override this if your MiSoSoup output
  uses a non-standard biomass reaction name, e.g.
  `biomassPattern = "(?i)^.*?MYBIO_"`.

- verbose:

  If `TRUE` (default), emit progress messages via
  [`cli::cli_inform()`](https://cli.r-lib.org/reference/cli_abort.html)
  and show a progress bar when importing a directory.

## Value

A
[ConsortiumMetabolismSet](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismSet.md)
object containing one
[ConsortiumMetabolism](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md)
per viable consortium found in the input. Zero-growth solutions are
silently skipped.

## Details

The function auto-detects format differences between MiSoSoup runs: the
biomass-reaction name (via a configurable regex that defaults to
matching `Growth_<sp>`, `R_Biomass_<sp>`, and `R_R_BIOMASS_<sp>`
case-insensitively — MiSoSoup derives this name from the underlying
model, so it varies), COBRApy-style `__NN__` escape sequences in
reaction IDs (e.g. `__40__` for `(`), and the second-level key of each
substrate. That second-level key is either `min` — indicating a **CMSC**
(complete minimal supplying community: no viable focal strain because no
single member grows alone) — or a focal-strain ID, indicating an **MSC**
(minimal supplying community for that strain).

Per-consortium metadata set on each returned CM:

- `metadata(cm)$media` — media-level exchange fluxes (reactions without
  a species suffix), preserving the original bounds.

- `metadata(cm)$communityGrowth` — the scalar `community_growth` summary
  row from the MiSoSoup solution, or `NA_real_` if absent.

- `metadata(cm)$misosoupMode` — `"CMSC"` or `"MSC"`.

- `metadata(cm)$focalStrain` — focal-strain ID (character) for MSC
  consortia, `NA_character_` for CMSC.

## See also

[`importSmetana()`](https://admarhi.github.io/ramen/reference/importSmetana.md)
for the sibling SMETANA import.

## Examples

``` r
# Build a minimal in-memory MiSoSoup-shaped list and import it. In
# practice you would point at a real YAML file or a directory of
# YAMLs from a MiSoSoup run.
raw <- list(
    glucose = list(
        min = list(
            list(
                community = list(y_sp1 = 1, y_sp2 = 1),
                solution = list(
                    R_Biomass_sp1 = 0.4,
                    R_Biomass_sp2 = 0.3,
                    R_EX_glc__D_e_sp1_i = -10,
                    R_EX_ac_e_sp1_i = 5,
                    R_EX_ac_e_sp2_i = -4,
                    R_EX_co2_e_sp2_i = 3,
                    community_growth = 0.7
                )
            )
        )
    )
)
cms <- importMisosoup(raw, name = "demo", verbose = FALSE)
#> 
#> ── Creating CMS "demo" ─────────────────────────────────────────────────────────
#> ℹ Validating 1 <ConsortiumMetabolism> object
#> ✔ Validating 1 <ConsortiumMetabolism> object [11ms]
#> 
#> ℹ Collecting metabolites from 1 consortia
#> ✔ Collecting metabolites from 1 consortia [26ms]
#> 
#> ℹ Re-indexing 3 unique metabolites
#> ✔ Re-indexing 3 unique metabolites [26ms]
#> 
#> ℹ Expanding 1 binary matrices to 3-dimensional space
#> ✔ Expanding 1 binary matrices to 3-dimensional space [21ms]
#> 
#> ℹ Computing 3 x 3 levels matrix
#> ✔ Computing 3 x 3 levels matrix [23ms]
#> 
#> ℹ Computing pairwise overlap (0 pairs via crossprod)
#> ✔ Computing pairwise overlap (0 pairs via crossprod) [21ms]
#> 
#> ℹ Assembling pathway data from 1 consortia
#> ✔ Assembling pathway data from 1 consortia [37ms]
#> 
#> ℹ Collecting 1 consortium graphs
#> CMS "demo" created: 1 consortia, 3 metabolites (0.2s)
#> ✔ Collecting 1 consortium graphs [78ms]
#> 
cms
#> 
#> ── ConsortiumMetabolismSet 
#> Name: "demo"
#> 1 consortia, 2 species, 3 metabolites.
#> Community size (species): min 2, mean 2, max 2.
#> Community size (metabolites): min 3, mean 3, max 3.

# Real-world usage (commented; needs an actual MiSoSoup YAML on disk):
# cms <- importMisosoup("path/to/misosoup.yaml")
# cms <- importMisosoup("path/to/misosoup_dir/", name = "experiment1")
```
