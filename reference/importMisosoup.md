# Import MiSoSoup YAML Output

Imports raw MiSoSoup YAML output into a
[ConsortiumMetabolismSet](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismSet.md)
object. MiSoSoup YAMLs describe many consortia per file (one per
`substrate` / second-level-key / `solution` triple), so a single call
always produces a CMS regardless of whether the input is a single file,
a directory of files, or a pre-loaded nested list.

## Usage

``` r
importMisosoup(data, name = NULL, normalize_ids = TRUE, verbose = TRUE)
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
growth-row key (`Growth_<species>` in older runs, `R_Biomass_<species>`
in newer runs), COBRApy-style `__NN__` escape sequences in reaction IDs
(e.g. `__40__` for `(`), and the presence or absence of a focal-strain
layer below each substrate.

Media-level exchange fluxes (reactions without a species suffix) are
stashed in `metadata(cm)$media` for each consortium, preserving the
original bounds information if downstream code needs it.

## See also

[`importSmetana()`](https://admarhi.github.io/ramen/reference/importSmetana.md)
for the sibling SMETANA import.

## Examples

``` r
# \donttest{
# Single YAML file -> CMS
# cms <- importMisosoup("path/to/misosoup.yaml")

# Directory of YAML files -> CMS with consortia from all files
# cms <- importMisosoup("path/to/misosoup_dir/", name = "experiment1")

# Pre-loaded list -> CMS (name is required here)
# raw <- yaml::read_yaml("path/to/misosoup.yaml")
# cms <- importMisosoup(raw, name = "experiment1")
# }
```
