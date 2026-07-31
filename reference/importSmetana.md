# Import SMETANA Detailed Output

Imports raw SMETANA `--detailed` output (one TSV file per consortium)
into a ConsortiumMetabolism (single file/data.frame) or
ConsortiumMetabolismSet (directory of files). The raw format has columns
`community`, `medium`, `receiver`, `donor`, `compound`, `scs`, `mus`,
`mps`, `smetana`, where each row describes a cross-feeding interaction:
`donor` produces `compound`, `receiver` consumes it. This function
collapses those directed interactions into ramen's (species, metabolite,
flux) edge-list representation.

## Usage

``` r
importSmetana(
  data,
  name = NULL,
  use_scores = TRUE,
  normalize_ids = TRUE,
  verbose = TRUE
)
```

## Arguments

- data:

  Either a file path to a single SMETANA TSV, a directory path
  containing multiple SMETANA TSVs, or a pre-loaded data.frame with the
  required columns (`receiver`, `donor`, `compound`, and `smetana` if
  `use_scores = TRUE`).

- name:

  Consortium (or set) name as a length-1 character scalar. If `NULL`
  (default), names are derived from the input filename(s) by stripping
  the `.tsv_detailed.tsv` suffix. Required when `data` is a data.frame.
  If `data` is a directory, `name` is passed through to
  [`ConsortiumMetabolismSet()`](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismSet.md).

- use_scores:

  If `FALSE`, flux is binary (+1 for production, -1 for consumption). If
  `TRUE` (default), flux magnitude is the SMETANA score, aggregated by
  [`max()`](https://rdrr.io/r/base/Extremes.html) when the same
  (species, compound) pair appears across multiple interaction partners.

- normalize_ids:

  If `TRUE` (default), normalize compound IDs via
  [`.normalizeBiggIds()`](https://admarhi.github.io/ramen/reference/dot-normalizeBiggIds.md)
  (strips `M_` prefix and `_e` compartment suffix).

- verbose:

  If `TRUE` (default), emit progress messages via cli progress messages.

## Value

A `ConsortiumMetabolism` object if `data` is a single file or
data.frame; a `ConsortiumMetabolismSet` if `data` is a directory.

## See also

[`ConsortiumMetabolism()`](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md),
[`ConsortiumMetabolismSet()`](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismSet.md),
[`importMisosoup()`](https://admarhi.github.io/ramen/reference/importMisosoup.md)

## Examples

``` r
if (FALSE) {
# Single file -> ConsortiumMetabolism
cm <- importSmetana("path/to/bq_0.tsv_detailed.tsv")

# Directory of files -> ConsortiumMetabolismSet
cms <- importSmetana("path/to/bq_subsample/")

# Weighted fluxes
cm <- importSmetana("path/to/bq_0.tsv_detailed.tsv", use_scores = TRUE)
}
```
