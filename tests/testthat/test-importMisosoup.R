# Tests for importMisosoup()
#
# Fixtures live in tests/testthat/fixtures/misosoup/ and contain real
# MiSoSoup YAML slices:
#   - single.yaml              (new format, 2obut/min with 2 solutions)
#   - multi/substrate_*.yaml   (new format, 3 files with 1 solution each)
#   - legacy.yaml              (OLD 202412 format, ac/A3R04 with 1 solution)

misosoup_dir <- testthat::test_path("fixtures/misosoup")
single_path <- file.path(misosoup_dir, "single.yaml")
legacy_path <- file.path(misosoup_dir, "legacy.yaml")
multi_dir <- file.path(misosoup_dir, "multi")


# ---------------------------------------------------------------------------
# Dispatch: input type handling
# ---------------------------------------------------------------------------

test_that("single file (new format) returns a CMS with one or more CMs", {
    cms <- importMisosoup(single_path, verbose = FALSE)
    expect_s4_class(cms, "ConsortiumMetabolismSet")
    expect_gte(length(cms@Consortia), 1L)
})

test_that("single file (legacy format) returns a CMS", {
    cms <- importMisosoup(legacy_path, verbose = FALSE)
    expect_s4_class(cms, "ConsortiumMetabolismSet")
    expect_gte(length(cms@Consortia), 1L)
})

test_that("directory input returns a CMS merging all files", {
    cms <- importMisosoup(multi_dir, verbose = FALSE)
    expect_s4_class(cms, "ConsortiumMetabolismSet")
    expect_equal(length(cms@Consortia), 3L)
})

test_that("pre-loaded list returns a CMS when name is supplied", {
    raw <- yaml::read_yaml(single_path)
    cms <- importMisosoup(raw, name = "test_list", verbose = FALSE)
    expect_s4_class(cms, "ConsortiumMetabolismSet")
    expect_equal(cms@Name, "test_list")
})

test_that("pre-loaded list without name errors", {
    raw <- yaml::read_yaml(single_path)
    expect_error(
        importMisosoup(raw, verbose = FALSE),
        regexp = "name.*required"
    )
})

test_that("invalid input type errors with informative message", {
    expect_error(
        importMisosoup(42L, verbose = FALSE),
        regexp = "file path.*directory.*list"
    )
})

test_that("empty directory errors with informative message", {
    empty_dir <- tempfile("empty_misosoup_")
    dir.create(empty_dir)
    on.exit(unlink(empty_dir, recursive = TRUE))
    expect_error(
        importMisosoup(empty_dir, verbose = FALSE),
        regexp = "No.*\\.yaml"
    )
})


# ---------------------------------------------------------------------------
# Name handling
# ---------------------------------------------------------------------------

test_that("CMS name is derived from filename for single-file input", {
    cms <- importMisosoup(single_path, verbose = FALSE)
    expect_equal(cms@Name, "single")
})

test_that("CMS name is derived from directory basename for dir input", {
    cms <- importMisosoup(multi_dir, verbose = FALSE)
    expect_equal(cms@Name, "multi")
})

test_that("user-supplied name overrides derived name", {
    cms <- importMisosoup(
        single_path,
        name = "my_experiment",
        verbose = FALSE
    )
    expect_equal(cms@Name, "my_experiment")
})


# ---------------------------------------------------------------------------
# cons_id format (both old and new)
# ---------------------------------------------------------------------------

test_that("new-format cons_id includes 'min' as second-level key", {
    cms <- importMisosoup(single_path, verbose = FALSE)
    cm_names <- vapply(cms@Consortia, \(cm) cm@Name, character(1))
    # single.yaml is 2obut/min with 2 solutions
    expect_true(all(grepl("^2obut_min_\\d+$", cm_names)))
})

test_that("legacy-format cons_id includes focal strain", {
    cms <- importMisosoup(legacy_path, verbose = FALSE)
    cm_names <- vapply(cms@Consortia, \(cm) cm@Name, character(1))
    # legacy.yaml is ac/A3R04 (1 solution)
    expect_equal(cm_names, "ac_A3R04_1")
})


# ---------------------------------------------------------------------------
# Growth extraction (auto-detect Growth_ and R_Biomass_)
# ---------------------------------------------------------------------------

test_that("growth is extracted from R_Biomass_ rows (new format)", {
    cms <- importMisosoup(single_path, verbose = FALSE)
    cm <- cms@Consortia[[1]]
    g <- growth(cm)
    expect_true(length(g) > 0)
    expect_true(all(g >= 0))
    # New format species include BS_L, BS_W (with underscores)
    expect_true(any(grepl("BS_", names(g))))
})

test_that("growth is extracted from Growth_ rows (legacy format)", {
    cms <- importMisosoup(legacy_path, verbose = FALSE)
    cm <- cms@Consortia[[1]]
    g <- growth(cm)
    expect_true(length(g) > 0)
    expect_true(all(g >= 0))
})


# ---------------------------------------------------------------------------
# Escape sequence decoding (end-to-end; unit test for .decodeBiggEscapes
# lives in test-helpers-import.R)
# ---------------------------------------------------------------------------

test_that("escaped metabolite IDs round-trip through the full pipeline", {
    cms <- importMisosoup(single_path, verbose = FALSE)
    cm <- cms@Consortia[[1]]
    # The single.yaml file contains at least one escaped metabolite
    # (acisnzd(e) encoded as __40__e__41__). After decoding + normalizing,
    # it should appear in metabolites as plain "acisnzd" (the (e) is
    # stripped by .normalizeBiggIds).
    expect_true("acisnzd" %in% cm@Metabolites)
    # No leftover escape sequences in any metabolite ID
    expect_false(any(grepl("__\\d+__", cm@Metabolites, perl = TRUE)))
})


# ---------------------------------------------------------------------------
# Species-list matching (robust to underscores in species IDs)
# ---------------------------------------------------------------------------

test_that("species with underscores (new format BS_L) parse correctly", {
    cms <- importMisosoup(single_path, verbose = FALSE)
    cm <- cms@Consortia[[1]]
    sp <- species(cm)
    # Expect BS_L, BS_W or similar underscore-containing strain IDs
    expect_true(any(grepl("_", sp)))
})

test_that("species without underscores (legacy A1R12 style) parse correctly", {
    cms <- importMisosoup(legacy_path, verbose = FALSE)
    cm <- cms@Consortia[[1]]
    sp <- species(cm)
    expect_true(length(sp) >= 1)
    expect_false(any(is.na(sp)))
})


# ---------------------------------------------------------------------------
# Media metadata
# ---------------------------------------------------------------------------

test_that("media-level fluxes are stashed in metadata", {
    cms <- importMisosoup(single_path, verbose = FALSE)
    cm <- cms@Consortia[[1]]
    media <- S4Vectors::metadata(cm)$media
    expect_s3_class(media, "data.frame")
    expect_named(media, c("metabolite", "flux"))
    expect_gt(nrow(media), 0)
})


# ---------------------------------------------------------------------------
# normalize_ids flag
# ---------------------------------------------------------------------------

test_that("normalize_ids = TRUE strips R_EX_ prefix and _e suffix", {
    cms <- importMisosoup(legacy_path, verbose = FALSE)
    cm <- cms@Consortia[[1]]
    expect_false(any(grepl("^R_EX_", cm@Metabolites)))
    expect_false(any(grepl("_e$", cm@Metabolites)))
})

test_that("normalize_ids = FALSE keeps R_EX_ prefix intact", {
    cms <- importMisosoup(
        legacy_path,
        normalize_ids = FALSE,
        verbose = FALSE
    )
    cm <- cms@Consortia[[1]]
    # With raw IDs, at least some metabolites should still carry the
    # R_EX_ prefix (the exchange reaction name form).
    expect_true(any(grepl("^R_EX_", cm@Metabolites)))
})


# ---------------------------------------------------------------------------
# Structural integrity
# ---------------------------------------------------------------------------

test_that("all returned CMs are weighted (flux magnitudes preserved)", {
    cms <- importMisosoup(single_path, verbose = FALSE)
    weighted <- vapply(cms@Consortia, \(cm) cm@Weighted, logical(1))
    expect_true(all(weighted))
})

test_that("returned CMS passes validObject", {
    cms <- importMisosoup(single_path, verbose = FALSE)
    expect_true(methods::validObject(cms))
})

test_that("every CM in the returned CMS passes validObject", {
    cms <- importMisosoup(multi_dir, verbose = FALSE)
    valid <- vapply(
        cms@Consortia,
        \(cm) isTRUE(methods::validObject(cm)),
        logical(1)
    )
    expect_true(all(valid))
})


# ---------------------------------------------------------------------------
# Zero-growth handling
# ---------------------------------------------------------------------------

test_that("zero-growth solutions are silently skipped", {
    # Construct an in-memory data structure with one viable and one
    # zero-growth solution. The zero-growth one (empty community) should
    # be skipped; the viable one should produce a CM.
    raw <- yaml::read_yaml(legacy_path)
    raw$ac$A3R04 <- c(
        raw$ac$A3R04,
        list(list(community = list(), solution = list()))
    )
    cms <- importMisosoup(raw, name = "mixed", verbose = FALSE)
    expect_equal(length(cms@Consortia), 1L)
})
