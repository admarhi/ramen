# Tests for importSmetana()
#
# Fixtures live in tests/testthat/fixtures/smetana/ and contain three real
# SMETANA `--detailed` output files copied from the Machado 2021 cooccurrence
# dataset:
#   - rq_72.tsv_detailed.tsv  (33 data rows, smallest)
#   - bq_1.tsv_detailed.tsv   (70 data rows, has bidirectional species-compound
#                              pairs — Corynebacterium_simulans_PES1 / M_nh4_e
#                              acts as both donor and receiver)
#   - rq_0.tsv_detailed.tsv   (150 data rows, larger)

smetana_dir <- testthat::test_path("fixtures/smetana")
bq_1_path <- file.path(smetana_dir, "bq_1.tsv_detailed.tsv")
rq_72_path <- file.path(smetana_dir, "rq_72.tsv_detailed.tsv")

# Small inline data.frame for unit-style tests. Mimics raw SMETANA columns.
.makeMiniSmetana <- function() {
    tibble::tibble(
        community = "test_comm",
        medium = "minimal",
        receiver = c("sA", "sA", "sB", "sC", "sC"),
        donor = c("sB", "sC", "sA", "sB", "sA"),
        compound = c("M_glc_e", "M_glc_e", "M_glc_e", "M_ac_e", "M_lac_e"),
        scs = c(0.5, 0.3, 0.7, 0.4, 0.2),
        mus = c(0.4, 0.4, 0.6, 0.5, 0.1),
        mps = c(1, 1, 1, 1, 1),
        smetana = c(0.20, 0.12, 0.42, 0.20, 0.02)
    )
}


# ---------------------------------------------------------------------------
# Dispatch: input type handling
# ---------------------------------------------------------------------------

test_that("data.frame input returns a ConsortiumMetabolism", {
    df <- .makeMiniSmetana()
    cm <- importSmetana(df, name = "mini", verbose = FALSE)
    expect_s4_class(cm, "ConsortiumMetabolism")
    expect_equal(cm@Name, "mini")
})

test_that("single file input returns a ConsortiumMetabolism", {
    cm <- importSmetana(rq_72_path, verbose = FALSE)
    expect_s4_class(cm, "ConsortiumMetabolism")
})

test_that("directory input returns a ConsortiumMetabolismSet", {
    cms <- importSmetana(smetana_dir, verbose = FALSE)
    expect_s4_class(cms, "ConsortiumMetabolismSet")
    expect_length(cms@Consortia, 3)
})

test_that("data.frame without name errors with informative message", {
    df <- .makeMiniSmetana()
    expect_error(
        importSmetana(df, verbose = FALSE),
        regexp = "name.*required|name.*data.frame"
    )
})

test_that("invalid input type errors with informative message", {
    expect_error(
        importSmetana(42L, verbose = FALSE),
        regexp = "file path.*directory.*data.frame"
    )
})

test_that("empty directory errors with informative message", {
    empty_dir <- tempfile("empty_smetana_")
    dir.create(empty_dir)
    on.exit(unlink(empty_dir, recursive = TRUE))
    expect_error(
        importSmetana(empty_dir, verbose = FALSE),
        regexp = "No.*\\.tsv|No matching"
    )
})


# ---------------------------------------------------------------------------
# Name handling
# ---------------------------------------------------------------------------

test_that("name is derived from filename when input is a file", {
    cm <- importSmetana(rq_72_path, verbose = FALSE)
    expect_equal(cm@Name, "rq_72")
})

test_that("user-supplied name is preserved when input is a file", {
    # Regression test: a previous bug nulled out user-supplied names because
    # the if/else for name derivation lacked an else branch.
    cm <- importSmetana(
        rq_72_path,
        name = "my_custom_name",
        verbose = FALSE
    )
    expect_equal(cm@Name, "my_custom_name")
})

test_that("name parameter is passed through to CMS in directory mode", {
    cms <- importSmetana(
        smetana_dir,
        name = "test_set",
        verbose = FALSE
    )
    expect_equal(cms@Name, "test_set")
})

test_that("CMS in directory mode contains correctly-named CMs", {
    cms <- importSmetana(smetana_dir, verbose = FALSE)
    cm_names <- vapply(cms@Consortia, \(cm) cm@Name, character(1))
    expect_setequal(cm_names, c("bq_1", "rq_0", "rq_72"))
})


# ---------------------------------------------------------------------------
# Column validation
# ---------------------------------------------------------------------------

test_that("missing required columns errors with informative message", {
    df <- .makeMiniSmetana() |> dplyr::select(-"donor")
    expect_error(
        importSmetana(df, name = "x", verbose = FALSE),
        regexp = "donor"
    )
})

test_that("missing smetana column errors when use_scores = TRUE", {
    df <- .makeMiniSmetana() |> dplyr::select(-"smetana")
    expect_error(
        importSmetana(df, name = "x", use_scores = TRUE, verbose = FALSE),
        regexp = "smetana"
    )
})

test_that("NA in required columns errors", {
    df <- .makeMiniSmetana()
    df$donor[1] <- NA
    expect_error(
        importSmetana(df, name = "x", verbose = FALSE),
        regexp = "missing|NA"
    )
})


# ---------------------------------------------------------------------------
# Core transformation: binary mode
# ---------------------------------------------------------------------------

test_that("binary mode produces an unweighted CM", {
    df <- .makeMiniSmetana()
    cm <- importSmetana(
        df,
        name = "mini",
        use_scores = FALSE,
        verbose = FALSE
    )
    expect_false(cm@Weighted)
})

test_that("binary mode produces +1/-1 fluxes in InputData", {
    df <- .makeMiniSmetana()
    cm <- importSmetana(
        df,
        name = "mini",
        use_scores = FALSE,
        verbose = FALSE
    )
    expect_setequal(unique(cm@InputData$flux), c(-1, 1))
})

test_that("binary mode deduplicates redundant donor-compound pairs", {
    # In .makeMiniSmetana(), (sB -> sA, M_glc_e) and (sC -> sA, M_glc_e) both
    # have sA as receiver of glc — should produce ONE consumption edge for sA.
    df <- .makeMiniSmetana()
    cm <- importSmetana(
        df,
        name = "mini",
        use_scores = FALSE,
        verbose = FALSE
    )
    sa_glc <- cm@InputData |>
        dplyr::filter(
            .data$species == "sA",
            .data$met == "glc",
            .data$flux < 0
        )
    expect_equal(nrow(sa_glc), 1L)
})


# ---------------------------------------------------------------------------
# Bidirectional species-compound pairs (critical case)
# ---------------------------------------------------------------------------

test_that("species as donor AND receiver for same compound yields two edges", {
    # In bq_1, Corynebacterium_simulans_PES1 both donates AND receives M_nh4_e.
    # After collapse, there should be TWO edges for that species+metabolite:
    # one with flux = +1 (production) and one with flux = -1 (consumption).
    cm <- importSmetana(bq_1_path, verbose = FALSE)
    nh4_edges <- cm@InputData |>
        dplyr::filter(
            .data$species == "Corynebacterium_simulans_PES1",
            .data$met == "nh4"
        )
    expect_equal(nrow(nh4_edges), 2L)
    expect_setequal(sign(nh4_edges$flux), c(-1, 1))
})

test_that("inline bidirectional pair is preserved", {
    # Build a minimal data.frame where sX donates AND receives the same
    # compound, to verify the collapse handles it without depending on a
    # specific real-data file.
    df <- tibble::tibble(
        community = "c",
        medium = "m",
        receiver = c("sX", "sY"),
        donor = c("sY", "sX"),
        compound = c("M_co2_e", "M_co2_e"),
        scs = c(0.5, 0.5),
        mus = c(0.5, 0.5),
        mps = c(1, 1),
        smetana = c(0.25, 0.25)
    )
    cm <- importSmetana(df, name = "x", verbose = FALSE)
    sx_co2 <- cm@InputData |>
        dplyr::filter(.data$species == "sX", .data$met == "co2")
    expect_equal(nrow(sx_co2), 2L)
    expect_setequal(sign(sx_co2$flux), c(-1, 1))
})


# ---------------------------------------------------------------------------
# Core transformation: weighted mode
# ---------------------------------------------------------------------------

test_that("weighted mode produces a weighted CM", {
    df <- .makeMiniSmetana()
    cm <- importSmetana(df, name = "mini", use_scores = TRUE, verbose = FALSE)
    expect_true(cm@Weighted)
})

test_that("weighted mode produces nonzero fractional fluxes", {
    df <- .makeMiniSmetana()
    cm <- importSmetana(df, name = "mini", use_scores = TRUE, verbose = FALSE)
    flux_abs <- abs(cm@InputData$flux)
    expect_true(all(flux_abs > 0))
    expect_true(any(flux_abs < 1)) # at least some fractional values
})

test_that("weighted mode aggregates duplicate species-compound pairs by sum", {
    # sA receives M_glc_e from sB (smetana=0.20) and sC (smetana=0.12).
    # Expected consumption flux for (sA, glc): -(0.20 + 0.12) = -0.32
    df <- .makeMiniSmetana()
    cm <- importSmetana(df, name = "mini", use_scores = TRUE, verbose = FALSE)
    sa_glc <- cm@InputData |>
        dplyr::filter(
            .data$species == "sA",
            .data$met == "glc",
            .data$flux < 0
        )
    expect_equal(nrow(sa_glc), 1L)
    expect_equal(sa_glc$flux, -0.32)
})


# ---------------------------------------------------------------------------
# Metabolite ID normalization
# ---------------------------------------------------------------------------

test_that("normalize_ids = TRUE strips M_ prefix and _e suffix", {
    df <- .makeMiniSmetana()
    cm <- importSmetana(df, name = "mini", verbose = FALSE)
    # Original IDs: M_glc_e, M_ac_e, M_lac_e -> glc, ac, lac
    expect_true(all(c("glc", "ac", "lac") %in% cm@Metabolites))
    expect_false(any(grepl("^M_", cm@Metabolites)))
    expect_false(any(grepl("_e$", cm@Metabolites)))
})

test_that("normalize_ids = FALSE leaves IDs untouched", {
    df <- .makeMiniSmetana()
    cm <- importSmetana(
        df,
        name = "mini",
        normalize_ids = FALSE,
        verbose = FALSE
    )
    expect_true(all(c("M_glc_e", "M_ac_e", "M_lac_e") %in% cm@Metabolites))
})


# ---------------------------------------------------------------------------
# CM/CMS structural integrity
# ---------------------------------------------------------------------------

test_that("resulting CM passes validity check", {
    cm <- importSmetana(rq_72_path, verbose = FALSE)
    expect_true(methods::validObject(cm))
})

test_that("resulting CMS passes validity check", {
    cms <- importSmetana(smetana_dir, verbose = FALSE)
    expect_true(methods::validObject(cms))
})
