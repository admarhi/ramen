## ---- align() dispatch tests -------------------------------------------------

test_that("align() dispatches CM,CM and returns CMA", {
    cm1 <- synCM("a", n_species = 3, max_met = 5)
    cm2 <- synCM("b", n_species = 3, max_met = 5)
    cma <- align(cm1, cm2)
    expect_s4_class(cma, "ConsortiumMetabolismAlignment")
    expect_equal(cma@Type, "pairwise")
})

test_that("align() generic dispatches to CMS,missing stub", {
    cm1 <- synCM("a", n_species = 3, max_met = 5)
    cm2 <- synCM("b", n_species = 3, max_met = 5)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2), name = "test"
    )
    expect_error(
        align(cms),
        "not yet implemented"
    )
})

## ---- Deprecated function tests ----------------------------------------------

test_that("deprecated functions raise errors", {
    expect_error(compareAlignments(), "deprecated")
    expect_error(plotAlignmentHeatmap(NULL), "deprecated")
    expect_error(plotAlignmentNetwork(NULL, 0.5), "deprecated")
})

## ---- Integration tests: align(CM, CM) --------------------------------------

test_that("align identical CMs returns score 1", {
    cm <- synCM("test", n_species = 4, max_met = 6, seed = 42)
    cma <- align(cm, cm)
    expect_s4_class(cma, "ConsortiumMetabolismAlignment")
    expect_equal(cma@PrimaryScore, 1)
    expect_equal(cma@Type, "pairwise")
    expect_equal(nrow(cma@UniqueQuery), 0)
    expect_equal(nrow(cma@UniqueReference), 0)
})

test_that("align is symmetric for FOS", {
    cm1 <- synCM("a", n_species = 4, max_met = 6, seed = 10)
    cm2 <- synCM("b", n_species = 4, max_met = 6, seed = 20)
    cma_ab <- align(cm1, cm2, method = "FOS")
    cma_ba <- align(cm2, cm1, method = "FOS")
    expect_equal(cma_ab@PrimaryScore, cma_ba@PrimaryScore)
})

test_that("align is symmetric for jaccard", {
    cm1 <- synCM("a", n_species = 4, max_met = 6, seed = 10)
    cm2 <- synCM("b", n_species = 4, max_met = 6, seed = 20)
    cma_ab <- align(cm1, cm2, method = "jaccard")
    cma_ba <- align(cm2, cm1, method = "jaccard")
    expect_equal(cma_ab@PrimaryScore, cma_ba@PrimaryScore)
})

test_that("all scores in [0, 1]", {
    cm1 <- synCM("a", n_species = 5, max_met = 8, seed = 42)
    cm2 <- synCM("b", n_species = 5, max_met = 8, seed = 43)
    cma <- align(cm1, cm2)
    for (nm in names(cma@Scores)) {
        if (!is.na(cma@Scores[[nm]])) {
            expect_true(
                cma@Scores[[nm]] >= 0,
                info = paste(nm, "< 0")
            )
            expect_true(
                cma@Scores[[nm]] <= 1,
                info = paste(nm, "> 1")
            )
        }
    }
})

test_that("MAAS method works", {
    cm1 <- synCM("a", n_species = 5, max_met = 8, seed = 42)
    cm2 <- synCM("b", n_species = 5, max_met = 8, seed = 43)
    cma <- align(cm1, cm2, method = "MAAS")
    expect_equal(cma@Metric, "MAAS")
    expect_true(!is.na(cma@PrimaryScore))
})

test_that("align with computePvalue returns valid p-value", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 42)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 43)
    cma <- align(cm1, cm2, computePvalue = TRUE,
                 nPermutations = 99L)
    expect_true(!is.na(cma@Pvalue))
    expect_true(cma@Pvalue >= 0 && cma@Pvalue <= 1)
})

test_that("SharedPathways + UniqueQuery + UniqueRef = total edges", {
    cm1 <- synCM("a", n_species = 4, max_met = 6, seed = 10)
    cm2 <- synCM("b", n_species = 4, max_met = 6, seed = 20)
    cma <- align(cm1, cm2)
    total_shared <- nrow(cma@SharedPathways)
    total_uq <- nrow(cma@UniqueQuery)
    total_ur <- nrow(cma@UniqueReference)
    expect_equal(
        total_shared + total_uq,
        sum(assays(cm1)$Binary)
    )
    expect_equal(
        total_shared + total_ur,
        sum(assays(cm2)$Binary)
    )
})

test_that("FOS from align matches .binMatOverlap via CMS", {
    cm1 <- synCM("a", n_species = 4, max_met = 8, seed = 42)
    cm2 <- synCM("b", n_species = 4, max_met = 8, seed = 43)
    cma <- align(cm1, cm2, method = "FOS")
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2), name = "test"
    )
    cms_fos <- 1 - cms@OverlapMatrix[1, 2]
    expect_equal(cma@PrimaryScore, cms_fos)
})

## ---- Edge case tests -------------------------------------------------------

test_that("align errors on invalid method", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    expect_error(align(cm1, cm2, method = "invalid"))
})

test_that("align handles single-edge CMs", {
    cm1 <- ConsortiumMetabolism(
        data = tibble::tibble(
            species = c("sp1", "sp1"),
            met = c("A", "B"),
            flux = c(-1, 1)
        ),
        name = "single1"
    )
    cm2 <- ConsortiumMetabolism(
        data = tibble::tibble(
            species = c("sp2", "sp2"),
            met = c("A", "B"),
            flux = c(-1, 1)
        ),
        name = "single2"
    )
    cma <- align(cm1, cm2)
    expect_equal(cma@PrimaryScore, 1)
})
