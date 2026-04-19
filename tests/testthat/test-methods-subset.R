## ---- Fixtures ---------------------------------------------------------------

.cms2 <- function() {
    cm1 <- synCM("a", n_species = 3, max_met = 6, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 6, seed = 2)
    ConsortiumMetabolismSet(list(cm1, cm2), name = "test", verbose = FALSE)
}

.cms3 <- function() {
    cm1 <- synCM("a", n_species = 3, max_met = 6, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 6, seed = 2)
    cm3 <- synCM("c", n_species = 3, max_met = 6, seed = 3)
    ConsortiumMetabolismSet(list(cm1, cm2, cm3), name = "test", verbose = FALSE)
}

## ---- CM: [ errors -----------------------------------------------------------

test_that("CM [ errors with informative message", {
    cm <- synCM("x", n_species = 3, max_met = 5, seed = 1)
    expect_error(cm[1:2], "not biologically meaningful")
    expect_error(cm[1:2, 1:2], "ConsortiumMetabolism")
})

## ---- CMS: basic subsetting --------------------------------------------------

test_that("CMS [ returns ConsortiumMetabolismSet", {
    cms <- .cms2()
    n <- nrow(cms@Metabolites)
    sub <- cms[seq_len(n - 1L)]
    expect_s4_class(sub, "ConsortiumMetabolismSet")
})

test_that("CMS [ Metabolites data.frame has only remaining mets", {
    cms <- .cms2()
    n <- nrow(cms@Metabolites)
    keep <- seq_len(n - 1L)
    sub <- cms[keep]
    expect_equal(nrow(sub@Metabolites), n - 1L)
    expect_true(all(sub@Metabolites$met %in% rownames(cms)[keep]))
})

test_that("CMS [ Metabolites met_ind is re-sequenced from 1", {
    cms <- .cms2()
    n <- nrow(cms@Metabolites)
    sub <- cms[seq_len(n - 1L)]
    expect_equal(sub@Metabolites$met_ind, seq_len(n - 1L))
})

test_that("CMS [ Pathways contains only remaining metabolites", {
    cms <- .cms2()
    n <- nrow(cms@Metabolites)
    sub <- cms[seq_len(n - 1L)]
    remaining <- sub@Metabolites$met
    if (nrow(sub@Pathways) > 0L) {
        expect_true(all(sub@Pathways$consumed %in% remaining))
        expect_true(all(sub@Pathways$produced %in% remaining))
    }
})

test_that("CMS [ BinaryMatrices dims match remaining metabolites", {
    cms <- .cms2()
    n <- nrow(cms@Metabolites)
    keep_n <- n - 1L
    sub <- cms[seq_len(keep_n)]
    dims <- vapply(sub@BinaryMatrices, nrow, integer(1L))
    expect_true(all(dims == keep_n))
})

test_that("CMS [ OverlapMatrix dims remain n_consortia x n_consortia", {
    cms <- .cms3()
    n <- nrow(cms@Metabolites)
    sub <- cms[seq_len(n - 1L)]
    expect_equal(dim(sub@OverlapMatrix), c(3L, 3L))
    expect_true(is.matrix(sub@OverlapMatrix))
})

test_that("CMS [ Consortia list length is unchanged", {
    cms <- .cms2()
    n <- nrow(cms@Metabolites)
    sub <- cms[seq_len(n - 1L)]
    expect_equal(length(sub@Consortia), 2L)
    expect_equal(
        vapply(sub@Consortia, \(x) x@Name, character(1L)),
        vapply(cms@Consortia, \(x) x@Name, character(1L))
    )
})

## ---- CMS: symmetric subsetting semantics ------------------------------------

test_that("CMS [i] and [i, i] produce identical results", {
    cms <- .cms2()
    keep <- seq_len(nrow(cms@Metabolites) - 1L)
    sub_i <- cms[keep]
    sub_ij <- cms[keep, keep]
    expect_equal(nrow(sub_i@Metabolites), nrow(sub_ij@Metabolites))
    expect_equal(sub_i@OverlapMatrix, sub_ij@OverlapMatrix)
})

test_that("CMS [i, j] with i != j errors", {
    cms <- .cms2()
    n <- nrow(cms@Metabolites)
    expect_error(cms[seq_len(n - 1L), seq_len(n - 2L)], "identical")
})

## ---- CMS: edge cases --------------------------------------------------------

test_that("CMS [ with empty index returns valid zero-metabolite object", {
    cms <- .cms2()
    sub <- cms[integer(0L)]
    expect_s4_class(sub, "ConsortiumMetabolismSet")
    expect_equal(nrow(sub@Metabolites), 0L)
    expect_equal(nrow(sub@Pathways), 0L)
    expect_equal(dim(sub@OverlapMatrix), c(2L, 2L))
})

test_that("CMS [ on single-consortium produces valid Dendrogram", {
    cm1 <- synCM("solo", n_species = 3, max_met = 6, seed = 1)
    cms <- ConsortiumMetabolismSet(list(cm1), name = "solo", verbose = FALSE)
    n <- nrow(cms@Metabolites)
    sub <- cms[seq_len(n - 1L)]
    expect_s4_class(sub, "ConsortiumMetabolismSet")
    expect_type(sub@Dendrogram, "list")
    expect_equal(nrow(sub@NodeData), 0L)
})

## ---- CMA: [ errors ----------------------------------------------------------

test_that("CMA [ errors with informative message", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    expect_error(cma[1L], "ConsortiumMetabolismAlignment")
    expect_error(cma[1L], "subsetting")
})
