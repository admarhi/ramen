## ---- CMA class construction and validity ------------------------------------

test_that("CMA class can be instantiated with defaults", {
    cma <- ConsortiumMetabolismAlignment()
    expect_s4_class(cma, "ConsortiumMetabolismAlignment")
    expect_true(is.na(cma@Name))
    expect_true(is.na(cma@Type))
    expect_equal(cma@Metric, "FOS")
    expect_true(is.na(cma@PrimaryScore))
    expect_true(is.na(cma@Pvalue))
})

test_that("CMA validity rejects invalid PrimaryScore", {
    expect_error(
        ConsortiumMetabolismAlignment(PrimaryScore = 1.5),
        "PrimaryScore"
    )
    expect_error(
        ConsortiumMetabolismAlignment(PrimaryScore = -0.1),
        "PrimaryScore"
    )
})

test_that("CMA validity rejects invalid Type", {
    expect_error(
        ConsortiumMetabolismAlignment(Type = "invalid"),
        "Type"
    )
})

test_that("CMA validity rejects pairwise without names", {
    expect_error(
        ConsortiumMetabolismAlignment(
            Type = "pairwise",
            ReferenceName = "ref"
        ),
        "QueryName"
    )
})

test_that("CMA accepts valid pairwise configuration", {
    cma <- ConsortiumMetabolismAlignment(
        Type = "pairwise",
        QueryName = "query",
        ReferenceName = "ref",
        PrimaryScore = 0.75
    )
    expect_equal(cma@Type, "pairwise")
    expect_equal(cma@PrimaryScore, 0.75)
})

## ---- CMA accessors: pairwise ---------------------------------------------

test_that("scores() returns Scores list for pairwise", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    s <- scores(cma)
    expect_true(is.list(s))
    expect_true("FOS" %in% names(s))
    expect_true("jaccard" %in% names(s))
})

test_that("pathways(type='shared') returns data.frame for pairwise", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    sp <- pathways(cma, type = "shared")
    expect_true(is.data.frame(sp))
    expect_named(sp, c("consumed", "produced"))
})

test_that("pathways(type='unique') returns list for pairwise", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    up <- pathways(cma, type = "unique")
    expect_true(is.list(up))
    expect_true("query" %in% names(up))
    expect_true("reference" %in% names(up))
    expect_true(is.data.frame(up$query))
})

test_that("metabolites() returns character vector for CMA", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    m <- metabolites(cma)
    expect_true(is.character(m))
    expect_true(length(m) > 0L)
})

test_that("pathways() returns data.frame for CMA", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    e <- pathways(cma)
    expect_true(is.data.frame(e))
    expect_named(e, c("consumed", "produced"))
    ev <- pathways(cma, verbose = TRUE)
    expect_true(is.data.frame(ev))
    expect_true(ncol(ev) >= ncol(e))
})

test_that("consortia() errors for CMA", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    expect_error(consortia(cma), "not applicable")
})

## ---- CMA accessors: type guards ------------------------------------------

test_that("pathways(type='shared') errors for multiple CMA", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    cma <- align(cms)
    expect_error(
        pathways(cma, type = "shared"), "pairwise"
    )
})

test_that("pathways(type='unique') errors for multiple CMA", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    cma <- align(cms)
    expect_error(
        pathways(cma, type = "unique"), "pairwise"
    )
})

test_that("similarityMatrix() errors for pairwise CMA", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    expect_error(similarityMatrix(cma), "multiple")
})

test_that("prevalence() errors for pairwise CMA", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    expect_error(prevalence(cma), "multiple")
})

## ---- CMA accessors: multiple ---------------------------------------------

test_that("scores() returns Scores list for multiple", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    cma <- align(cms)
    s <- scores(cma)
    expect_true(is.list(s))
    expect_true("median" %in% names(s))
})

test_that("similarityMatrix() returns matrix for multiple", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    cma <- align(cms)
    sm <- similarityMatrix(cma)
    expect_true(is.matrix(sm))
    expect_equal(nrow(sm), 2L)
})

test_that("prevalence() returns data.frame for multiple", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    cma <- align(cms)
    p <- prevalence(cma)
    expect_true(is.data.frame(p))
    expect_true("nConsortia" %in% names(p))
})

test_that("pathways(type='consensus') returns data.frame for multiple", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    cma <- align(cms)
    cp <- pathways(cma, type = "consensus")
    expect_true(is.data.frame(cp))
    expect_true("nConsortia" %in% names(cp))
    expect_true("proportion" %in% names(cp))
})

test_that("pathways(type='consensus') errors for pairwise CMA", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    expect_error(
        pathways(cma, type = "consensus"), "multiple"
    )
})

test_that("pathways(type='shared', verbose=TRUE) returns full data", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    sp_concise <- pathways(cma, type = "shared")
    sp_verbose <- pathways(
        cma, type = "shared", verbose = TRUE
    )
    expect_true(is.data.frame(sp_verbose))
    expect_true(ncol(sp_verbose) >= ncol(sp_concise))
    expect_true("querySpecies" %in% names(sp_verbose))
})
