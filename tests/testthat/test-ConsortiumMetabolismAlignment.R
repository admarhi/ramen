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

test_that("sharedPathways() returns data.frame for pairwise", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    sp <- sharedPathways(cma)
    expect_true(is.data.frame(sp))
    expect_true("consumed" %in% names(sp))
    expect_true("produced" %in% names(sp))
})

test_that("uniquePathways() returns list for pairwise", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    up <- uniquePathways(cma)
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
})

test_that("consortia() errors for CMA", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cma <- align(cm1, cm2)
    expect_error(consortia(cma), "not applicable")
})

## ---- CMA accessors: type guards ------------------------------------------

test_that("sharedPathways() errors for multiple CMA", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    cma <- align(cms)
    expect_error(sharedPathways(cma), "pairwise")
})

test_that("uniquePathways() errors for multiple CMA", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    cma <- align(cms)
    expect_error(uniquePathways(cma), "pairwise")
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

test_that("consensusPathways() returns data.frame for multiple", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    cma <- align(cms)
    cp <- consensusPathways(cma)
    expect_true(is.data.frame(cp))
    expect_identical(cp, prevalence(cma))
})
