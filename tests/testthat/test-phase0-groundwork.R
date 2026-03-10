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

test_that("align() generic dispatches to CM,CM stub", {
    cm1 <- synCM("a", n_species = 3, max_met = 5)
    cm2 <- synCM("b", n_species = 3, max_met = 5)
    expect_error(
        align(cm1, cm2),
        "not yet implemented"
    )
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

test_that("CMS BinaryMatrices slot is populated", {
    cm1 <- synCM("a", n_species = 3, max_met = 5)
    cm2 <- synCM("b", n_species = 3, max_met = 5)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2), name = "test"
    )
    expect_true(length(cms@BinaryMatrices) == 2L)
    ## All matrices same dimension
    dims <- vapply(
        cms@BinaryMatrices,
        nrow,
        integer(1L)
    )
    expect_true(all(dims == dims[1L]))
})

test_that("CMS OverlapMatrix unchanged after refactor", {
    set.seed(42)
    cm1 <- synCM("a", n_species = 4, max_met = 8, seed = 42)
    cm2 <- synCM("b", n_species = 4, max_met = 8, seed = 43)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2), name = "test"
    )
    ## OverlapMatrix should be a 2x2 matrix
    expect_equal(nrow(cms@OverlapMatrix), 2L)
    expect_equal(ncol(cms@OverlapMatrix), 2L)
    expect_true(is.matrix(cms@OverlapMatrix))
})

test_that("deprecated functions raise errors", {
    expect_error(compareAlignments(), "deprecated")
    expect_error(plotAlignmentHeatmap(NULL), "deprecated")
    expect_error(plotAlignmentNetwork(NULL, 0.5), "deprecated")
})

test_that("single-consortium CMS cannot be created (needs >= 2)", {
    cm1 <- synCM("a", n_species = 3, max_met = 5)
    ## CMS constructor requires >= 2 consortia for hclust
    expect_error(
        ConsortiumMetabolismSet(list(cm1), name = "single"),
        "n >= 2"
    )
})
