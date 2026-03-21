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
