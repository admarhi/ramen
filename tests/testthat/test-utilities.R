test_that("pivotCM transforms data correctly", {
    # Create test data in uptake/secretion format
    test_data <- tibble::tibble(
        species = c("s1", "s2"),
        uptake = c("m1", "m2"),
        secretion = c("m2", "m3"),
        flux = c(1, 1)
    )

    # Pivot the data
    pivoted <- pivotCM(
        test_data,
        species = "species",
        from = "uptake",
        to = "secretion",
        flux = "flux"
    )

    # Check structure
    expect_s3_class(pivoted, "data.frame")
    expect_true("species" %in% colnames(pivoted))
    expect_true("met" %in% colnames(pivoted))
    expect_true("flux" %in% colnames(pivoted))

    # Should have both consumption (negative) and production (positive) rows
    expect_true(any(pivoted$flux < 0))
    expect_true(any(pivoted$flux > 0))
})

test_that("pivotCM handles custom column names", {
    test_data <- tibble::tibble(
        org = c("s1"),
        consumed = c("m1"),
        produced = c("m2"),
        flow = c(2.5)
    )

    pivoted <- pivotCM(
        test_data,
        species = "org",
        from = "consumed",
        to = "produced",
        flux = "flow"
    )

    expect_s3_class(pivoted, "data.frame")
    expect_equal(nrow(pivoted), 2) # One row for consumption, one for production
})

test_that("pivotCM errors on non-data.frame input", {
    expect_error(
        pivotCM(list(a = 1), "a", "b", "c", "d"),
        "must be a data.frame"
    )
})

test_that("pivotCM errors on missing columns", {
    tb <- tibble::tibble(species = "s1", flux = 1)
    expect_error(
        pivotCM(tb, "species", "uptake", "secretion", "flux"),
        "uptake.*secretion|secretion.*uptake"
    )
})

test_that("synCM errors on invalid name", {
    expect_error(synCM(123, 3, 5), "single character string")
    expect_error(synCM(c("a", "b"), 3, 5), "single character string")
})

test_that("synCM errors on non-positive n_species", {
    expect_error(synCM("x", 0, 5), "positive integer")
    expect_error(synCM("x", -1, 5), "positive integer")
})

test_that("synCM errors on max_met < 2", {
    expect_error(synCM("x", 3, 1), "integer >= 2")
    expect_error(synCM("x", 3, 0), "integer >= 2")
})

test_that("synCM errors on non-positive scale_fac", {
    expect_error(synCM("x", 3, 5, scale_fac = 0), "positive integer")
})

test_that("name<- works for ConsortiumMetabolism", {
    test_data <- tibble::tibble(
        species = c("s1", "s1"),
        met = c("m1", "m2"),
        flux = c(-1, 1)
    )

    cm <- ConsortiumMetabolism(test_data, name = "original")
    name(cm) <- "new_name"

    expect_equal(name(cm), "new_name")
})

test_that("description<- works for ConsortiumMetabolism", {
    test_data <- tibble::tibble(
        species = c("s1", "s1"),
        met = c("m1", "m2"),
        flux = c(-1, 1)
    )

    cm <- ConsortiumMetabolism(test_data, name = "test")
    description(cm) <- "Test description"

    # ConsortiumMetabolism doesn't have Description slot, but method should work
    expect_s4_class(cm, "ConsortiumMetabolism")
})

test_that("consortia returns consortium data", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        met = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )

    cm <- ConsortiumMetabolism(test_data, name = "test")
    co_data <- consortia(cm)

    expect_s3_class(co_data, "data.frame")
    expect_true(nrow(co_data) > 0)
})
