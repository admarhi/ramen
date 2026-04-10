test_that("ConsortiumMetabolism constructor works with basic data", {
    # Create simple test data
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )

    # Test constructor
    cm <- ConsortiumMetabolism(test_data, name = "test_cm")

    # Basic assertions
    expect_s4_class(cm, "ConsortiumMetabolism")
    expect_equal(cm@Name, "test_cm")
    expect_true(length(cm@Metabolites) > 0)
    expect_s3_class(cm@Pathways, "data.frame")
})

test_that("ConsortiumMetabolism handles weighted vs unweighted networks", {
    # Unweighted data (all fluxes = 1)
    unweighted_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )

    cm_unweighted <- ConsortiumMetabolism(unweighted_data, name = "unweighted")
    expect_false(cm_unweighted@Weighted)

    # Weighted data (varying fluxes)
    weighted_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-2, 3, -1.5, 2.5)
    )

    cm_weighted <- ConsortiumMetabolism(weighted_data, name = "weighted")
    expect_true(cm_weighted@Weighted)
})

test_that("species returns correct species", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2", "s3", "s3"),
        metabolite = c("m1", "m2", "m1", "m3", "m2", "m4"),
        flux = c(-1, 1, -1, 1, -1, 1)
    )

    cm <- ConsortiumMetabolism(test_data, name = "test")
    sp <- species(cm)

    expect_type(sp, "character")
    expect_true(length(sp) >= 2)
    expect_true(all(c("s1", "s2", "s3") %in% sp))
})

test_that("metabolites returns metabolites", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )

    cm <- ConsortiumMetabolism(test_data, name = "test")
    mets <- metabolites(cm)

    expect_type(mets, "character")
    expect_true(length(mets) > 0)
})

test_that("pathways returns concise output by default", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )

    cm <- ConsortiumMetabolism(test_data, name = "test")
    pw <- pathways(cm)

    expect_s3_class(pw, "data.frame")
    expect_true(nrow(pw) > 0)
    expect_named(pw, c("consumed", "produced", "n_species"))
})

test_that("pathways verbose returns full detail", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )

    cm <- ConsortiumMetabolism(test_data, name = "test")
    pw <- pathways(cm, verbose = TRUE)

    expect_s3_class(pw, "data.frame")
    expect_true(ncol(pw) > 3)
    expect_true("c_sum" %in% names(pw))
    expect_true("c_eff" %in% names(pw))
    expect_true("data" %in% names(pw))
})

test_that("CM constructor errors on non-data.frame input", {
    expect_error(
        ConsortiumMetabolism("not a df", name = "x"),
        "must be a data.frame"
    )
    expect_error(
        ConsortiumMetabolism(list(a = 1), name = "x"),
        "must be a data.frame"
    )
})

## ---- growth parameter -------------------------------------------------------

test_that("CM constructor accepts growth parameter", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )
    gr <- c(s1 = 0.5, s2 = 0.3)
    cm <- ConsortiumMetabolism(
        test_data,
        name = "test",
        growth = gr
    )
    expect_equal(growth(cm), gr)
})

test_that("growth() returns NULL when no growth supplied", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )
    cm <- ConsortiumMetabolism(test_data, name = "test")
    expect_null(growth(cm))
})

test_that("growth validation: must be numeric", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )
    expect_error(
        ConsortiumMetabolism(
            test_data,
            name = "test",
            growth = c(s1 = "high", s2 = "low")
        ),
        "must be a numeric"
    )
})

test_that("growth validation: must be named", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )
    expect_error(
        ConsortiumMetabolism(
            test_data,
            name = "test",
            growth = c(0.5, 0.3)
        ),
        "named numeric"
    )
})

test_that("growth validation: names must match species", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )
    expect_error(
        ConsortiumMetabolism(
            test_data,
            name = "test",
            growth = c(s1 = 0.5, s99 = 0.3)
        ),
        "not found"
    )
})

test_that("growth validation: no duplicate names", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )
    gr <- c(s1 = 0.5, s1 = 0.3)
    expect_error(
        ConsortiumMetabolism(
            test_data,
            name = "test",
            growth = gr
        ),
        "duplicate"
    )
})

test_that("growth validation: non-negative values", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )
    expect_error(
        ConsortiumMetabolism(
            test_data,
            name = "test",
            growth = c(s1 = 0.5, s2 = -0.1)
        ),
        "non-negative"
    )
})

test_that("growth accepts subset of species", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )
    gr <- c(s1 = 0.5)
    cm <- ConsortiumMetabolism(
        test_data,
        name = "test",
        growth = gr
    )
    expect_equal(growth(cm), gr)
})

test_that("synCM generates synthetic communities", {
    # Test synthetic community generation
    syn <- synCM(name = "synthetic", n_species = 5, max_met = 10)

    expect_s4_class(syn, "ConsortiumMetabolism")
    expect_equal(syn@Name, "synthetic")
    expect_true(nrow(syn@Pathways) > 0)
})
