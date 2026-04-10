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

test_that("synCM returns valid CM object", {
    cm <- synCM("syn", n_species = 4, max_met = 6, seed = 42L)
    expect_s4_class(cm, "ConsortiumMetabolism")
})

test_that("synCM seed reproducibility", {
    a <- synCM("s", 4, 8, seed = 99L, cm = FALSE)
    b <- synCM("s", 4, 8, seed = 99L, cm = FALSE)
    expect_identical(a, b)
})

test_that("synCM cross-feeding connectivity", {
    tbl <- synCM("c", n_species = 5, max_met = 10, seed = 7L, cm = FALSE)
    spp <- unique(tbl$species)
    # every species shares at least one metabolite with another
    for (sp in spp) {
        sp_mets <- tbl$metabolites[tbl$species == sp]
        others <- tbl$metabolites[tbl$species != sp]
        expect_true(
            length(intersect(sp_mets, others)) >= 1L,
            info = paste0("species ", sp, " is isolated")
        )
    }
})

test_that("synCM species have both consumed and produced", {
    tbl <- synCM(
        "d",
        n_species = 5,
        max_met = 10,
        seed = 11L,
        dead_ends = FALSE,
        cm = FALSE
    )
    spp <- unique(tbl$species)
    for (sp in spp) {
        fl <- tbl$fluxes[tbl$species == sp]
        expect_true(any(fl > 0), info = paste0(sp, " has no production"))
        expect_true(any(fl < 0), info = paste0(sp, " has no consumption"))
    }
})

test_that("synCM flux magnitudes are strictly positive", {
    tbl <- synCM("f", n_species = 4, max_met = 8, seed = 13L, cm = FALSE)
    expect_true(all(abs(tbl$fluxes) > 0))
})

test_that("synCM handles n_species = 2 edge case", {
    cm <- synCM("edge", n_species = 2, max_met = 4, seed = 21L)
    expect_s4_class(cm, "ConsortiumMetabolism")
})

test_that("synCM cm=FALSE returns tibble", {
    tbl <- synCM("t", n_species = 3, max_met = 5, seed = 31L, cm = FALSE)
    expect_s3_class(tbl, "tbl_df")
    expect_true(all(
        c("species", "metabolites", "fluxes") %in%
            colnames(tbl)
    ))
})

test_that("name<- works for ConsortiumMetabolism", {
    test_data <- tibble::tibble(
        species = c("s1", "s1"),
        metabolite = c("m1", "m2"),
        flux = c(-1, 1)
    )

    cm <- ConsortiumMetabolism(test_data, name = "original")
    name(cm) <- "new_name"

    expect_equal(name(cm), "new_name")
})

test_that("description<- works for ConsortiumMetabolism", {
    test_data <- tibble::tibble(
        species = c("s1", "s1"),
        metabolite = c("m1", "m2"),
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
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )

    cm <- ConsortiumMetabolism(test_data, name = "test")
    co_data <- consortia(cm)

    expect_s3_class(co_data, "data.frame")
    expect_true(nrow(co_data) > 0)
})
