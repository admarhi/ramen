test_that("ConsortiumMetabolism constructor works with basic data", {
    # Create simple test data
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        met = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )

    # Test constructor
    cm <- ConsortiumMetabolism(test_data, name = "test_cm")

    # Basic assertions
    expect_s4_class(cm, "ConsortiumMetabolism")
    expect_equal(cm@Name, "test_cm")
    expect_true(length(cm@Metabolites) > 0)
    expect_s3_class(cm@Edges, "data.frame")
})

test_that("ConsortiumMetabolism handles weighted vs unweighted networks", {
    # Unweighted data (all fluxes = 1)
    unweighted_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        met = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )

    cm_unweighted <- ConsortiumMetabolism(unweighted_data, name = "unweighted")
    expect_false(cm_unweighted@Weighted)

    # Weighted data (varying fluxes)
    weighted_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        met = c("m1", "m2", "m1", "m3"),
        flux = c(-2, 3, -1.5, 2.5)
    )

    cm_weighted <- ConsortiumMetabolism(weighted_data, name = "weighted")
    expect_true(cm_weighted@Weighted)
})

test_that("getSpecies returns correct species", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2", "s3", "s3"),
        met = c("m1", "m2", "m1", "m3", "m2", "m4"),
        flux = c(-1, 1, -1, 1, -1, 1)
    )

    cm <- ConsortiumMetabolism(test_data, name = "test")
    species <- getSpecies(cm)

    expect_type(species, "character")
    expect_true(length(species) >= 2)
    expect_true(all(c("s1", "s2", "s3") %in% species))
})

test_that("getMet returns metabolites", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        met = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )

    cm <- ConsortiumMetabolism(test_data, name = "test")
    mets <- getMet(cm)

    expect_type(mets, "character")
    expect_true(length(mets) > 0)
})

test_that("getEdges returns edge data", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        met = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )

    cm <- ConsortiumMetabolism(test_data, name = "test")
    edges <- getEdges(cm)

    expect_s3_class(edges, "data.frame")
    expect_true(nrow(edges) > 0)
})

test_that("synCM generates synthetic communities", {
    # Test synthetic community generation
    syn <- synCM(name = "synthetic", n_species = 5, max_met = 10)

    expect_s4_class(syn, "ConsortiumMetabolism")
    expect_equal(syn@Name, "synthetic")
    expect_true(nrow(syn@Edges) > 0)
})
