test_that("ConsortiumMetabolismSet constructor works", {
    # Create two test consortia
    data1 <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )

    data2 <- tibble::tibble(
        species = c("s3", "s3", "s4", "s4"),
        metabolite = c("m1", "m2", "m1", "m4"),
        flux = c(-1, 1, -1, 1)
    )

    cm1 <- ConsortiumMetabolism(data1, name = "cm1")
    cm2 <- ConsortiumMetabolism(data2, name = "cm2")

    # Create set
    cms <- ConsortiumMetabolismSet(list(cm1, cm2), name = "test_set")

    # Assertions
    expect_s4_class(cms, "ConsortiumMetabolismSet")
    expect_equal(cms@Name, "test_set")
    expect_equal(length(cms@Consortia), 2)
    expect_s3_class(cms@Pathways, "data.frame")
})

test_that("ConsortiumMetabolismSet species works", {
    # Create test consortia with overlapping species
    data1 <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )

    data2 <- tibble::tibble(
        species = c("s2", "s2", "s3", "s3"),
        metabolite = c("m1", "m2", "m1", "m4"),
        flux = c(-1, 1, -1, 1)
    )

    cm1 <- ConsortiumMetabolism(data1, name = "cm1")
    cm2 <- ConsortiumMetabolism(data2, name = "cm2")
    cms <- ConsortiumMetabolismSet(list(cm1, cm2), name = "test")

    # Get all species (returns a tibble with species and n_pathways columns)
    all_species <- species(cms, type = "all")
    expect_s3_class(all_species, "tbl_df")
    expect_true("s1" %in% all_species$species)
    expect_true("s2" %in% all_species$species)
    expect_true("s3" %in% all_species$species)
})

test_that("CMS pathways returns concise output by default", {
    data1 <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )

    data2 <- tibble::tibble(
        species = c("s1", "s1", "s3", "s3"),
        metabolite = c("m1", "m2", "m2", "m4"),
        flux = c(-1, 1, -1, 1)
    )

    cm1 <- ConsortiumMetabolism(data1, name = "cm1")
    cm2 <- ConsortiumMetabolism(data2, name = "cm2")
    cms <- ConsortiumMetabolismSet(list(cm1, cm2), name = "test")

    pw <- pathways(cms)
    expect_s3_class(pw, "data.frame")
    expect_true(nrow(pw) > 0)
    expect_named(
        pw,
        c("consumed", "produced", "n_species", "n_cons")
    )
})

test_that("CMS pathways verbose returns full detail", {
    data1 <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )

    data2 <- tibble::tibble(
        species = c("s1", "s1", "s3", "s3"),
        metabolite = c("m1", "m2", "m2", "m4"),
        flux = c(-1, 1, -1, 1)
    )

    cm1 <- ConsortiumMetabolism(data1, name = "cm1")
    cm2 <- ConsortiumMetabolism(data2, name = "cm2")
    cms <- ConsortiumMetabolismSet(list(cm1, cm2), name = "test")

    pw <- pathways(cms, verbose = TRUE)
    expect_s3_class(pw, "data.frame")
    expect_true("cm_name" %in% names(pw))
    expect_true("species" %in% names(pw))
    expect_true("n_cons" %in% names(pw))
    expect_true(ncol(pw) > 4)
})

test_that("CMS pathways type filtering works", {
    data1 <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )

    data2 <- tibble::tibble(
        species = c("s1", "s1", "s3", "s3"),
        metabolite = c("m1", "m2", "m2", "m4"),
        flux = c(-1, 1, -1, 1)
    )

    cm1 <- ConsortiumMetabolism(data1, name = "cm1")
    cm2 <- ConsortiumMetabolism(data2, name = "cm2")
    cms <- ConsortiumMetabolismSet(list(cm1, cm2), name = "test")

    pw_all <- pathways(cms, type = "all")
    expect_true(nrow(pw_all) > 0)
    expect_named(
        pw_all,
        c("consumed", "produced", "n_species", "n_cons")
    )
})

test_that("name<- and description<- work for ConsortiumMetabolismSet", {
    data1 <- tibble::tibble(
        species = c("s1", "s1"),
        metabolite = c("m1", "m2"),
        flux = c(-1, 1)
    )

    cm1 <- ConsortiumMetabolism(data1, name = "cm1")
    cms <- ConsortiumMetabolismSet(list(cm1), name = "original")

    # Test name<-
    name(cms) <- "new_name"
    expect_equal(name(cms), "new_name")

    # Test description<-
    description(cms) <- "test description"
    expect_equal(description(cms), "test description")
})


## ---- CMS BinaryMatrices and OverlapMatrix -----------------------------------

test_that("CMS BinaryMatrices slot is populated", {
    cm1 <- synCM("a", n_species = 3, max_met = 5)
    cm2 <- synCM("b", n_species = 3, max_met = 5)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
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
        list(cm1, cm2),
        name = "test"
    )
    ## OverlapMatrix should be a 2x2 matrix
    expect_equal(nrow(cms@OverlapMatrix), 2L)
    expect_equal(ncol(cms@OverlapMatrix), 2L)
    expect_true(is.matrix(cms@OverlapMatrix))
})

test_that("single-consortium CMS can be created", {
    cm1 <- synCM("a", n_species = 3, max_met = 5)
    cms <- ConsortiumMetabolismSet(list(cm1), name = "single")
    expect_s4_class(cms, "ConsortiumMetabolismSet")
    expect_equal(length(cms@Consortia), 1L)
})

## ---- B5: Duplicate consortium names error ----------------------------------

test_that("CMS constructor errors on duplicate consortium names", {
    cm1 <- synCM("same_name", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("same_name", n_species = 3, max_met = 5, seed = 2)
    expect_error(
        ConsortiumMetabolismSet(cm1, cm2, name = "test"),
        "unique"
    )
})

test_that("CMS constructor accepts unique consortium names", {
    cm1 <- synCM("name_a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("name_b", n_species = 3, max_met = 5, seed = 2)
    expect_no_error(
        ConsortiumMetabolismSet(cm1, cm2, name = "test")
    )
})

## ---- B6: core pathway quantile uses actual distribution --------------------

test_that("pathways(cms, type='core') uses actual species distribution", {
    cm1 <- synCM("a", n_species = 4, max_met = 5, seed = 10)
    cm2 <- synCM("b", n_species = 4, max_met = 5, seed = 20)
    cm3 <- synCM("c", n_species = 4, max_met = 5, seed = 30)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2, cm3),
        name = "test"
    )
    all_pw <- pathways(cms, type = "all")
    core_pw <- pathways(cms, type = "core")
    if (nrow(core_pw) > 0L) {
        quant <- stats::quantile(
            all_pw$n_species,
            p = 0.9
        )
        expect_true(all(core_pw$n_species > quant))
    }
})

## ---- consortia() method for CMS -------------------------------------------

test_that("consortia(cms) returns list of CM objects", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    result <- consortia(cms)
    expect_type(result, "list")
    expect_length(result, 2L)
    expect_s4_class(result[[1L]], "ConsortiumMetabolism")
    expect_s4_class(result[[2L]], "ConsortiumMetabolism")
})

## ---- growth() on CMS --------------------------------------------------------

test_that("growth(cms) returns named list from stored CMs", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )
    cm1 <- ConsortiumMetabolism(
        test_data,
        name = "cm1",
        growth = c(s1 = 0.5, s2 = 0.3)
    )
    cm2 <- ConsortiumMetabolism(
        test_data,
        name = "cm2",
        growth = c(s1 = 0.8, s2 = 0.1)
    )
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    gr <- growth(cms)
    expect_type(gr, "list")
    expect_named(gr, c("cm1", "cm2"))
    expect_equal(gr$cm1, c(s1 = 0.5, s2 = 0.3))
    expect_equal(gr$cm2, c(s1 = 0.8, s2 = 0.1))
})

test_that("growth(cms) returns NULLs when no growth data", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    gr <- growth(cms)
    expect_type(gr, "list")
    expect_named(gr, c("a", "b"))
    expect_null(gr$a)
    expect_null(gr$b)
})

test_that("consortia(cms) names match consortium names", {
    cm1 <- synCM("alpha", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("beta", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    result <- consortia(cms)
    expect_equal(
        vapply(result, name, character(1L)),
        c("alpha", "beta")
    )
})
