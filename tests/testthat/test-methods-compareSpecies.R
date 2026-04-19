## ---- compareSpecies (same CM) -----------------------------------------------

test_that("compareSpecies same-CM returns named list", {
    cm <- synCM("test", n_species = 3, max_met = 5, seed = 1)
    sp <- species(cm)
    result <- compareSpecies(cm, sp[1], sp[2])
    expect_type(result, "list")
    expect_named(
        result,
        c("fos", "jaccard", "n_shared", "n_unique_sp1", "n_unique_sp2")
    )
})

test_that("compareSpecies same-CM scores in [0, 1]", {
    cm <- synCM("test", n_species = 4, max_met = 6, seed = 2)
    sp <- species(cm)
    result <- compareSpecies(cm, sp[1], sp[2])
    expect_true(result$fos >= 0 && result$fos <= 1)
    expect_true(result$jaccard >= 0 && result$jaccard <= 1)
})

test_that("compareSpecies same-CM identical species gives fos=1", {
    cm <- synCM("test", n_species = 3, max_met = 5, seed = 3)
    sp <- species(cm)
    result <- compareSpecies(cm, sp[1], sp[1])
    expect_equal(result$fos, 1)
    expect_equal(result$jaccard, 1)
    expect_equal(result$n_unique_sp1, 0L)
    expect_equal(result$n_unique_sp2, 0L)
})

test_that("compareSpecies same-CM n_shared + unique = total", {
    cm <- synCM("test", n_species = 3, max_met = 5, seed = 4)
    sp <- species(cm)
    result <- compareSpecies(cm, sp[1], sp[2])
    n1 <- result$n_shared + result$n_unique_sp1
    n2 <- result$n_shared + result$n_unique_sp2
    pw <- pathways(cm, verbose = FALSE)

    get_sp_paths <- function(cm_obj, s) {
        cm_obj@Pathways |>
            tidyr::unnest("data") |>
            dplyr::filter(.data$species == s) |>
            dplyr::select("consumed", "produced") |>
            dplyr::distinct()
    }
    expect_equal(n1, nrow(get_sp_paths(cm, sp[1])))
    expect_equal(n2, nrow(get_sp_paths(cm, sp[2])))
})

test_that("compareSpecies errors on unknown species", {
    cm <- synCM("test", n_species = 3, max_met = 5, seed = 1)
    sp <- species(cm)
    expect_error(
        compareSpecies(cm, sp[1], "zzz_not_there"),
        "not found"
    )
})

## ---- compareSpecies (cross-CM) ----------------------------------------------

test_that("compareSpecies cross-CM returns named list", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    sp1 <- species(cm1)
    sp2 <- species(cm2)
    result <- compareSpecies(cm1, cm2, sp1[1], sp2[1])
    expect_type(result, "list")
    expect_named(
        result,
        c("fos", "jaccard", "n_shared", "n_unique_sp1", "n_unique_sp2")
    )
})

test_that("compareSpecies cross-CM scores in [0, 1]", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    sp1 <- species(cm1)
    sp2 <- species(cm2)
    result <- compareSpecies(cm1, cm2, sp1[1], sp2[1])
    expect_true(result$fos >= 0 && result$fos <= 1)
    expect_true(result$jaccard >= 0 && result$jaccard <= 1)
})

test_that("compareSpecies cross-CM errors on unknown species", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    sp1 <- species(cm1)
    expect_error(
        compareSpecies(cm1, cm2, sp1[1], "zzz_not_there"),
        "not found"
    )
})
