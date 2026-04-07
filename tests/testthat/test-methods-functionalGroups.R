## ---- functionalGroups() return structure -----------------------------------

test_that("functionalGroups returns correct list structure", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 4, max_met = 6, seed = 2)
    cm3 <- synCM("c", n_species = 3, max_met = 5, seed = 3)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2, cm3),
        name = "test"
    )
    fg <- functionalGroups(cms)

    expect_true(is.list(fg))
    expected_names <- c(
        "dendrogram",
        "similarity_matrix",
        "incidence_matrix",
        "reactions_per_species"
    )
    expect_named(fg, expected_names, ignore.order = FALSE)
    ## No $plot element in refactored output
    expect_null(fg$plot)
})

## ---- Similarity matrix properties -----------------------------------------

test_that("similarity matrix has correct properties", {
    cm1 <- synCM("a", n_species = 4, max_met = 6, seed = 10)
    cm2 <- synCM("b", n_species = 3, max_met = 6, seed = 20)
    cm3 <- synCM("c", n_species = 3, max_met = 6, seed = 30)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2, cm3),
        name = "test"
    )
    fg <- functionalGroups(cms)
    sim <- fg$similarity_matrix

    ## Square matrix
    expect_equal(nrow(sim), ncol(sim))

    ## Diagonal is all 1s
    expect_equal(unname(diag(sim)), rep(1, nrow(sim)))

    ## Symmetric
    expect_equal(sim, t(sim))

    ## Values in [0, 1]
    expect_true(all(sim >= 0))
    expect_true(all(sim <= 1))

    ## Dimension matches total unique species count
    all_species <- unique(cms@Pathways$species)
    expect_equal(nrow(sim), length(all_species))
})

## ---- Dendrogram class -----------------------------------------------------

test_that("dendrogram has correct class", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cm3 <- synCM("c", n_species = 3, max_met = 5, seed = 3)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2, cm3),
        name = "test"
    )
    fg <- functionalGroups(cms)

    expect_true(is(fg$dendrogram, "dendrogram"))
})

## ---- Incidence matrix properties ------------------------------------------

test_that("incidence matrix is sparse with correct dims", {
    cm1 <- synCM("a", n_species = 4, max_met = 6, seed = 10)
    cm2 <- synCM("b", n_species = 3, max_met = 6, seed = 20)
    cm3 <- synCM("c", n_species = 4, max_met = 6, seed = 30)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2, cm3),
        name = "test"
    )
    fg <- functionalGroups(cms)
    inc <- fg$incidence_matrix

    ## Sparse matrix
    expect_true(is(inc, "sparseMatrix"))

    ## Rows = species, columns = pathways
    all_species <- sort(unique(cms@Pathways$species))
    expect_equal(nrow(inc), length(all_species))
    expect_equal(rownames(inc), all_species)

    ## All values are 0 or 1
    vals <- unique(as.vector(
        Matrix::summary(inc)$x
    ))
    expect_true(all(vals %in% c(0, 1)))

    ## Number of columns matches unique pathways
    unique_pathways <- unique(paste0(
        cms@Pathways$consumed,
        "-",
        cms@Pathways$produced
    ))
    expect_equal(ncol(inc), length(unique_pathways))
})

## ---- Linkage argument effect ----------------------------------------------

test_that("linkage argument changes dendrogram", {
    cm1 <- synCM("a", n_species = 4, max_met = 8, seed = 42)
    cm2 <- synCM("b", n_species = 4, max_met = 8, seed = 43)
    cm3 <- synCM("c", n_species = 4, max_met = 8, seed = 44)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2, cm3),
        name = "test"
    )

    fg_complete <- functionalGroups(
        cms, linkage = "complete"
    )
    fg_average <- functionalGroups(
        cms, linkage = "average"
    )

    ## Cophenetic distances should differ between methods
    coph_complete <- stats::cophenetic(
        fg_complete$dendrogram
    )
    coph_average <- stats::cophenetic(
        fg_average$dendrogram
    )

    ## Not identical (unless pathologically degenerate)
    expect_false(identical(coph_complete, coph_average))
})

## ---- Linkage validation ---------------------------------------------------

test_that("invalid linkage value errors", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )

    expect_error(
        functionalGroups(cms, linkage = "invalid_method")
    )
})

## ---- Deprecation warning for old viz args ---------------------------------

test_that("passing k triggers deprecation warning", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cm3 <- synCM("c", n_species = 3, max_met = 5, seed = 3)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2, cm3),
        name = "test"
    )

    expect_warning(
        functionalGroups(cms, k = 2),
        "plotFunctionalGroups"
    )
})

## ---- plotFunctionalGroups returns ggplot -----------------------------------

test_that("plotFunctionalGroups returns ggplot", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cm3 <- synCM("c", n_species = 3, max_met = 5, seed = 3)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2, cm3),
        name = "test"
    )

    fg <- functionalGroups(cms)
    p <- plotFunctionalGroups(fg, k = 2)
    expect_s3_class(p, "ggplot")
})

## ---- plotFunctionalGroups with label_colours ------------------------------

test_that("plotFunctionalGroups works with label_colours", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cm3 <- synCM("c", n_species = 3, max_met = 5, seed = 3)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2, cm3),
        name = "test"
    )

    fg <- functionalGroups(cms)

    ## Build label_colours tibble matching all species
    all_species <- rownames(fg$similarity_matrix)
    colours_tb <- tibble::tibble(
        label = all_species,
        colour = rep_len(
            c("red", "blue", "green"),
            length(all_species)
        )
    )

    expect_no_error(
        plotFunctionalGroups(
            fg,
            k = 2,
            label_colours = colours_tb
        )
    )

    p <- plotFunctionalGroups(
        fg,
        k = 2,
        label_colours = colours_tb
    )
    expect_s3_class(p, "ggplot")
})

## ---- plotFunctionalGroups errors on bad input -----------------------------

test_that("plotFunctionalGroups errors on non-list input", {
    expect_error(
        plotFunctionalGroups("not_a_list", k = 2),
        "fg"
    )
})

test_that("plotFunctionalGroups errors on list without dendrogram", {
    bad_input <- list(
        similarity_matrix = matrix(1, 2, 2),
        incidence_matrix = Matrix::sparseMatrix(
            i = 1, j = 1, x = 1, dims = c(2, 2)
        )
    )
    expect_error(
        plotFunctionalGroups(bad_input, k = 2),
        "dendrogram"
    )
})
