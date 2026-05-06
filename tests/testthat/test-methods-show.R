## ---- show(CM) per-species pathway range ------------------------------------

test_that("show(cm) prints per-species pathway range when n_species >= 2", {
    cm <- synCM("a", n_species = 4, max_met = 8, seed = 1)
    out <- paste(
        capture.output(show(cm), type = "message"),
        collapse = " "
    )
    expect_match(out, "Pathways per species")
    expect_match(out, "min")
    expect_match(out, "mean")
    expect_match(out, "max")
})

test_that("show(cm) omits per-species range for single-species CM", {
    cm <- synCM("solo", n_species = 1, max_met = 4, seed = 1)
    out <- paste(
        capture.output(show(cm), type = "message"),
        collapse = " "
    )
    expect_false(grepl("Pathways per species", out))
})

## ---- show(CMS) per-class breakdown -----------------------------------------

test_that("show(cms) prints pathway and species class breakdowns", {
    cms <- ConsortiumMetabolismSet(
        list(
            synCM("a", n_species = 4, max_met = 8, seed = 1),
            synCM("b", n_species = 3, max_met = 6, seed = 2),
            synCM("c", n_species = 5, max_met = 10, seed = 3),
            synCM("d", n_species = 4, max_met = 8, seed = 4)
        ),
        name = "demo",
        verbose = FALSE
    )
    out <- paste(
        capture.output(show(cms), type = "message"),
        collapse = " "
    )
    expect_match(out, "pan-cons")
    expect_match(out, "niche")
    expect_match(out, "core")
    expect_match(out, "aux")
    expect_match(out, "generalists")
    expect_match(out, "specialists")
    expect_match(out, "quantile")
})

test_that("show(cms) omits per-class lines for single-consortium set", {
    cms <- ConsortiumMetabolismSet(
        list(synCM("only", n_species = 4, max_met = 8, seed = 1)),
        name = "solo",
        verbose = FALSE
    )
    out <- paste(
        capture.output(show(cms), type = "message"),
        collapse = " "
    )
    expect_false(grepl("pan-cons", out))
})

## ---- show(CMA) levels ------------------------------------------------------

test_that("show(cma) pairwise prints shared / unique pathway counts", {
    cm1 <- synCM("A", n_species = 3, max_met = 6, seed = 1)
    cm2 <- synCM("B", n_species = 3, max_met = 6, seed = 2)
    cma <- align(cm1, cm2)
    out <- paste(
        capture.output(show(cma), type = "message"),
        collapse = " "
    )
    expect_match(out, "shared")
    expect_match(out, "query-only")
    expect_match(out, "reference-only")
})

test_that("show(cma) multiple prints level histogram and core size", {
    cms <- ConsortiumMetabolismSet(
        list(
            synCM("a", n_species = 3, max_met = 6, seed = 1),
            synCM("b", n_species = 3, max_met = 6, seed = 2),
            synCM("c", n_species = 3, max_met = 6, seed = 3),
            synCM("d", n_species = 3, max_met = 6, seed = 4)
        ),
        name = "demo",
        verbose = FALSE
    )
    cma <- align(cms)
    out <- paste(
        capture.output(show(cma), type = "message"),
        collapse = " "
    )
    expect_match(out, "Levels:")
    expect_match(out, "n=1")
    expect_match(out, "n=4")
    expect_match(out, "Core")
})
