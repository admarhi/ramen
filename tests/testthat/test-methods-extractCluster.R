## extractCluster input validation + happy path
## (this file did not exist previously; covr reported 0.00%
## coverage for R/methods-extractCluster.R)

.make_cms <- function(n = 3L) {
    cms_list <- lapply(seq_len(n), function(i) {
        synCM(
            paste0("c", i),
            n_species = 3L,
            max_met = 5L,
            seed = i
        )
    })
    ConsortiumMetabolismSet(
        cms_list,
        name = "set",
        verbose = FALSE
    )
}

test_that("extractCluster happy path returns a CMS with expected leaves", {
    cms <- .make_cms(3L)
    sub <- suppressMessages(extractCluster(cms, node_id = 1L))
    expect_s4_class(sub, "ConsortiumMetabolismSet")
    leaf_names <- vapply(sub@Consortia, name, character(1L))
    expect_true(length(leaf_names) >= 1L)
    expect_true(all(
        leaf_names %in%
            vapply(cms@Consortia, name, character(1L))
    ))
})

test_that("extractCluster errors on out-of-bounds node_id", {
    cms <- .make_cms(3L)
    expect_error(
        extractCluster(cms, node_id = 999L),
        "between 1 and"
    )
})

test_that("extractCluster errors on NA / NULL / negative / character node_id", {
    cms <- .make_cms(3L)
    expect_error(
        extractCluster(cms, node_id = NA),
        "single integer"
    )
    expect_error(
        extractCluster(cms, node_id = NA_integer_),
        "single integer"
    )
    expect_error(
        extractCluster(cms, node_id = NULL),
        "single integer"
    )
    expect_error(
        extractCluster(cms, node_id = -1L),
        "single integer"
    )
    expect_error(
        extractCluster(cms, node_id = "1"),
        "single integer"
    )
    expect_error(
        extractCluster(cms, node_id = c(1L, 2L)),
        "single integer"
    )
})

test_that("extractCluster default name uses 'Cluster N from <set>'", {
    cms <- .make_cms(3L)
    sub <- suppressMessages(extractCluster(cms, node_id = 1L))
    expect_equal(sub@Name, "Cluster 1 from set")
})

test_that("extractCluster errors on single-CMS (empty NodeData)", {
    cm <- synCM("solo", n_species = 3L, max_met = 5L, seed = 1L)
    cms <- ConsortiumMetabolismSet(
        list(cm),
        name = "lone",
        verbose = FALSE
    )
    expect_error(
        extractCluster(cms, node_id = 1L),
        "at least 2 consortia"
    )
})

test_that("extractCluster errors on non-character name / description", {
    cms <- .make_cms(3L)
    expect_error(
        extractCluster(cms, node_id = 1L, name = 1L),
        "name.*character"
    )
    expect_error(
        extractCluster(
            cms,
            node_id = 1L,
            description = list("x")
        ),
        "description.*character"
    )
})
