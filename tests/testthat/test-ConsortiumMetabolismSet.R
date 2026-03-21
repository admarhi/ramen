test_that("ConsortiumMetabolismSet constructor works", {
  # Create two test consortia
  data1 <- tibble::tibble(
    species = c("s1", "s1", "s2", "s2"),
    met = c("m1", "m2", "m1", "m3"),
    flux = c(-1, 1, -1, 1)
  )

  data2 <- tibble::tibble(
    species = c("s3", "s3", "s4", "s4"),
    met = c("m1", "m2", "m1", "m4"),
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
  expect_s3_class(cms@Edges, "data.frame")
})

test_that("ConsortiumMetabolismSet getSpecies works", {
  # Create test consortia with overlapping species
  data1 <- tibble::tibble(
    species = c("s1", "s1", "s2", "s2"),
    met = c("m1", "m2", "m1", "m3"),
    flux = c(-1, 1, -1, 1)
  )

  data2 <- tibble::tibble(
    species = c("s2", "s2", "s3", "s3"),
    met = c("m1", "m2", "m1", "m4"),
    flux = c(-1, 1, -1, 1)
  )

  cm1 <- ConsortiumMetabolism(data1, name = "cm1")
  cm2 <- ConsortiumMetabolism(data2, name = "cm2")
  cms <- ConsortiumMetabolismSet(list(cm1, cm2), name = "test")

  # Get all species (returns a tibble with species and n_edges columns)
  all_species <- getSpecies(cms, type = "all")
  expect_s3_class(all_species, "tbl_df")
  expect_true("s1" %in% all_species$species)
  expect_true("s2" %in% all_species$species)
  expect_true("s3" %in% all_species$species)
})

test_that("ConsortiumMetabolismSet getEdges works", {
  data1 <- tibble::tibble(
    species = c("s1", "s1", "s2", "s2"),
    met = c("m1", "m2", "m1", "m3"),
    flux = c(-1, 1, -1, 1)
  )

  data2 <- tibble::tibble(
    species = c("s1", "s1", "s3", "s3"),
    met = c("m1", "m2", "m2", "m4"),
    flux = c(-1, 1, -1, 1)
  )

  cm1 <- ConsortiumMetabolism(data1, name = "cm1")
  cm2 <- ConsortiumMetabolism(data2, name = "cm2")
  cms <- ConsortiumMetabolismSet(list(cm1, cm2), name = "test")

  # Get all edges
  edges <- getEdges(cms, type = "all")
  expect_s3_class(edges, "data.frame")
  expect_true(nrow(edges) > 0)
})

test_that("setName and setDesc work for ConsortiumMetabolismSet", {
  data1 <- tibble::tibble(
    species = c("s1", "s1"),
    met = c("m1", "m2"),
    flux = c(-1, 1)
  )

  cm1 <- ConsortiumMetabolism(data1, name = "cm1")
  cms <- ConsortiumMetabolismSet(list(cm1), name = "original")

  # Test setName
  cms_renamed <- setName(cms, "new_name")
  expect_equal(cms_renamed@Name, "new_name")

  # Test setDesc
  cms_described <- setDesc(cms, "test description")
  expect_equal(cms_described@Description, "test description")
})

test_that("cluster method works for ConsortiumMetabolismSet", {
  # Create test data
  data1 <- tibble::tibble(
    species = c("s1", "s1", "s2", "s2"),
    met = c("m1", "m2", "m1", "m3"),
    flux = c(-1, 1, -1, 1)
  )

  data2 <- tibble::tibble(
    species = c("s3", "s3", "s4", "s4"),
    met = c("m1", "m2", "m1", "m4"),
    flux = c(-1, 1, -1, 1)
  )

  cm1 <- ConsortiumMetabolism(data1, name = "cm1")
  cm2 <- ConsortiumMetabolism(data2, name = "cm2")
  cms <- ConsortiumMetabolismSet(list(cm1, cm2), name = "test")

  # Test clustering
  cms_clustered <- cluster(cms)

  expect_s4_class(cms_clustered, "ConsortiumMetabolismSet")
  expect_true(length(cms_clustered@Dendrogram) > 0)
})

## ---- CMS BinaryMatrices and OverlapMatrix -----------------------------------

test_that("CMS BinaryMatrices slot is populated", {
    cm1 <- synCM("a", n_species = 3, max_met = 5)
    cm2 <- synCM("b", n_species = 3, max_met = 5)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2), name = "test"
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
        list(cm1, cm2), name = "test"
    )
    ## OverlapMatrix should be a 2x2 matrix
    expect_equal(nrow(cms@OverlapMatrix), 2L)
    expect_equal(ncol(cms@OverlapMatrix), 2L)
    expect_true(is.matrix(cms@OverlapMatrix))
})

test_that("single-consortium CMS cannot be created (needs >= 2)", {
    cm1 <- synCM("a", n_species = 3, max_met = 5)
    ## CMS constructor requires >= 2 consortia for hclust
    expect_error(
        ConsortiumMetabolismSet(list(cm1), name = "single"),
        "n >= 2"
    )
})
