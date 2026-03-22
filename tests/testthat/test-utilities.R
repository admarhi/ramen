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

test_that("setName works for ConsortiumMetabolism", {
  test_data <- tibble::tibble(
    species = c("s1", "s1"),
    met = c("m1", "m2"),
    flux = c(-1, 1)
  )

  cm <- ConsortiumMetabolism(test_data, name = "original")
  cm_renamed <- setName(cm, "new_name")

  expect_equal(cm_renamed@Name, "new_name")
})

test_that("setDesc works for ConsortiumMetabolism", {
  test_data <- tibble::tibble(
    species = c("s1", "s1"),
    met = c("m1", "m2"),
    flux = c(-1, 1)
  )

  cm <- ConsortiumMetabolism(test_data, name = "test")
  cm_described <- setDesc(cm, "Test description")

  # ConsortiumMetabolism doesn't have Description slot, but method should work
  expect_s4_class(cm_described, "ConsortiumMetabolism")
})

test_that("getCo returns consortium data", {
  test_data <- tibble::tibble(
    species = c("s1", "s1", "s2", "s2"),
    met = c("m1", "m2", "m1", "m3"),
    flux = c(-1, 1, -1, 1)
  )

  cm <- ConsortiumMetabolism(test_data, name = "test")
  co_data <- getCo(cm)

  expect_s3_class(co_data, "data.frame")
  expect_true(nrow(co_data) > 0)
})
