test_that("rcs returns error when input is not a numeric vector of length > 1", {
  expect_error(rcs("a"))
  expect_error(rcs(1))
})

