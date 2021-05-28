test_that("rcs returns error when input is not a numeric vector of length > 1", {
  expect_error(rcs("a"))
  expect_error(rcs(1))
})

test_that("knots is numeric, has length >= 3, and is in the correct order ", {
  expect_error(rcs(1:10, 1))
  expect_error(rcs(1:10, 3:1))
})

test_that("warning on duplicate knots", {
  expect_warning(rcs(1:10, c(1,1,2,3)))
})
