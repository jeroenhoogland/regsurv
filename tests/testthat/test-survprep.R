test_that("expect output of class survprep for loghazard models", {
  prep <- survprep(tte=simdata$eventtime,
                   delta=simdata$status,
                   X=as.matrix(simdata[ ,grep("x", names(simdata))]),
                   model.scale="loghazard",
                   time.scale="time",
                   spline.type="rcs",
                   ntimebasis=4,
                   qpoints=9)
  expect_s3_class(prep, "survprep")
  expect_equal(dim(prep$mm.scaled$d), c(250, 14))
})

test_that("expect output of class survprep for logHazard models", {
  prep <- survprep(tte=simdata$eventtime,
                   delta=simdata$status,
                   X=as.matrix(simdata[ ,grep("x", names(simdata))]),
                   model.scale="logHazard",
                   time.scale="logtime",
                   spline.type="rcs",
                   ntimebasis=4,
                   tv=1,nitimebasis=2,
                   qpoints=9)
  expect_s3_class(prep, "survprep")
  expect_equal(prep$mm.scaled$nibasis, 4)
})

test_that("expect output of class survprep for logHazard models", {
  prep <- survprep(tte=simdata$eventtime,
                   delta=simdata$status,
                   X=as.matrix(simdata[ ,grep("x", names(simdata))]),
                   model.scale="logHazard",
                   time.scale="logtime",
                   spline.type="rcs",
                   ntimebasis=4,
                   tv=1:9,
                   nitimebasis=1,
                   qpoints=9)
  expect_s3_class(prep, "survprep")
  expect_equal(prep$mm.scaled$nibasis, 9)
})

test_that("expect output of class survprep for logHazard models", {
  prep <- survprep(tte=simdata$eventtime,
                   delta=simdata$status,
                   X=as.matrix(simdata[ ,grep("x", names(simdata))]),
                   model.scale="logHazard",
                   time.scale="logtime",
                   spline.type="rcs",
                   ntimebasis=1,
                   tv=1,
                   nitimebasis=1,
                   qpoints=9)
  expect_s3_class(prep, "survprep")
  expect_equal(prep$mm.scaled$nbasis, 1)
  expect_equal(prep$mm.scaled$nibasis, 1)
})

test_that("expect output of class survprep for logHazard models", {
  prep <- survprep(tte=simdata$eventtime,
                   delta=simdata$status,
                   X=as.matrix(simdata[ ,grep("x", names(simdata))]),
                   model.scale="logHazard",
                   time.scale="logtime",
                   spline.type="rcs",
                   ntimebasis=7,
                   tv=1:9,
                   nitimebasis=0,
                   qpoints=9)
  expect_s3_class(prep, "survprep")
  expect_equal(prep$mm.scaled$nbasis, 7)
  expect_equal(prep$mm.scaled$nibasis, 0)
})
