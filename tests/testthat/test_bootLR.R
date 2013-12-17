# Tests for package bootLR


# Options
# expect_that(x, is_true())   expect_true(x)
# expect_that(x, is_false()) 	expect_false(x)
# expect_that(x, is_a(y)) 	expect_is(x, y)
# expect_that(x, equals(y)) 	expect_equal(x, y)
# expect_that(x, is_equivalent_to(y)) 	expect_equivalent(x, y)
# expect_that(x, is_identical_to(y)) 	expect_identical(x, y)
# expect_that(x, matches(y)) 	expect_match(x, y)
# expect_that(x, prints_text(y)) 	expect_output(x, y)
# expect_that(x, shows_message(y)) 	expect_message(x, y)
# expect_that(x, gives_warning(y)) 	expect_warning(x, y)
# expect_that(x, throws_error(y)) 	expect_error(x, y)

context( "Output is as expected" )

test_that( "Sens 100/100 and Spec 60/100", {
  
} )

context( "Inputs are validated" )

test_that( "totalNeg or totalPos == 0", {
  expect_error( BayesianLR.test( truePos=98, totalPos=0, trueNeg=60, totalNeg=100 ) )
  expect_error( BayesianLR.test( truePos=98, totalPos=100, trueNeg=60, totalNeg=0 ) )
} )


test_that( "trueNeg or truePos greater than than totalNeg or totalPos", {
  expect_error( BayesianLR.test( truePos=98, totalPos=95, trueNeg=60, totalNeg=100 ) )
  expect_error( BayesianLR.test( truePos=98, totalPos=100, trueNeg=60, totalNeg=55 ) )
} )


test_that("dummy test", { expect_true(FALSE) } )