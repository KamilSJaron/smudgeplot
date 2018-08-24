context("utilities")

test_that("transform_x project corectly coordinates", {
        expect_that(transform_x(1, F), equals(0))
        expect_that(transform_x(30), equals(0.48))
        expect_that(transform_x(30, F), equals(0.5))
    }
)