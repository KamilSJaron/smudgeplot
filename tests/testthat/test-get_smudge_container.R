context("smudge container generation")

x <- c(2,2,2,3,3,3,3,3,3,3)
y <- c(1,1,2,2,2,2,3,3,3,3)

cont <- get_smudge_container(x, y, .nbins = 3, .xlim = c(0, 3), .ylim = c(0, 3))
exp <- matrix(c(0,0,0,2,1,0,0,3,4), ncol = 3, byrow = T)

test_that("generation of smudge matrix categorize as expected", {
        expect_true(all(cont$dens == exp))
    }
)

minor_cov <- c(0.5, 0.33, 0.32, 0.25, 0.125, 0.31)
pair_cov <- c(40, 60, 61, 83, 164, 123)
cont <- get_smudge_container(minor_cov, pair_cov, .nbins = 4)
exp <- matrix(c(0,0,0,1, 0,0,1,0, 0,2,1,0, 1,0,0,0), ncol = 4, byrow = T)

test_that("generation of smudge matrix categorize as expected", {
        expect_true(all(cont$dens == exp))
    }
)
