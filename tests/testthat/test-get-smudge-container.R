context("core")

nbins <- 3
x.bin <- seq(, length=nbins)
y.bin <- seq(, length=nbins)
x <- c(2,2,2,3,3,3,3,3,3,3)
y <- c(1,1,2,2,2,2,3,3,3,3)

cont <- get_smudge_container(x, y, .nbins = 3, .xlim = c(0, 3), .ylim = c(0, 3))
exp <- matrix(c(0,0,0,2,1,0,0,3,4), ncol = 3, byrow = T)

test_that("generation of smudge matrix categorize as expected", {
        expect_true(all(cont$dens == exp))
    }
)

