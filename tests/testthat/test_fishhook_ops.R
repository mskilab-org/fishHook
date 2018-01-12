library(fishhook)

context('test fishhook operations')











test_that('qq_pval', {
	pvals = c(0.0001, 0.0001, 0.032, 0.005, 0.9, 0.15)
	## check default params
	foo = qq_pval(pvals)
    expect_match(names(foo)[1], 'rect')
    expect_match(names(foo)[2], 'text')
    expect_equal(round(foo$rect$w, 3), round(0.8597488, 3))  ## rounding; otherwise, 'double' types not exactly the same
    expect_equal(round(foo$rect$h, 3), round(0.4891113, 3))
    expect_equal(round(foo$rect$left, 3), round(3.820251, 3))
    expect_equal(round(foo$rect$top, 3), round(0.3091113, 3))
    expect_equal(round(foo$text$x, 3), round(4.073376, 3))
    expect_equal(round(foo$text$y, 3), round(0.008372093, 3))
    ## exp: Error: length of exp must be = length(obs
    expect_error(qq_pval(pvals, exp=(c(1, 2, 3))))
    
})
