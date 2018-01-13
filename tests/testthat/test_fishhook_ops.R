library(fishhook)


Sys.setenv(DEFAULT_BSGENOME = 'BSgenome.Hsapiens.UCSC.hg19::Hsapiens')

# Sample Events
events = readRDS('/home/travis/build/mskilab/fish.hook/data/events.rds')

# Sample Targets
targets = readRDS('/home/travis/build/mskilab/fish.hook/data/targets.rds')

# Sample Covariate
replication_timing = readRDS('/home/travis/build/mskilab/fish.hook/data/covariate.rds')

# Same Eligible Subset
eligible = readRDS('/home/travis/build/mskilab/fish.hook/data/eligible.rds')

# indexed pathways
indexed_pathways = readRDS('/home/travis/build/mskilab/fish.hook/data/indexed_pathways.rds')

segs = readRDS('/home/travis/build/mskilab/fish.hook/data/jabba_segs_11517.rds')


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
    ## not sure how to test 'lwd'
    ## or these other args
})

### from Zoran's tests

test_that('Cov', {
	rept = Cov$new(Covariate = replication_timing, type = 'numeric', name = 'rept')
	expect_equal(length(names(rept)), 19)
	expect_match(rept$type, 'numeric')
	expect_match(rept$name, 'rept')
	expect_true(inherits(rept$Covariate, 'GRanges'))
	expect_equal(length(rept$Covariate), 2385966)
    Covariates = c(rept, rept, rept, rept, rept)
    x = Covariates[2:4]
    y = x$merge(x, x)
    expect_equal(length(x$type), 3)
    test1 = x$toList()
    expect_equal(length(test1), 3)
    z = y$toList()
    expect_equal(length(z), 9)
    expect_true(all(x$names == rep('rept', 3)))
    x2 = x[1]$merge(x[2:3], x)
    x3 = c(x[1], x[2:3], x) 
    x4 = x3  
    expect_false(all(x4$chr()))
    expect_equal(length(x4$seqlevels()[[1]]), 25)
    r = rept$convert2Arr()
    r1 = c(rept, x, rept, rept, rept, x)
    r2 = c(x, rept, x, x, x, rept)
    targets$pathways = NULL
})


test_that('FishHook', {
	fish1 = FishHook$new(targets = targets, events = events, eligible = eligible)
	expect_equal(length(names(fish1)), 32)
	expect_error(fish1$annotate(), NA)  ### one way to check that an error does NOT occur
	## testing patient ID
	events$id = events$patient_code
	events$patient_code = NULL
	fish2 = FishHook$new(targets = targets, events = events, eligible = eligible)
	fish2$annotate(maxpatientpergene = 1, ptidcol = "id")
	fish2$score()
	testList = GRangesList(events[1:100], events[101:200], events[201:300])
	expect_equal(fish2$aggregated, NULL)
    fish2$aggregated = testList
    expect_equal(length(fish2$aggregated[[1]]), 100)
    expect_true(inherits(fish2$aggregated[[2]] , 'GRanges'))
	## testing aggregation
    tiles = gr.tile(hg_seqlengths(), 100000)
    expect_equal(nrow(segs), 134231)
    colnames(segs) = c('ID', 'chr', 'start', 'end', 'd', 'cn')
    segs = segs[cn >= 1.2]
    expect_equal(nrow(segs), 12373)
    exptect_true(inherits(segs, 'data.table'))
    segs = dt2gr(segs)
    exptect_true(inherits(segs, 'GRanges'))
})




