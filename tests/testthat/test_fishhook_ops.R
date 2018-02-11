
library(fishHook)
library(testthat)

Sys.setenv(DEFAULT_BSGENOME = 'BSgenome.Hsapiens.UCSC.hg19::Hsapiens')

# Sample Events
events = readRDS('/home/travis/build/mskilab/fishHook/data/events.rds')
## events = readRDS('events.rds')

# Sample Targets
targets = readRDS('/home/travis/build/mskilab/fishHook/data/targets.rds')
## targets = readRDS('targets.rds')


# Sample Covariate
replication_timing = readRDS('/home/travis/build/mskilab/fishHook/data/covariate.rds')
## replication_timing = readRDS('covariate.rds')


# Same Eligible Subset
eligible = readRDS('/home/travis/build/mskilab/fishHook/data/eligible.rds')
## eligible  = readRDS('eligible.rds')


# indexed pathways
indexed_pathways = readRDS('/home/travis/build/mskilab/fishHook/data/indexed_pathways.rds')
## indexed_pathways = readRDS('indexed_pathways.rds')


segs = readRDS('/home/travis/build/mskilab/fishHook/data/jabba_segs_11517.rds')
## segs = readRDS('jabba_segs_11517.rds')




context('unit testing fishhook operations')



## annotate.targets = function(targets, covered = NULL, events = NULL,  mc.cores = 1, na.rm = TRUE, pad = 0, verbose = TRUE, max.slice = 1e3, 
##    ff.chunk = 1e6, max.chunk = 1e11, out.path = NULL, covariates = list(), maxpatientpergene = Inf, ptidcol = NULL, weightEvents = FALSE, ...)

test_that('annotate.targets', {
	## default args
	expect_true(is(annotate.targets(targets), 'GRanges'))   ## annotate.targets(targets, weightEvents=TRUE)
    expect_equal(length(annotate.targets(targets)), 19688)
    foo = annotate.targets(targets)
    expect_equal(max(foo$query.id), 19688)
    expect_equal(max(foo$coverage), 2304638)
    expect_equal(max(foo$count), 0)
    ## covered
    ## events
    ## mc.cores
    expect_true(is(annotate.targets(targets, mc.cores=2), 'GRanges'))
    expect_equal(length(annotate.targets(targets, mc.cores=2)), 19688)
    expect_equal(max(annotate.targets(targets, mc.cores=2)$query.id), 19688)
    expect_equal(max(annotate.targets(targets, mc.cores=2)$coverage), 2304638)
    expect_equal(max(annotate.targets(targets, mc.cores=2)$count), 0)
    ## na.rm
    ## pad
    ## verbose
    expect_true(is(annotate.targets(targets, verbose=FALSE), 'GRanges'))
    expect_equal(length(annotate.targets(targets, verbose=FALSE)), 19688)
    expect_equal(max(annotate.targets(targets, verbose=FALSE)$query.id), 19688)
    expect_equal(max(annotate.targets(targets, verbose=FALSE)$coverage), 2304638)
    expect_equal(max(annotate.targets(targets, verbose=FALSE)$count), 0)
    ## max.slice
    ## annotate.targets(targets, max.slice=1e8)
    ## ff.chunk
    ## 
    ## max.chunk
    ## out.path
    ## covariates
    ## maxpatientpergene
    ## ptidcol
    ## weightEvents
    expect_equal(is.na(all(annotate.targets(targets, weightEvents=TRUE)$count)), TRUE)
    ## if (is.character(targets)){
    expect_equal(length(annotate.targets('/home/travis/build/mskilab/fishHook/data/targets.rds')), 19688)
    ## if (length(targets)==0){
    expect_error(annotate.targets(GRanges()))
    ## if (!is.null(out.path)){
    expect_equal(length(annotate.targets('/home/travis/build/mskilab/fishHook/data/targets.rds', out.path = '/home/travis/build/mskilab/fishHook/data/output.RDS')), 19688)
    ## if (!all(cov.types %in% COV.TYPES) & !(all(cov.classes %in% COV.CLASSES))){
    expect_error(annotate.targets(targets, covariates = grl2))
    ## if events != NULL
    expect_equal(max(annotate.targets(targets, events=events)$count), 9750)


})



## aggregate.targets

test_that('aggregate.targets', {

    expect_error(aggregate.targets(targets))  ## Error: argument "by" must be specified and same length as targets or "rolling" must be non NULL
    foo = aggregate.targets(targets, by='gene_name')
    expect_equal(length(foo$gene_name), 16352)
    ##  if (is.null(by) & is.character(targets)){
    expect_error(aggregate.targets('/home/travis/build/mskilab/fishHook/data/targets.rds'))

})












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

##test_that('Cov', {
##	rept = Cov$new(Covariate = replication_timing, type = 'numeric', name = 'rept')
##	expect_equal(length(names(rept)), 19)
##	expect_match(rept$type, 'numeric')
##	expect_match(rept$name, 'rept')
##	expect_true(inherits(rept$Covariate, 'GRanges'))
##	expect_equal(length(rept$Covariate), 2385966)
##    Covariates = c(rept, rept, rept, rept, rept)
##   x = Covariates[2:4]
##    y = x$merge(x, x)
##    expect_equal(length(x$type), 3)
##    test1 = x$toList()
##    expect_equal(length(test1), 3)
##    z = y$toList()
##    expect_equal(length(z), 9)
##    expect_true(all(x$names == rep('rept', 3)))
##    x2 = x[1]$merge(x[2:3], x)
##    x3 = c(x[1], x[2:3], x) 
##    x4 = x3  
##    expect_false(all(x4$chr()))
##    expect_equal(length(x4$seqlevels()[[1]]), 25)
##    r = rept$convert2Arr()
##    r1 = c(rept, x, rept, rept, rept, x)
##    r2 = c(x, rept, x, x, x, rept)
##    targets$pathways = NULL
##})


#test_that('FishHook', {
#	fish1 = FishHook$new(targets = targets, events = events, eligible = eligible)
##	expect_equal(length(names(fish1)), 32)
##	expect_error(fish1$annotate(), NA)  ### one way to check that an error does NOT occur
##	## testing patient ID
##	events$id = events$patient_code
##	events$patient_code = NULL
##	fish2 = FishHook$new(targets = targets, events = events, eligible = eligible)
##	fish2$annotate(maxpatientpergene = 1, ptidcol = "id")
##	fish2$score()
##	testList = GRangesList(events[1:100], events[101:200], events[201:300])
##	expect_equal(fish2$aggregated, NULL)
##   fish2$aggregated = testList
##    expect_equal(length(fish2$aggregated[[1]]), 100)
##    expect_true(inherits(fish2$aggregated[[2]] , 'GRanges'))
##	## testing aggregation
##    tiles = gr.tile(hg_seqlengths(), 100000)
##    expect_equal(nrow(segs), 134231)
##    colnames(segs) = c('ID', 'chr', 'start', 'end', 'd', 'cn')
##    segs = segs[cn >= 1.2]
##    expect_equal(nrow(segs), 12373)
##    expect_true(inherits(segs, 'data.table'))
##    segs = dt2gr(segs)
##    expect_true(inherits(segs, 'GRanges'))
##})




