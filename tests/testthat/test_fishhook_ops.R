
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

eligible = readRDS('/home/travis/build/mskilab/fishHook/data/eligible.rds')

annocov = readRDS('/home/travis/build/mskilab/fishHook/data/anno_covs.rds')



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
    ## expect_equal(length(annotate.targets('/home/travis/build/mskilab/fishHook/data/targets.rds')), 19688)
    ## if (length(targets)==0){
    expect_error(annotate.targets(GRanges()))
    ## if (!is.null(out.path)){
    ## expect_equal(length(annotate.targets('/home/travis/build/mskilab/fishHook/data/targets.rds', out.path = '/home/travis/build/mskilab/fishHook/data/output.RDS')), 19688)
    ## if (!all(cov.types %in% COV.TYPES) & !(all(cov.classes %in% COV.CLASSES))){
    expect_error(annotate.targets(targets, covariates = grl2))
    ## if events != NULL
    expect_equal(max(annotate.targets(targets, events=events)$count), 9750)
    ##
    annotate.targets('/home/travis/build/mskilab/fishHook/data/targets.rds')
    expect_true(is(annotate.targets(targets), 'GRanges')) 

})



## aggregate.targets

test_that('aggregate.targets', {

    expect_error(aggregate.targets(targets))  ## Error: argument "by" must be specified and same length as targets or "rolling" must be non NULL
    foo = aggregate.targets(targets, by='gene_name')
    expect_equal(length(foo$gene_name), 16352)
    ##  if (is.null(by) & is.character(targets)){
    ## expect_error(aggregate.targets('/home/travis/build/mskilab/fishHook/data/targets.rds'))
    ## annotate targets
    annotated = annotate.targets(targets, events=events)
    ## by
    expect_equal(length(aggregate.targets(annotated, by = 'gene_name')[[1]]), 16352)
    ## fields
    expect_equal(length(aggregate.targets(annotated, fields = 'count', rolling = 200)), 15032)
    ## rolling
    #### Error: Malformed input, only NA ranges produced. Reduce value of running
    expect_error(aggregate.targets(annotated, rolling = 5e5))
    #### Column 1 of result for group 2 is type 'logical' but expecting type 'double'. Column types must be consistent for each group.
    expect_error(aggregate.targets(annotated, rolling = 1e3))
    #### works
    expect_equal(length(aggregate.targets(annotated, rolling = 200)), 15032)
    ## disjoint
    expect_equal(length(aggregate.targets(annotated, rolling = 200, disjoint = FALSE)), 15032)
    ## na.rm
    expect_equal(length(aggregate.targets(annotated, rolling = 200, na.rm = TRUE)), 15032)   
    ## FUN 
    ## verbose
    expect_equal(length(aggregate.targets(annotated, by = 'gene_name', verbose = FALSE)[[1]]), 16352)
    ## 
    ##  if (is.null(by) & is.character(targets)){
    expect_error(aggregate.targets('/home/travis/build/mskilab/fishHook/data/targets.rds'))  ## Coverage missing for input targets
    ##  if (is.null(by) & is.character(targets)){ (continued)
    expect_error(aggregate.targets('/home/travis/build/mskilab/fishHook/data/annotated_cov.rds'))


})






## score.targets

test_that('score.targets', {

    annotated = annotate.targets(targets, events=events)
    ## targets
    ## covariates = names(values(targets))
    ## model = NULL
    ## return.model = FALSE
    expect_equal(as.integer(round(score.targets(annotated, return.model = TRUE)$coefficients)), as.integer(-5))
    ## nb = TRUE
    expect_equal(round(max(score.targets(annotated, nb=TRUE)$effectsize)), 5)
    expect_equal(round(max(score.targets(annotated, nb=FALSE)$effectsize)), 7)
    ## verbose = TRUE
    expect_equal(dim(as.data.table(score.targets(annotated, verbose=FALSE)))[1], 19688)
    expect_equal(dim(as.data.table(score.targets(annotated, verbose=FALSE)))[2], 16)
    ## iter = 200
    expect_equal(dim(as.data.table(score.targets(annotated, iter=4000)))[1], 19688)
    expect_equal(dim(as.data.table(score.targets(annotated, iter=4000)))[2], 16)
    ## subsample = 1e5
    expect_equal(dim(as.data.table(score.targets(annotated, subsample = 5)))[1], 19688)
    expect_equal(dim(as.data.table(score.targets(annotated, subsample = 5)))[2], 16)
    ## seed = 42
    ## p.randomized = TRUE
    expect_equal(dim(as.data.table(score.targets(annotated, p.randomized = FALSE)))[1], 19688)
    ## classReturn = FALSE
    expect_equal(dim(as.data.table(score.targets(annotated, classReturn = TRUE)))[1], 19688)

})





## Cov()

test_that('Cov', {

    expect_error(Cov())
    foo = Cov$new(Covariate=replication_timing, type='numeric')
    ## seqlevels
    expect_equal(length(foo$seqlevels()), 25)
    ## toString
    expect_match(foo$toString(),"Name: \ntype: numeric\tsignature: NA\nfield: NA\tpad: NA\nna.rm: NA\tgrep: NA\nCovariate: GRanges\n")
    ## convert2Arr
    expect_equal(foo$convert2Arr()$type, 'numeric') ## outputs whatever the user above input
    ## print
    expect_equal(print(foo$print()), NULL)   
    ## chr
    expect_equal(foo$chr(), FALSE)
    ## toList
    expect_equal(length(foo$toList()$track), 2385966)
    expect_equal(foo$toList()$type, 'numeric')
    expect_equal(foo$toList()$signature, NA)
    expect_equal(foo$toList()$pad, NA)
    expect_equal(foo$toList()$na.rm, NA)
    expect_equal(foo$toList()$field, NA)
    expect_equal(foo$toList()$grep, NA)
    ## 
    ## if(is.null(Covariate) | is.null(type)){
    expect_error(Cov$new())
    ## 

})



## Cov_Arr

test_that('Cov_Arr', {

    ## initialize = function(..., name = NULL, cvs = NULL, pad = NULL, type = NULL, signature = NULL, field = NULL, na.rm = NULL, grep = NULL)
    ## default
    foobar = Cov_Arr$new(Cov$new(Covariate=replication_timing, type='numeric'))
    foobar2 = Cov_Arr$new(Cov$new(Covariate=replication_timing[50:100], type='numeric'))
    ## merge
    #### check it runs
    expect_error(foobar$merge(foobar2), NA)
    ## chr
    expect_false(foobar$chr())
    ## seqlevels
    expect_equal(length(foobar$seqlevels()[[1]]), 25)
    ## subset
    expect_equal(length(foobar$subset(1:4)), 4)
    ## toList
    expect_equal(length(foobar$toList()[[1]]$track), 2385966)
    expect_equal(foobar$toList()[[1]]$type, 'numeric')
    expect_equal(foobar$toList()[[1]]$signature, NA)
    expect_equal(foobar$toList()[[1]]$pad, NA)
    expect_equal(foobar$toList()[[1]]$na.rm, NA)
    expect_equal(foobar$toList()[[1]]$field, NA)
    expect_equal(foobar$toList()[[1]]$grep, NA)
    ## print
    expect_equal(foobar$print()[[1]], NULL)

})




## FishHook

## initialize = function(targets = NULL, out.path = NULL, eligible = NULL, ... ,events = NULL, covariates = NULL,
##     use_local_mut_density = FALSE, local_mut_density_bin = 1e6, genome = 'BSgenome.Hsapiens.UCSC.hg19::Hsapiens',
##     mc.cores = 1, na.rm = TRUE, pad = 0, verbose = TRUE, max.slice = 1e3, ff.chunk = 1e6, max.chunk = 1e11, ptidcol = NULL,
##     maxpatientpergene = Inf, weightEvents = FALSE, nb = TRUE)

test_that('FishHook', {

    ## default
    fish1 = FishHook$new(targets = targets, events = events)
    expect_true(fish1$nb)
    expect_equal(fish1$state, 'Initialized')
    expect_false(fish1$weightEvents)
    anno = fish1$annotate(mc.cores=1);
    print('created anno')
    ## toList
    ## > fish1$toList()
    ## Error in fish1$toList() : argument "x" is missing, with no default
    ## print
    expect_equal(print(fish1$print()), NULL)
    ## annotate
    expect_equal(anno, "Annotated")
    ## aggregate
    ## > agg = fish1$aggregate(mc.cores=2)
    ## Error in fish1$aggregate(mc.cores = 2) : unused argument (mc.cores = 2)
    expect_error(fish1$aggregate())
    agg = fish1$aggregate(by='gene_name')
    expect_equal(print(agg), "Clear Completed")
    ## score
    ##
    ## Error in score.targets(targ, covariates = names(values(private$panno)),  : 
    ## Error: "score.targets" input malformed --> count does not vary!
    ## qq_plot
    expect_error(fish1$qq_plot())
    ## clear
    expect_equal(print(fish1$clear()), "Clear Completed")
    ## with eligible
    ## foobar = FishHook$new(targets = targets, events = events, eligible = eligible)
    fish2 = FishHook$new(targets = targets, events = events)
    ## test active bindings
    ## cvs
    expect_equal(fish2$csv, NULL)
    ## eligible
    expect_equal(fish2$eligible, NULL)
    anno2 = fish2$annotate(mc.cores=1);
    print('created anno2')
    ## anno
    expect_equal(max(fish2$anno$count), 9750)
    ## targets 
    expect_equal(length(fish2$targets), 19688)



})


### Annotate

## anno= Annotated$new(targets=targets, events=events, covariates=replication_timing)
## SLOW



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
    foobar = qq_pval(pvals, plotly=TRUE)
    expect_match(names(foobar)[1], 'x')
    expect_match(names(foobar)[2], 'width')
    expect_match(names(foobar)[3], 'height')

})



