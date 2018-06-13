library(fishHook)
library(testthat)
library(zoo)

Sys.setenv(DEFAULT_BSGENOME = 'BSgenome.Hsapiens.UCSC.hg19::Hsapiens')

## hypotheses BED
#hypothesesbed = system.file("extdata", "targets.bed", package = 'fishHook')
hypothesesbed = system.file("data", "targets.bed", package = 'fishHook')

## Sample Covariate
#replication_timing = readRDS(system.file("extdata", "covariate.rds", package = 'fishHook'))
replication_timing = readRDS(system.file("data", "covariate.rds", package = 'fishHook'))

## Same Eligible Subset
## eligible = readRDS(system.file("extdata", "eligible.rds", package = 'fishHook'))

## Sample Hypotheses
#hypotheses = readRDS(system.file("extdata", "targets.rds", package = 'fishHook'))
hypotheses  = readRDS(system.file("data", "targets.rds", package = 'fishHook'))


## indexed pathways
#indexed_pathways = readRDS(system.file("extdata", "indexed_pathways.rds", package = 'fishHook'))
indexed_pathways = readRDS(system.file("data", "indexed_pathways.rds", package = 'fishHook'))

## segs
segs = readRDS(system.file("data", "jabba_segs_11517.rds", package = 'fishHook'))
#segs = readRDS(system.file("extdata", "jabba_segs_11517.rds", package = 'fishHook'))

## eligible
## eligible = readRDS(system.file("extdata", "eligible.rds", package = 'fishHook'))

# Sample annotate
#anno = readRDS(system.file("extdata", "anno.rds", package = 'fishHook'))
anno = readRDS(system.file("data", "anno.rds", package = 'fishHook'))

context('unit testing fishhook operations')


## annotate.hypotheses = function(hypotheses, covered = NULL, events = NULL,  mc.cores = 1, na.rm = TRUE, pad = 0, verbose = TRUE, max.slice = 1e3, 
##    ff.chunk = 1e6, max.chunk = 1e11, out.path = NULL, covariates = list(), maxpatientpergene = Inf, ptidcol = NULL, weightEvents = FALSE, ...)

test_that('annotate.hypotheses', {

    ## default args
    t1 = hypotheses[1]
    e1 = events[1]
    el1 = eligible[1]
    r1 = replication_timing[1]
    expect_error({annotate.hypotheses(hypotheses = hypotheses[1], events = events[1], idcol = 'hello')})
    c1 = Covariate$new(data = system.file("data", "covariate.rds", package = 'fishHook'), type = 'numeric', name = 'c1', pad = 10)
    c2 = Covariate$new(data = system.file("data", "covariate.rds", package = 'fishHook'), type = 'interval', name = 'c2', pad = 10, na.rm = TRUE)
    c3 = Covariate$new(data = list(RleList()), type = 'interval', name = 'c3', pad = 10, na.rm = TRUE)
#    c1 = Covariate$new(data = list(system.file("extdata", "covariate.rds", package = 'fishHook')), type = 'numeric', name = 'c1', pad = 10)
#   c2 = Covariate$new(data = list(system.file("extdata", "covariate.rds", package = 'fishHook')), type = 'interval', name = 'c2', pad = 10, na.rm = TRUE)
#    c3 = Covariate$new(data = list(RleList()), type = 'interval', name = 'c3', pad = 10, na.rm = TRUE)
    anno = annotate.hypotheses(hypotheses = t1, events = e1, covariates = c1$toList())
    anno2 = annotate.hypotheses(hypotheses = t1, events = e1, covariates = c2$toList())
    #anno3 = annotate.hypotheses(hypotheses = t1, events = e1, covariates = c3$toList())
    anno4 = annotate.hypotheses(hypotheses = e1, events = e1, eligible = e1)
    expect_true(is(annotate.hypotheses(hypotheses), 'GRanges'))   ## annotate.hypotheses(hypotheses, weightEvents=TRUE)
    expect_equal(length(annotate.hypotheses(hypotheses)), 19688)
    foo = annotate.hypotheses(hypotheses)
    expect_equal(max(foo$query.id), 19688)
    expect_equal(max(foo$coverage), 2304638)
    expect_equal(max(foo$count), 0)
    ## covered
    ## events
    ## mc.cores
    expect_true(is(annotate.hypotheses(hypotheses, mc.cores=2), 'GRanges'))
    expect_equal(length(annotate.hypotheses(hypotheses, mc.cores=2)), 19688)
    expect_equal(max(annotate.hypotheses(hypotheses, mc.cores=2)$query.id), 19688)
    expect_equal(max(annotate.hypotheses(hypotheses, mc.cores=2)$coverage), 2304638)
    expect_equal(max(annotate.hypotheses(hypotheses, mc.cores=2)$count), 0)
    ## na.rm
    ## pad
    ## verbose
    expect_true(is(annotate.hypotheses(hypotheses, verbose=FALSE), 'GRanges'))
    expect_equal(length(annotate.hypotheses(hypotheses, verbose=FALSE)), 19688)
    expect_equal(max(annotate.hypotheses(hypotheses, verbose=FALSE)$query.id), 19688)
    expect_equal(max(annotate.hypotheses(hypotheses, verbose=FALSE)$coverage), 2304638)
    expect_equal(max(annotate.hypotheses(hypotheses, verbose=FALSE)$count), 0)
    ## max.slice
    ## annotate.hypotheses(hypotheses, max.slice=1e8)
    ## ff.chunk
    ## 
    ## max.chunk
    ## out.path
    ## covariates
    ## maxpatientpergene
    ## ptidcol
    ## weightEvents
    expect_equal(is.na(all(annotate.hypotheses(hypotheses, weightEvents=TRUE)$count)), TRUE)
    ## if (is.character(hypotheses)){
    ## expect_equal(length(annotate.hypotheses('/home/travis/build/mskilab/fishHook/data/hypotheses.rds')), 19688)
    ## if (length(hypotheses)==0){
    expect_error(annotate.hypotheses(GRanges()))
    ## if (!is.null(out.path)){
    ## expect_equal(length(annotate.hypotheses('/home/travis/build/mskilab/fishHook/data/hypotheses.rds', out.path = '/home/travis/build/mskilab/fishHook/data/output.RDS')), 19688)
    ## if (!all(cov.types %in% COV.TYPES) & !(all(cov.classes %in% COV.CLASSES))){
    expect_error(annotate.hypotheses(hypotheses, covariates = grl2))
    ## if events != NULL
    expect_equal(max(annotate.hypotheses(hypotheses, events=events)$count), 9750)
    ##
    annotate.hypotheses(system.file("data", "targets.rds", package = 'fishHook'))
    expect_true(is(annotate.hypotheses(hypotheses), 'GRanges'))
    ##Interval Covariates
    int = Covariate$new(data = replication_timing[1], name = 'int', type = 'interval', pad = 10)
    fish3 = FishHook$new(hypotheses = hypotheses, events = events[1], covariates = int)
    fish3$score()
    ## 
    ## else if (grepl('(\\.bed$)', hypotheses[1])){
    expect_equal(length(annotate.hypotheses(hypothesesbed)), 19688)
})



## aggregate.hypotheses

test_that('aggregate.hypotheses', {

    expect_error(aggregate.hypotheses(system.file("data", "targets.rds", package = 'fishHook'), rolling = 1))
    agg = aggregate.hypotheses(c(system.file("data", "anno.rds", package = 'fishHook'), system.file("data", "anno.rds", package = 'fishHook')), rolling = 1, na.rm = T)
    expect_error({agg = aggregate.hypotheses(anno, rolling = -1)})
    expect_error(aggregate.hypotheses(hypotheses))  ## Error: argument "by" must be specified and same length as hypotheses or "rolling" must be non NULL
    foo = aggregate.hypotheses(hypotheses, by='gene_name')
    expect_equal(length(foo$gene_name), 16352)
    ##  if (is.null(by) & is.character(hypotheses)){
    ## expect_error(aggregate.hypotheses('/home/travis/build/mskilab/fishHook/data/hypotheses.rds'))
    ## annotate hypotheses
    annotated = annotate.hypotheses(hypotheses, events=events)
    ## by
    expect_equal(length(aggregate.hypotheses(annotated, by = 'gene_name')[[1]]), 16352)
    ## fields
    expect_equal(length(aggregate.hypotheses(annotated, fields = 'count', rolling = 200)), 15032)
    ## rolling
    #### Error: Malformed input, only NA ranges produced. Reduce value of running
    expect_error(aggregate.hypotheses(annotated, rolling = 5e5))
    #### Column 1 of result for group 2 is type 'logical' but expecting type 'double'. Column types must be consistent for each group.
    expect_error(aggregate.hypotheses(annotated, rolling = 1e3))
    #### works
    expect_equal(length(aggregate.hypotheses(annotated, rolling = 200)), 15032)
    ## disjoint
    expect_equal(length(aggregate.hypotheses(annotated, rolling = 200, disjoint = FALSE)), 15032)
    ## na.rm
    expect_equal(length(aggregate.hypotheses(annotated, rolling = 200, na.rm = TRUE)), 15032)   
    ## FUN 
    ## verbose
    expect_equal(length(aggregate.hypotheses(annotated, by = 'gene_name', verbose = FALSE)[[1]]), 16352)
    ## 
    ##  if (is.null(by) & is.character(hypotheses)){
    expect_error(aggregate.hypotheses(system.file("data", "hypotheses.rds", package = 'fishHook')))  ## Coverage missing for input hypotheses
    ##  if (is.null(by) & is.character(hypotheses)){ (continued)
    expect_error(aggregate.hypotheses(system.file("data", "annotated_cov.rds", package = 'fishHook')))
    ##testing rolling aggregation
    start = c(1,1001,2001,3001,4001,5001)
    end = c(1000,2000,3000,4000,5000,6000)
    chr = 1
    test = dt2gr(as.data.table(cbind(start,end,chr)))
    mcols(test) = NULL
    fishAgg = FishHook$new(hypotheses = test, events = test)
    fishAgg$aggregate(rolling = 3)
    expect_error(fishAgg$aggregate(rolling = 0))
    expect_error(fishAgg$aggregate(rolling = 4.5))
    expect_equal(length(fishAgg$aggregated), 4)
    expect_equal(any(!width(fishAgg$aggregated) == rep(3000,4)), FALSE)
})






## score.hypotheses

test_that('score.hypotheses', {

    annotated = annotate.hypotheses(hypotheses, events=events)
    ## hypotheses
    ## covariates = names(values(hypotheses))
    ## model = NULL
    ## return.model = FALSE
    expect_equal(as.integer(round(score.hypotheses(annotated, return.model = TRUE)$coefficients)), as.integer(-5))
    ## nb = TRUE
    expect_equal(round(max(score.hypotheses(annotated, nb=TRUE)$effectsize)), 5)
    expect_equal(round(max(score.hypotheses(annotated, nb=FALSE)$effectsize)), 7)
    ## verbose = TRUE
    expect_equal(dim(as.data.table(score.hypotheses(annotated, verbose=FALSE)))[1], 19688)
    expect_equal(dim(as.data.table(score.hypotheses(annotated, verbose=FALSE)))[2], 16)
    ## iter = 200
    expect_equal(dim(as.data.table(score.hypotheses(annotated, iter=4000)))[1], 19688)
    expect_equal(dim(as.data.table(score.hypotheses(annotated, iter=4000)))[2], 16)
    ## subsample = 1e5
    expect_equal(dim(as.data.table(score.hypotheses(annotated, subsample = 5)))[1], 19688)
    expect_equal(dim(as.data.table(score.hypotheses(annotated, subsample = 5)))[2], 16)
    ## seed = 42
    ## p.randomized = TRUE
    expect_equal(dim(as.data.table(score.hypotheses(annotated, p.randomized = FALSE)))[1], 19688)
    ## classReturn = FALSE
    expect_equal(dim(as.data.table(score.hypotheses(annotated, classReturn = TRUE)$res))[1], 19688)

    ## MISC testing
    anno = annotated
    expect_error({s1 = score.hypotheses(anno, covariates = c('hello'))})
    anno$cov2 = c(1:length(anno))
    anno$cov2 = as.character(anno$cov2)
    expect_error({s2 = score.hypotheses(anno)})
    anno3 = anno
    anno3$count = NULL
    anno3$cov2 = NULL
    expect_error({s3 = score.hypotheses(anno3)})
    anno3$count = 1
    expect_error({s4 = score.hypotheses(anno3)})
    anno3$dummy = NA
    anno3$dummy = as.numeric(anno3$dummy)
    anno3$count = c(1:length(anno3))
    expect_error({s5 = score.hypotheses(anno3)})
    anno$cov2 = NULL
    anno$count = round(runif(length(anno))*5)
    s6 = score.hypotheses(anno, subsample = 0.7)
    anno$cov = NULL
    anno$p = NULL
    m1 = score.hypotheses(anno, return.model = TRUE)
    s8 = score.hypotheses(anno[1:1000], model = m1)
    anno$fac = c(1:2)
    anno$fac = as.factor(anno$fac)
    anno$cov = as.numeric(c(1:length(anno)))
    anno4 = anno[1:1000]
    s9 = score.hypotheses(anno4, verbose = TRUE,  covariates = c('cov', 'fac'))
})




## Covariate

test_that('Covariate', {

    ## initialize = function(..., name = NULL, data = NULL, pad = NULL, type = NULL, signature = NULL, field = NULL, na.rm = NULL, grep = NULL)
    ## default
    foobar = Covariate$new(data = replication_timing, type='numeric')
    foobar2 = Covariate$new(data=replication_timing[50:100], type='numeric')
    foobar3 = Covariate$new(data=list(replication_timing,replication_timing), type = c('numeric','numeric'), name = c('1','1'), na.rm = F)
    ## merge
    #### check it runs
    ##expect_error(foobar$merge(foobar2), NA)
    ## chr
    ## expect_false(foobar$chr())
    ## seqlevels
    expect_equal(length(foobar$seqlevels()[[1]]), 25)
    ## subset
    foobar3$subset(1)
    expect_equal(length(foobar3$data), 1)
    ## toList
    ##expect_equal(length(foobar$toList()[[1]]$track), 2385966)
    ##expect_equal(foobar$toList()[[1]]$type, 'numeric')
    ##expect_equal(foobar$toList()[[1]]$signature, NA)
    ##expect_equal(foobar$toList()[[1]]$pad, NA)
    ##expect_equal(foobar$toList()[[1]]$na.rm, NA)
    ##expect_equal(foobar$toList()[[1]]$field, NA)
    ##expect_equal(foobar$toList()[[1]]$grep, NA)
    ## print
    ##expect_equal(foobar$print()[[1]], NULL)
    foobar4 = Covariate$new(data = list(replication_timing, replication_timing), type = c('numeric','numeric'))
    expect_equal(any(foobar4$chr()), FALSE)
    expect_equal(length(foobar4$toList()), 2)
    expect_equal(foobar4$toList()[[2]]$type, 'numeric')
    expect_output(foobar4$print())
    foobar5 = c(foobar4,foobar4)
    expect_equal(length(foobar5$data), 4)
    foobar6 = foobar5[c(1,3)]
    expect_equal(length(foobar6$data), 2)
    ##Testing print empty cov_arr
    empty = Covariate$new()
    ##Testing active binding errors in cov_arr    
    ##Names
    expect_error({foobar3$names = 1})
    expect_error({foobar3$names = c('1','2')})
    ##Type    
    foobar3 = c(foobar3, foobar3, foobar3, foobar3)
    expect_error({foobar3$type = 1})
    foobar3$type = c('sequence')
    expect_error({foobar3$type = 'hello'})
    foobar3$type = 'interval'
    expect_error({foobar3$type = rep('interval', 10)})
    ##pad
    expect_error({foobar3$pad = '1'})
    foobar3$pad = 1
    expect_equal(foobar3$pad[3], 1)
    expect_error({foobar3$pad = c(1,1,1,1,1,1,1,1,1,1,1,1)})
    ##na.rm
    expect_error({foobar3$na.rm = '1'})
    expect_error({foobar3$na.rm = c('1','2')})
    expect_error({foobar3$na.rm = c(T,T,T,T,T,T,T,T,T,T,T,T)})
    ##grep
    expect_error({foobar3$grep = 1})
    expect_error({foobar3$grep = c('1','2','1','2','3','4','5','6')})
    ##Sequence
    #expect_error({foobar3$signature = 1})
    ##Field
    expect_error({foobar3$field = 1})
    foobar3$field = '1'
    expect_error({foobar3$field = c('1','1','1','1','1','1','1','1','1')})
    expect_equal(foobar3$field[1], '1')
    ##Covariate concatentation error on non cov_arr inputs
    expect_error({foobar3$signature = c(1,2,3,4,5,6,7,8)})
    ##covariate merging
    foobar7 = foobar5$merge(foobar5)
    expect_equal(length(foobar7$data), 8)
    ##chr function
    foobar8 = Covariate$new()
    expect_equal(foobar8$chr(), NULL)
    foobar9 = expect_error(Covariate$new(data = list(replication_timing, 'hello'), name = c('r1','r2')))
    ##seqlevels function
    expect_equal(foobar8$seqlevels(), NULL)    
    ##empty output print
    expect_equal(foobar8$print(), NULL)
    ##c with item that is not Covariate
    expect_error({foobar10 = c(foobar2, 1)})       
})




## FishHook

## initialize = function(hypotheses = NULL, out.path = NULL, eligible = NULL, ... ,events = NULL, covariates = NULL,
##     use_local_mut_density = FALSE, local_mut_density_bin = 1e6, genome = 'BSgenome.Hsapiens.UCSC.hg19::Hsapiens',
##     mc.cores = 1, na.rm = TRUE, pad = 0, verbose = TRUE, max.slice = 1e3, ff.chunk = 1e6, max.chunk = 1e11, ptidcol = NULL,
##     maxpatientpergene = Inf, weightEvents = FALSE, nb = TRUE)

test_that('FishHook', {

    ##FishHook initialization warnings:
  ##Mismatched seqlevels
  set.seed(42)
  events = sample(events, 1e5)
    r1 = replication_timing[1]
    r2 = replication_timing[1]
    r2 = gr2dt(r2)
    r2$seqnames = 'chr1'
    r2 = dt2gr(r2)
    c1 = Covariate$new(data = r1, type = 'interval', name = 'c1')
    c2 = Covariate$new(data = r2, type = 'interval', name = 'c2')
    t1 = hypotheses[1]
    t2 = t1
    t2 = gr2dt(t2)
    t2$seqnames = 'chr1'
    t2 = dt2gr(t2)
    expect_warning({fish = FishHook$new(hypotheses = t1, events = t1,  covariates = c(c1,c2))})
    expect_warning({fish = FishHook$new(hypotheses = t1, events = t1, covariates = c(c2))})
    #expect_error({fish = FishHook$new(hypotheses = t1, events = t1, eligible = t2)})
    #expect_error({fish = FishHook$new(hypotheses = t1, events = t2)})
    ##Improper Cov class
    expect_error({fish = FishHook$new(hypotheses = t1, events = t1, covariates = t1)})
    t1 = hypotheses[1]
    expect_error({fish = FishHook$new(hypotheses = 'hello', events = t1)})
    expect_error({fish = FishHook$new(hypotheses = t1, events = 'hello')})
    expect_error({fish = FishHook$new(hypotheses = t1, events = t1, eligible = 'hello')})
    ##local mut_density with covariates present
    r1 = replication_timing[1]
    c1 = Covariate$new(data = r1, type = 'interval', name = 'c1')
    fish = FishHook$new(hypotheses = hypotheses, events = events, eligible = eligible, covariates = c(c1), use_local_mut_density = T)
    ##Trying to score after we have already scored
    fish = FishHook$new(hypotheses = hypotheses, events = events, eligible = eligible)
    fish$score()
    fish$aggregate(rolling = 4)
    expect_error({fish$aggregated = fish$anno})
    fish$score()    
    ##Printing functions, all regions eligible and no covs
    t1 = hypotheses[1]
    fish = FishHook$new(hypotheses = t1, events = t1)
    expect_output(fish$print())
    ##QQ_plot
    fish = FishHook$new(hypotheses = hypotheses, events = events, eligible = eligible)
    fish$score()
    p = fish$qqp(columns = c('gene_name', 'strand'))
    expect_equal(grepl('gene_name', p$x[[3]][[1]]$text[1]) && grepl('strand', p$x[[3]][[1]]$text[1]), T)   
    ## default
    fish1 = FishHook$new(hypotheses = hypotheses, events = events, eligible =eligible, use_local_mut_density = T)
    expect_equal(length(fish1$toList(fish1$covariates)), 1)
    expect_output(fish1$print())
    expect_true(fish1$nb)
    expect_equal(fish1$state, 'Annotated')
    expect_false(fish1$weightEvents)
    anno = fish1$state
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
    ## score
    ##
    ## Error in score.hypotheses(targ, covariates = names(values(private$panno)),  : 
    ## Error: "score.hypotheses" input malformed --> count does not vary!
    ## qq_plot
    expect_error(fish1$qqp())
    ## clear
    ## with eligible
    ## foobar = FishHook$new(hypotheses = hypotheses, events = events, eligible = eligible)
    fish2 = FishHook$new(hypotheses = hypotheses, events = events,
                         covariates = Covariate$new(data = replication_timing,
                                                  name = 'rept', type = 'numeric',
                                                  field = NA, pad = 10))

    fish3 = fish2$clone()
    replication_timing$bla = runif(length(replication_timing))
    fish3$merge(Covariate$new(data = replication_timing[1:1000], field = 'bla', name = 'c.new'))
    fish3$merge(fish2)
    ## test set analysis
    expect_equal(dim(fish3)[2],3)
    set.seed(42)
    fish2$sets = split(sample(length(fish2), 100), 1:20)
    expect_equal(nrow(fish2$setres), 20)
    res = fishHook::score(fish2[1:200,], fish2[1:200, ])
    ## test active bindings
    ## data
    expect_equal(length(fish2$covariates), 1)
    ## eligible
    expect_equal(fish2$eligible, NULL)
    anno2 = fish2$state
    print('created anno2')
    ## anno
    expect_equal(max(fish2$data$count), 500)
    ## hypotheses 
    expect_equal(length(fish2$hypotheses), 19688)
    ##Scoring
    fish2$score()
    #expect_error(fish2$score())
    expect_equal(ncol(fish2$res), 15)
    ##Clearing
    fish2$clear('Annotated')
    expect_equal(fish2$state, 'Annotated')
    expect_error(fish2$res)
    fish2$score()
    plot = fish2$qqp()
    expect_equal(length(plot[[1]]), 7)           
    ##testing rolling aggregation
    start = c(1,1001,2001,3001,4001,5001)
    end = c(1000,2000,3000,4000,5000,6000)
    chr = 1
    test = dt2gr(as.data.table(cbind(start,end,chr)))
    mcols(test) = NULL
    fishAgg = FishHook$new(hypotheses = test, events = test)
    fishAgg$aggregate(rolling = 3)
    fishAgg$clear('')
    expect_error(fishAgg$aggregate(rolling = 0))
    expect_error(fishAgg$aggregate(rolling = 4.5))
    expect_equal(length(fishAgg$aggregated), 4)
    expect_equal(any(!width(fishAgg$aggregated) == rep(3000,4)), FALSE)
    ##Interval Covariates
    int = Covariate$new(data = replication_timing[1], name = 'int', type = 'interval', pad = 10)
    fish3 = FishHook$new(hypotheses = hypotheses, events = events[1], covariates = int)
    fish3$score()
    ##Mismatched warning
    mis_eve = gr2dt(events[1])
    mis_eve$seqnames = 'f'
    mis_eve = dt2gr(mis_eve)
    #expect_error({fish4 = FishHook$new(hypotheses = hypotheses, events = mis_eve, eligible = hypotheses[1])})
    #expect_error({fish4 = FishHook$new(hypotheses = hypotheses, events = events[1],
    #                                  eligible = mis_eve,
    #                                  covariates = Covariate$new(csv = c(mis_eve,events[1]) ,
    #                                              type = c('interval'), name = c('mis','eve')))})
    #expect_error({fish5 = FishHook$new(hypotheses = hypotheses, events = events[1], eligible = mis_eve)})
    ##Testing FishHook subsetting
    fish2$eligible = eligible[1:100]
    expect_equal(fish2$state, 'Annotated')
    fish2$covariates = c(int,int,int,int)
    expect_equal(fish2$state, 'Annotated')
    fish5 = fish2
    fish6 = fish5[c(1:10000), c(1,2)]
    expect_equal(length(fish6$hypotheses), 10000)
    expect_equal(length(fish6$covariates$data), 2)
    ##Testing fishHook active Bindings
    #data
    expect_error({fish2$data = 'Hello World'})
    #eligible
    expect_error({fish2$eligible = 'Hello World'})
    #hypotheses
    expect_error({fish2$hypotheses = 42})
    ##Need to fix issue with loading bed files
    #fish2$hypotheses = '~/git/fishHook/data/hypotheses.bed'
      ##Events
    expect_error({fish2$events = 'Hello World'})
    expect_equal(length(fish2$events), 1e5)
    ##Out path
    expect_error({fish2$out.path = 12})
#    fish2$out.path = '/home/travis/build/mskilab/fishHook/data/out'
#    expect_error({fish2$out.path = '~/git/fishHook/data/out'})
 #   expect_equal(fish2$out.path, '/home/travis/build/mskilab/fishHook/data/out')
    ##anno
    expect_error({fish2$anno = 'hello world'})
    expect_warning({fish2$data = fish2$data[, 1:4]})
    ##scores
    expect_error({fish2$scores = 'Hello World'})
    x = as.data.table('a')
    ##Model
    expect_equal(fish2$model, NULL)
    expect_error({fish2$model = 'hello'})
    ##mc.cores
    expect_error({fish2$mc.cores = 'hello'})
    fish2$mc.cores = 2
    expect_equal(fish2$mc.cores, 2)
    ##na.rm
    expect_error({fish2$na.rm = 'hello'})
    fish2$na.rm = T
    expect_equal(fish2$na.rm, T)
    ##pad
    expect_error({fish2$pad = 'hello'})
    fish2$pad = 2
    expect_equal(fish2$pad, 2)
    ##verbose
    expect_error({fish2$verbose = 'hello'})
    fish2$verbose = T
    expect_equal(fish2$verbose, T)
    ##max.slice
    expect_error({fish2$max.slice = 'hello'})
    fish2$max.slice = 2
    expect_equal(fish2$max.slice, 2)
    ##ff.chunk
    expect_error({fish2$ff.chunk = 'hello'})
    fish2$ff.chunk = 2
    expect_equal(fish2$ff.chunk, 2)
    ##max.chunk
    expect_error({fish2$max.chunk = 'hello'})
    fish2$max.chunk = 2
    expect_equal(fish2$max.chunk, 2)
    ##ptidcol
    expect_error({fish2$idcol = 1})
    expect_error({fish2$idcol = "hello"})
    fish2$max.chunk = 1e11
    fish2$max.slice = 100000
    fish2$idcol = 'id'
    expect_equal(fish2$idcol, "id")
    ##maxpatientpergene
    expect_error({fish2$idcap = 'hello'})
    fish2$idcap = 2
    expect_equal(fish2$idcap, 2)
    ##weightEvents
    expect_error({fish2$weightEvents = 'hello'})
    fish2$ff.chunk = 1e7
    fish2 = fish2[1:1000, ]
    fish2$weightEvents = T
    expect_equal(fish2$weightEvents, T)
    ##nb
    expect_error({fish2$nb = 'hello'})
    fish2$nb = T
    expect_equal(fish2$nb, T)
    ##all
    expect_error({fish2$all = 'hello'})
    ##state
    expect_error({fish2$state = 'hello'})
    ##agg
    expect_error({fish2$aggregated = 'hello'})
    grl = GRangesList(events[1], events[1], events[1])
    expect_error({fish2$aggregated = grl})
})


### Annotate

## anno= Annotated$new(hypotheses=hypotheses, events=events, covariates=replication_timing)
## SLOW



test_that('qqp', {
    pvals = c(0.0001, 0.0001, 0.032, 0.005, 0.9, 0.15)
    ## check default params
    foo = qqp(pvals, plotly = FALSE)
    expect_match(names(foo)[1], 'rect')
    expect_match(names(foo)[2], 'text')
    expect_equal(round(foo$rect$w, 2), round(0.8597488, 2))  ## rounding; otherwise, 'double' types not exactly the same
    expect_equal(round(foo$rect$h, 2), round(0.4891113, 2))
    expect_equal(round(foo$rect$left, 2), round(3.820251, 2))
    expect_equal(round(foo$rect$top, 2), round(0.3091113, 2))
    expect_equal(round(foo$text$x, 2), round(4.073376, 2))
    expect_equal(round(foo$text$y, 2), round(0.008372093, 2))
    ## exp: Error: length of exp must be = length(obs
    expect_error(qqp(pvals, exp=(c(1, 2, 3))))
    ## not sure how to test 'lwd'
    ## or these other args
    foobar = qqp(pvals, plotly=TRUE)
    expect_match(names(foobar)[1], 'x')
    expect_match(names(foobar)[2], 'width')
    expect_match(names(foobar)[3], 'height')
    dist = rnorm(1000, mean = 10)/1000
    p1 = qqp(dist)
    p2 = qqp(dist, exp = rnorm(1000, mean = 8)/1000)
    p3 = qqp(dist, exp = rnorm(1000, mean = 8)/1000, highlight = rep(c(T,F),500))
    p4 = qqp(dist, exp = rnorm(1000, mean = 8)/1000, highlight = rep(c(T,F),500), max.x = 10, max.y = 10)
    p4 = qqp(dist, exp = rnorm(1000, mean = 8)/1000, highlight = rep(c(T,F),500), max.x = 10, max.y = 10, subsample = 100)
    p5 = qqp(plotly = T, obs = dist, exp = rnorm(1000, mean = 8)/1000, highlight = rep(c(T,F),500))
    p6 = qqp(plotly = T, obs = dist, exp = rnorm(1000, mean = 8)/1000, highlight = rep(c(T,F),500), gradient = list('arb' = c(1:1000)))
    p7 = qqp(plotly = T, obs = rnorm(10000,mean = 10)/1000, exp = rnorm(10000, mean = 8)/10000, highlight = rep(c(T,F),5000), gradient = list('arb' = c(1:10000)))
    p8 = qqp(plotly = T, obs = rnorm(10000,mean = 10)/1000, exp = rnorm(10000, mean = 8)/10000, highlight = rep(c(T,F),5000))
    p9 = qqp(plotly = T, obs = rnorm(10000,mean = 10)/1000, exp = rnorm(10000, mean = 8)/10000, highlight = 'hello')
    expect_error({p10 = qqp(plotly = T, obs = rnorm(10000,mean = 10)/1000, exp = rnorm(10000, mean = 8)/10000, highlight = rep(c(T,F), 1000000))})
    p11 = qqp(obs = rnorm(10000,mean = 10)/1000, exp = rnorm(10000, mean = 8)/10000, label = c(1:1000))
    p12 = qqp(obs = rnorm(10020,mean = 10)/1000, exp = rnorm(10020, mean = 8)/10020)
    expect_error({p13 = qqp(obs = rnorm(10000,mean = 10)/1000, exp = rnorm(10000, mean = 8)/10000, highlight = rep(c(T,F), 100000000))})
    p14 = qqp(obs = rnorm(10000,mean = 10)/1000, exp = rnorm(10000, mean = 8)/10000, highlight = c(1:5000))
})



