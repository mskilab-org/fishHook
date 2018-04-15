
library(fishHook)
library(testthat)
##ZG Testing Paths


## load('~/git/fishHook/data/events.rda')
## load('~/git/fishHook/data/targets.rda')
## load('~/git/fishHook/data/replication_timing_cov.rda')
## eligible = load('~/git/fishHook/data/eligible.rda')
## anno = readRDS('~/git/fishHook/data/anno.rds')





Sys.setenv(DEFAULT_BSGENOME = 'BSgenome.Hsapiens.UCSC.hg19::Hsapiens')
# Sample Events


## Events
## events = readRDS(system.file("data", "events.rds", package = 'fishHook'))

## Sample Targets
targets = readRDS(system.file("data", "targets.rds", package = 'fishHook'))

## targets BED
targetsbed = system.file("data", "targets.bed", package = 'fishHook')

## Sample Covariate
replication_timing = readRDS(system.file("data", "covariate.rds", package = 'fishHook'))

## Same Eligible Subset
## eligible = readRDS(system.file("data", "eligible.rds", package = 'fishHook'))

## indexed pathways
indexed_pathways = readRDS(system.file("data", "indexed_pathways.rds", package = 'fishHook'))

## segs 
segs = readRDS(system.file("data", "jabba_segs_11517.rds", package = 'fishHook'))

## eligible
## eligible = readRDS(system.file("data", "eligible.rds", package = 'fishHook'))

# Sample annotate
anno = readRDS(system.file("data", "anno.rds", package = 'fishHook'))



context('unit testing fishhook operations')



## annotate.targets = function(targets, covered = NULL, events = NULL,  mc.cores = 1, na.rm = TRUE, pad = 0, verbose = TRUE, max.slice = 1e3, 
##    ff.chunk = 1e6, max.chunk = 1e11, out.path = NULL, covariates = list(), maxpatientpergene = Inf, ptidcol = NULL, weightEvents = FALSE, ...)

test_that('annotate.targets', {

    ## default args
    t1 = targets[1]
    e1 = events[1]
    el1 = eligible[1]
    r1 = replication_timing[1]
    expect_error({annotate.targets(targets = targets[1], events = events[1], maxpatientpergene = 'hello')})
    c1 = Cov_Arr$new(cvs = list(system.file("data", "covariate.rds", package = 'fishHook')), type = 'numeric', name = 'c1', pad = 10)
    c2 = Cov_Arr$new(cvs = list(system.file("data", "covariate.rds", package = 'fishHook')), type = 'interval', name = 'c2', pad = 10, na.rm = TRUE)
    c3 = Cov_Arr$new(cvs = list(RleList()), type = 'interval', name = 'c3', pad = 10, na.rm = TRUE)
    anno = annotate.targets(targets = t1, events = e1, covariates = c1$toList())
    anno2 = annotate.targets(targets = t1, events = e1, covariates = c2$toList())
    #anno3 = annotate.targets(targets = t1, events = e1, covariates = c3$toList())
    anno4 = annotate.targets(targets = e1, events = e1, eligible = e1)
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
    annotate.targets(system.file("data", "targets.rds", package = 'fishHook'))
    expect_true(is(annotate.targets(targets), 'GRanges'))
    ##Interval Covariates
    int = Cov_Arr$new(cvs = replication_timing[1], name = 'int', type = 'interval', pad = 10)
    fish3 = FishHook$new(targets = targets, events = events[1], covariates = int)
    fish3$annotate()
    fish3$score()
    ## 
    ## else if (grepl('(\\.bed$)', targets[1])){
    expect_equal(length(annotate.targets(targetsbed)), 19688)
    ## if (!is.null(out.path)){
    ## expect_equal(length(annotate.targets(targetsbed, out.path='/home/travis/build/mskilab/fishHook/data/foobar.rds')))
    

})



## aggregate.targets

c(system.file("data", "anno.rds", package = 'fishHook'), system.file("data", "anno.rds", package = 'fishHook'))

test_that('aggregate.targets', {

    expect_error(aggregate.targets(system.file("data", "targets.rds", package = 'fishHook'), rolling = 1))
    expect_error(aggregate.targets(system.file("data", "targets.rds", package = 'fishHook'), rolling = 1))
    agg = aggregate.targets(c(system.file("data", "anno.rds", package = 'fishHook'), system.file("data", "anno.rds", package = 'fishHook')), rolling = 1, na.rm = T)
    expect_error({agg = aggregate.targets(anno, rolling = -1)})
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
    expect_error(aggregate.targets(system.file("data", "targets.rds", package = 'fishHook')))  ## Coverage missing for input targets
    ##  if (is.null(by) & is.character(targets)){ (continued)
    expect_error(aggregate.targets(system.file("data", "annotated_cov.rds", package = 'fishHook')))
    ##testing rolling aggregation
    start = c(1,1001,2001,3001,4001,5001)
    end = c(1000,2000,3000,4000,5000,6000)
    chr = 1
    test = dt2gr(as.data.table(cbind(start,end,chr)))
    mcols(test) = NULL
    fishAgg = FishHook$new(targets = test, events = test)
    fishAgg$aggregate(rolling = 3)
    expect_error(fishAgg$aggregate(rolling = 0))
    expect_error(fishAgg$aggregate(rolling = 4.5))
    expect_equal(length(fishAgg$aggregated), 4)
    expect_equal(any(!width(fishAgg$aggregated) == rep(3000,4)), FALSE)

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

    ## MISC testing
    expect_error({s1 = score.targets(anno, covariates = c('hello'))})
    anno$cov2 = c(1:length(anno))
    anno$cov2 = as.character(anno$cov2)
    expect_error({s2 = score.targets(anno)})
    anno3 = anno
    anno3$count = NULL
    anno3$cov2 = NULL
    expect_error({s3 = score.targets(anno3)})
    anno3$count = 1
    expect_error({s4 = score.targets(anno3)})
    anno3$dummy = NA
    anno3$dummy = as.numeric(anno3$dummy)
    anno3$count = c(1:length(anno3))
    expect_error({s5 = score.targets(anno3)})
    anno$cov2 = NULL
    s6 = score.targets(anno, subsample = 0.7)
    anno$cov = NULL
    anno$p = NULL
    m1 = score.targets(anno, return.model = TRUE)
    s8 = score.targets(anno[1:1000], model = m1)
    anno$fac = c(1:2)
    anno$fac = as.factor(anno$fac)
    anno$cov = as.numeric(c(1:length(anno)))
    anno4 = anno[1:1000]
    s9 = score.targets(anno4, verbose = TRUE,  covariates = c('cov', 'fac'))
})




## Cov_Arr

test_that('Cov_Arr', {

    ## initialize = function(..., name = NULL, cvs = NULL, pad = NULL, type = NULL, signature = NULL, field = NULL, na.rm = NULL, grep = NULL)
    ## default
    foobar = Cov_Arr$new(cvs = replication_timing, type='numeric')
    foobar2 = Cov_Arr$new(cvs=replication_timing[50:100], type='numeric')
    foobar3 = Cov_Arr$new(cvs=list(replication_timing,replication_timing), type = c('numeric','numeric'), name = c('1','1'), is.na = F)
    ## merge
    #### check it runs
    ##expect_error(foobar$merge(foobar2), NA)
    ## chr
    ## expect_false(foobar$chr())
    ## seqlevels
    expect_equal(length(foobar$seqlevels()[[1]]), 25)
    ## subset
    foobar3$subset(1)
    expect_equal(length(foobar3$cvs), 1)
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
    foobar4 = Cov_Arr$new(cvs = list(replication_timing, replication_timing), type = c('numeric','numeric'))
    expect_equal(any(foobar4$chr()), FALSE)
    expect_equal(length(foobar4$toList()), 2)
    expect_equal(foobar4$toList()[[2]]$type, 'numeric')
    expect_output(foobar4$print())
    foobar5 = c(foobar4,foobar4)
    expect_equal(length(foobar5$cvs), 4)
    foobar6 = foobar5[c(1,3)]
    expect_equal(length(foobar6$cvs), 2)
    ##Testing print empty cov_arr
    empty = Cov_Arr$new()
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
    ##Cov_Arr concatentation error on non cov_arr inputs
    expect_error({foobar3$signature = c(1,2,3,4,5,6,7,8)})
    ##covariate merging
    foobar7 = foobar5$merge(foobar5)
    expect_equal(length(foobar7$cvs), 8)
    ##chr function
    foobar8 = Cov_Arr$new()
    expect_equal(foobar8$chr(), NULL)
    foobar9 = Cov_Arr$new(cvs = list(replication_timing, 'hello'), name = c('r1','r2'))
    expect_equal(foobar9$chr()[2], NA)
    expect_equal(foobar9$chr()[1], F)
    ##seqlevels function
    expect_equal(foobar8$seqlevels(), NULL)
    expect_equal(foobar9$seqlevels()[[2]], NA)
    expect_equal(length(foobar9$seqlevels()[[1]]), 25)
    ##empty output print
    expect_output(foobar8$print())
    expect_equal(foobar8$print(), NULL)
    ##c with item that is not Cov_Arr
    expect_error({foobar10 = c(foobar2, 1)})
    


    
})




## FishHook

## initialize = function(targets = NULL, out.path = NULL, eligible = NULL, ... ,events = NULL, covariates = NULL,
##     use_local_mut_density = FALSE, local_mut_density_bin = 1e6, genome = 'BSgenome.Hsapiens.UCSC.hg19::Hsapiens',
##     mc.cores = 1, na.rm = TRUE, pad = 0, verbose = TRUE, max.slice = 1e3, ff.chunk = 1e6, max.chunk = 1e11, ptidcol = NULL,
##     maxpatientpergene = Inf, weightEvents = FALSE, nb = TRUE)

test_that('FishHook', {

    ##FishHook initialization warnings:
    ##Mismatched seqlevels
    r1 = replication_timing[1]
    r2 = replication_timing[1]
    r2 = gr2dt(r2)
    r2$seqnames = 'chr1'
    r2 = dt2gr(r2)
    c1 = Cov_Arr$new(cvs = r1, type = 'interval', name = 'c1')
    c2 = Cov_Arr$new(cvs = r2, type = 'interval', name = 'c2')
    t1 = targets[1]
    t2 = t1
    t2 = gr2dt(t2)
    t2$seqnames = 'chr1'
    t2 = dt2gr(t2)
    expect_warning({fish = FishHook$new(targets = t1, events = t1,  covariates = c(c1,c2))})
    expect_warning({fish = FishHook$new(targets = t1, events = t1, covariates = c(c2))})
    #expect_error({fish = FishHook$new(targets = t1, events = t1, eligible = t2)})
    #expect_error({fish = FishHook$new(targets = t1, events = t2)})
    ##Improper Cov class
    expect_error({fish = FishHook$new(targets = t1, events = t1, covariates = t1)})
    t1 = targets[1]
    expect_error({fish = FishHook$new(targets = 'hello', events = t1)})
    expect_error({fish = FishHook$new(targets = t1, events = 'hello')})
    expect_error({fish = FishHook$new(targets = t1, events = t1, eligible = 'hello')})
    ##local mut_density with covariates present
    r1 = replication_timing[1]
    c1 = Cov_Arr$new(cvs = r1, type = 'interval', name = 'c1')
    fish = FishHook$new(targets = targets, events = events, eligible = eligible, covariates = c(c1), use_local_mut_density = T)
    ##Trying to score after we have already scored
    fish = FishHook$new(targets = targets, events = events, eligible = eligible)
    fish$score()
    expect_error(fish$annotate())
    fish$clear('Annotated')
    fish$aggregate(rolling = 4)
    fish$aggregated = fish$anno
    fish$score()    
    ##Printing functions, all regions eligible and no covs
    t1 = targets[1]
    fish = FishHook$new(targets = t1, events = t1)
    expect_output(fish$print())
    ##QQ_plot
    fish = FishHook$new(targets = targets, events = events, eligible = eligible)
    fish$score()
    p = fish$qq_plot(columns = c('gene_name', 'strand'))
    expect_equal(grepl('gene_name', p$x[[3]][[1]]$text[1]) && grepl('strand', p$x[[3]][[1]]$text[1]), T)   
    ## default
    fish1 = FishHook$new(targets = targets, events = events, eligible =eligible, use_local_mut_density = T)
    expect_equal(length(fish1$toList(fish1$cvs)), 1)
    expect_output(fish1$print())
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
    fish2 = FishHook$new(targets = targets, events = events,
                         covariates = Cov_Arr$new(cvs = replication_timing,
                                                  name = 'rept', type = 'numeric',
                                                  field = NA, pad = 10))
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
    ##Scoring
    fish2$score()
    #expect_error(fish2$score())
    expect_equal(ncol(fish2$scores), 17)
    ##Clearing
    fish2$clear('Annotated')
    expect_equal(fish2$state, 'Annotated')
    expect_equal(fish2$scores, NULL)
    fish2$score()
    plot = fish2$qq_plot()
    expect_equal(length(plot[[1]]), 7)           
    ##testing rolling aggregation
    start = c(1,1001,2001,3001,4001,5001)
    end = c(1000,2000,3000,4000,5000,6000)
    chr = 1
    test = dt2gr(as.data.table(cbind(start,end,chr)))
    mcols(test) = NULL
    fishAgg = FishHook$new(targets = test, events = test)
    fishAgg$aggregate(rolling = 3)
    fishAgg$clear('')
    expect_error(fishAgg$aggregate(rolling = 0))
    expect_error(fishAgg$aggregate(rolling = 4.5))
    expect_equal(length(fishAgg$aggregated), 4)
    expect_equal(any(!width(fishAgg$aggregated) == rep(3000,4)), FALSE)
    ##Interval Covariates
    int = Cov_Arr$new(cvs = replication_timing[1], name = 'int', type = 'interval', pad = 10)
    fish3 = FishHook$new(targets = targets, events = events[1], covariates = int)
    fish3$annotate()
    fish3$score()
    ##Mismatched warning
    mis_eve = gr2dt(events[1])
    mis_eve$seqnames = 'f'
    mis_eve = dt2gr(mis_eve)
    #expect_error({fish4 = FishHook$new(targets = targets, events = mis_eve, eligible = targets[1])})
    #expect_error({fish4 = FishHook$new(targets = targets, events = events[1],
    #                                  eligible = mis_eve,
    #                                  covariates = Cov_Arr$new(csv = c(mis_eve,events[1]) ,
    #                                              type = c('interval'), name = c('mis','eve')))})
    #expect_error({fish5 = FishHook$new(targets = targets, events = events[1], eligible = mis_eve)})
    ##Testing FishHook subsetting
    fish2$eligible = eligible
    expect_equal(fish2$state, 'Initialized')
    fish2$cvs = c(int,int,int,int)
    expect_equal(fish2$state, 'Initialized')
    fish5 = fish2
    fish6 = fish5[c(1:10000), c(1:100000), c(1,2), c(1:100)]
    expect_equal(length(fish6$targets), 10000)
    expect_equal(length(fish6$events), 100000)
    expect_equal(length(fish6$cvs$cvs), 2)
    expect_equal(length(fish6$eligible), 100)
    ##Testing fishHook active Bindings
    #cvs
    expect_error({fish2$cvs = 'Hello World'})
    #eligible
    expect_error({fish2$eligible = 'Hello World'})
    #targets
    expect_error({fish2$targets = 42})
    ##Need to fix issue with loading bed files
    #fish2$targets = '~/git/fishHook/data/targets.bed'
    fish2$targets = system.file("data", "targets.rds", package = 'fishHook')
    empty_ranges = GRanges()
    expect_error({fish2$targets = empty_ranges})
    ##Events
    expect_error({fish2$events = 'Hello World'})
    expect_equal(length(fish2$events) , 1985704)
    ##Out path
    expect_error({fish2$out.path = 12})
    fish2$out.path = '/home/travis/build/mskilab/fishHook/data/out'
#    expect_error({fish2$out.path = '~/git/fishHook/data/out'})
    expect_equal(fish2$out.path, '/home/travis/build/mskilab/fishHook/data/out')
    ##anno
    expect_error({fish2$anno = 'hello world'})
    expect_warning({fish2$anno = eligible})
    ##scores
    expect_error({fish2$scores = 'Hello World'})
    x = as.data.table('a')
    expect_warning({fish2$scores = x})
    ##Model
    expect_equal(fish2$model, NULL)
    expect_warning({fish2$model = 'hello'})
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
    expect_error({fish2$ptidcol = 1})
    fish2$ptidcol = "hello"
    expect_equal(fish2$ptidcol, "hello")
    ##maxpatientpergene
    expect_error({fish2$maxpatientpergene = 'hello'})
    fish2$maxpatientpergene = 2
    expect_equal(fish2$maxpatientpergene, 2)
    ##weightEvents
    expect_error({fish2$weightEvents = 'hello'})
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
    expect_warning({fish2$aggregated = grl})    
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
    expect_equal(round(foo$rect$w, 2), round(0.8597488, 2))  ## rounding; otherwise, 'double' types not exactly the same
    expect_equal(round(foo$rect$h, 2), round(0.4891113, 2))
    expect_equal(round(foo$rect$left, 2), round(3.820251, 2))
    expect_equal(round(foo$rect$top, 2), round(0.3091113, 2))
    expect_equal(round(foo$text$x, 2), round(4.073376, 2))
    expect_equal(round(foo$text$y, 2), round(0.008372093, 2))
    ## exp: Error: length of exp must be = length(obs
    expect_error(qq_pval(pvals, exp=(c(1, 2, 3))))
    ## not sure how to test 'lwd'
    ## or these other args
    foobar = qq_pval(pvals, plotly=TRUE)
    expect_match(names(foobar)[1], 'x')
    expect_match(names(foobar)[2], 'width')
    expect_match(names(foobar)[3], 'height')
    dist = rnorm(1000, mean = 10)/1000
    p1 = qq_pval(dist)
    p2 = qq_pval(dist, exp = rnorm(1000, mean = 8)/1000)
    p3 = qq_pval(dist, exp = rnorm(1000, mean = 8)/1000, highlight = rep(c(T,F),500))
    p4 = qq_pval(dist, exp = rnorm(1000, mean = 8)/1000, highlight = rep(c(T,F),500), max.x = 10, max.y = 10)
    p4 = qq_pval(dist, exp = rnorm(1000, mean = 8)/1000, highlight = rep(c(T,F),500), max.x = 10, max.y = 10, subsample = 100)
    p5 = qq_pval(plotly = T, obs = dist, exp = rnorm(1000, mean = 8)/1000, highlight = rep(c(T,F),500))
    p6 = qq_pval(plotly = T, obs = dist, exp = rnorm(1000, mean = 8)/1000, highlight = rep(c(T,F),500), gradient = list('arb' = c(1:1000)))
    p7 = qq_pval(plotly = T, obs = rnorm(10000,mean = 10)/1000, exp = rnorm(10000, mean = 8)/10000, highlight = rep(c(T,F),5000), gradient = list('arb' = c(1:10000)))
    p8 = qq_pval(plotly = T, obs = rnorm(10000,mean = 10)/1000, exp = rnorm(10000, mean = 8)/10000, highlight = rep(c(T,F),5000))
    p9 = qq_pval(plotly = T, obs = rnorm(10000,mean = 10)/1000, exp = rnorm(10000, mean = 8)/10000, highlight = 'hello')
    expect_error({p10 = qq_pval(plotly = T, obs = rnorm(10000,mean = 10)/1000, exp = rnorm(10000, mean = 8)/10000, highlight = rep(c(T,F), 1000000))})
    expect_warning({p11 = qq_pval(obs = rnorm(10000,mean = 10)/1000, exp = rnorm(10000, mean = 8)/10000, label = c(1:1000))})
    p12 = qq_pval(obs = rnorm(10020,mean = 10)/1000, exp = rnorm(10020, mean = 8)/10020)
    expect_error({p13 = qq_pval(obs = rnorm(10000,mean = 10)/1000, exp = rnorm(10000, mean = 8)/10000, highlight = rep(c(T,F), 100000000))})
    p14 = qq_pval(obs = rnorm(10000,mean = 10)/1000, exp = rnorm(10000, mean = 8)/10000, highlight = c(1:5000))
})



