#############################################################################
## Marcin Imielinski
## The Broad Institute of MIT and Harvard / Cancer program.
## marcin@broadinstitute.org

## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.

## You should have received a copy of the GNU Lesser General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################


#' annotate.targets
#'
#' Takes input of GRanges targets, an optional set of "covered" intervals, and an indefinite list of covariates which can be R objects
#' (GRanges, ffTrack, Rle) or file paths to .rds, .bw, .bed files, and an annotated target intervals GRanges with covariates computed
#' for each interval.   These target intervals can be further annotated with mutation counts and plugged into a generalized linear regression
#' (or other) model downstream.
#' 
#' There are three types of covariates: numeric, sequence, interval.  The covariates are computed as follows:
#' numeric covariates: the mean value
#' sequence covarites: fraction of bases satisfying $signature
#' interval covariates: fraction of bases overlapping feature
#' 
#' @param targets path to bed or rds containing genomic target regions with optional target name
#' @param covered  optional path to bed or rds containing  granges object containing "covered" genomic regions
#' @param covered  optional path to bed or rds containing ranges corresponding to events (ie mutations etc)
#' @param out.path  out.path to save variable to
#' @param ... paths to sequence covariates whose output names will be their argument names, and each consists of a list with
#' $track field corresponding to a GRanges, RleList, ffTrack object (or path to rds containing that object), $type which can
#' have one of three values "numeric", "sequence", "interval".
#' Numeric tracks must have $score field if they are GRanges), and can have a $na.rm logical field describing how to treat NA values
#' (set to na.rm argument by default)
#' Sequence covariates must be ffTrack objects (or paths to ffTrack rds) and require an additional variables $signatures, which
#' will be used as input to fftab, and can have optional logical argument $grep to specify inexact matches (see fftab)
#' Interval covariates must be Granges (or paths to GRanges rds) or paths to bed files
#' @param out.path  out.path to save variable to
#' @return GRanges of input targets annotated with covariate statistics (+/- constrained to the subranges in optional argument covered)
#' @author Marcin Imielinski
#' @importFrom ffTrack fftab
#' @export
annotate.targets = function(targets, ## path to bed or rds containing genomic target regions with optional target name
    covered = NULL,
    events = NULL,  
    ...,
    mc.cores = 1,
    na.rm = TRUE,
    pad = 0,
    verbose = TRUE,
    max.slice = 1e3, ## max slice of intervals to evaluate with  gr.val
    ff.chunk = 1e6, ## max chunk to evaluate with fftab
    max.chunk = 1e11, ## gr.findoverlaps parameter
    out.path = NULL,
    covariates = list())
    {
        if (is.character(targets))
            if (grepl('\\.rds$', targets[1]))
                targets = readRDS(targets[1])
            else if (grepl('(\\.bed$)', targets[1]))
                targets = import.ucsc(targets[1])

        if (length(targets)==0)
            stop('Must provide non-empty targets')
        
        if (!is.null(out.path))
            tryCatch(saveRDS(targets, paste(gsub('.rds', '', out.path), '.source.rds', sep = '')), error = function(e) warning(sprintf('Error writing to file %s', out.file)))
        
        COV.TYPES = c('numeric', 'sequence', 'interval')
        COV.CLASSES = c('GRanges', 'RleList', 'ffTrack', 'character')
        ## covariates = list(...)          
        
        cov.types = sapply(covariates, function(x) if (!is.null(x$type)) x$type else NA)
        cov.classes = lapply(covariates, function(x) if (!is.null(x$track)) class(x$track) else NA)
        
        if (is.list(cov.types))
            cov.type = ''
        
        if (is.list(cov.classes))
            cov.classes = ''
        
        cov.types[is.na(cov.types)] = ''
        
        if (!all(cov.types %in% COV.TYPES) & !(all(cov.classes %in% COV.CLASSES)))
            stop(sprintf('Malformed covariate input: each covariate must be a list with fields $tracks and $type, $track must be of class GRanges, ffTrack, Rle, object, or character path to rds object with the latter or .bw, .bed file, $type must have value %s', paste(COV.TYPES, collapse = ',')))


        if (any(ix <- (cov.types == 'sequence')))
            {
                for (cov in covariates[ix])
                    {
                        if (is.character(cov$track))                      
                            cov$track = tryCatch(readRDS(cov$track), error = function(e) 'not ffTrack')                        
                        
                        if (class(cov$track) != 'ffTrack')
                            stop('sequence tracks must have ffTrack object as $track field or $track must be a path to an ffTrack object rds file')
                    }
            }
        
        if (any(ix <- (cov.classes == 'ffTrack' & cov.types == 'sequence')))
            {
                if (!all(sapply(covariates, function(x) !is.null(x$signature))))
                    stop('sequence tracks must be ffTracks and have a $signature field specified (see fftab in ffTrack)')                    
            }

        if (verbose)
            cat('Overlapping with covered intervals\n')
        if (!is.null(covered))
            ov = gr.findoverlaps(targets, covered, verbose = verbose, max.chunk = max.chunk, mc.cores = mc.cores)
        else
            {
                ov = targets[, c()]
                ov$query.id = ov$subject.id = 1:length(targets)
            }
        
       if (verbose)
           cat('Finished overlapping with covered intervals\n')

        if (length(ov)>0)
            {
                if (!is.null(events))
                    {
                        if (is(events, 'GRanges'))
                            {
                                ev = gr.fix(events[gr.in(events, ov)])                                

                                ## weighing each event by width means that each event will get one
                                ## total count, and if an event is split between two tile windows
                                ## then it will contribute a fraction of event proprotional to the number
                                ## oof bases overlapping
                                counts = coverage(ev, weight = 1/width(ev))
                                oix = which(gr.in(ov, events))
                            }
                        else ## assume it is an Rle of event counts along the genome
                            {
                                counts = events
                                oix = 1:length(ov)
                            }

                        if (verbose)
                            cat('Computing event counts\n')

                        ov$count = 0

                        if (length(oix)>0)
                            ov$count[oix] = fftab(counts, ov[oix], chunksize = ff.chunk, na.rm = TRUE, mc.cores = mc.cores, verbose = verbose)$score

                        if (!is.null(out.path))
                            tryCatch(saveRDS(ov, paste(out.path, '.intermediate.rds', sep = '')), error = function(e) warning(sprintf('Error writing to file %s', out.file)))

                        if (verbose)
                            cat('Finished counting events\n')
                    }
                
                for (nm in names(covariates))
                    {
                        cov = covariates[[nm]]
                        if (verbose)
                            cat('Starting track', nm, '\n')
                        
                        if (cov$type == 'sequence')
                            {
                                if (is.null(cov$grep))
                                    cov$grep = FALSE

                                if (verbose)
                                    cat('Starting fftab for track', nm, '\n')

                                if (!is.list(cov$signature))
                                    cov$signature = list(cov$signature)

                                if (is.null(names(cov$signature)))
                                    if (length(cov$signature)>1)
                                        names(cov$signature) = 1:length(cov$signature)

                                if (!is.null(names(cov$signature)))
                                    names(cov$signature) = paste(nm, names(cov$signature), sep = '.')
                                else
                                    names(cov$signature) = nm

                                if (is.null(cov$pad))
                                    cov$pad = pad
                                
                                val = fftab(cov$track, ov + cov$pad, cov$signature, chunksize = ff.chunk, FUN = mean, na.rm = TRUE, grep = cov$grep, mc.cores = mc.cores)
                                values(ov) = values(val)
                                
                                if (verbose)
                                    cat('Finished fftab for track', nm, '\n')

                                if (!is.null(out.path))
                                    tryCatch(saveRDS(ov, paste(out.path, '.intermediate.rds', sep = '')), error = function(e) warning(sprintf('Error writing to file %s', out.file)))
                            }
                        else if (cov$type == 'numeric')
                            {
                                if (is.character(cov$track))
                                    if (grepl('.rds$', cov$track))
                                        cov$track = readRDS(cov$track)
                                    else ## assume it is a UCSC format
                                        {
                                            require(rtracklayer)
                                            cov$track = import(cov$track)
                                        }

                                if (is.null(cov$pad))
                                    cov$pad = pad

                                if (is(cov$track, 'ffTrack') | is(cov$track, 'RleList'))
                                    {
                                        val = fftab(cov$track, ov + cov$pad, signature = cov$signature, FUN = sum, chunksize = ff.chunk, grep = cov$grep, mc.cores = mc.cores)
                                        values(ov) = values(val)
                                    }
                                else ## then must be GRanges
                                    {
                                        if (is.null(cov$field))
                                            cov$field = 'score'

                                        if (is.null(cov$na.rm))
                                            cov$na.rm = na.rm
                                        new.col = data.frame(val = values(gr.val(ov + cov$pad, cov$track, cov$field, mc.cores = mc.cores, max.slice = max.slice, max.chunk = max.chunk, mean = TRUE, na.rm = cov$na.rm))[, cov$field])
                                        names(new.col) = nm
                                        values(ov) = cbind(values(ov), new.col)                                
                                    }
                                
                                if (!is.null(out.path))
                                    tryCatch(saveRDS(ov, paste(out.path, '.intermediate.rds', sep = '')), error = function(e) warning(sprintf('Error writing to file %s', out.file)))
                            }
                        else if (cov$type == 'interval')
                            {
                                if (is.character(cov$track))
                                    if (grepl('.rds$', cov$track))
                                        cov$track = readRDS(cov$track)
                                    else ## assume it is a UCSC format
                                        {
                                            require(rtracklayer)
                                            cov$track = import(cov$track)
                                        }
                                
                                if (is(cov, 'GRanges'))
                                    stop('Interval tracks must be GRanges')

                                if (is.null(cov$pad))
                                    cov$pad = pad
                                
                                if (is.null(cov$na.rm))
                                    cov$na.rm = na.rm

                                cov$track = reduce(cov$track)
                                
                                new.col = data.frame(val = gr.val(ov + cov$pad, cov$track[, c()], mean = FALSE, weighted = TRUE, mc.cores = mc.cores, max.slice = max.slice, max.chunk = max.chunk, na.rm = TRUE)$value/(width(ov)+2*cov$pad))
                                new.col$val = ifelse(is.na(new.col$val), 0, new.col$val)
                                names(new.col) = nm
                                values(ov) = cbind(values(ov), new.col)

                                if (!is.null(out.path))
                                    tryCatch(saveRDS(ov, paste(out.path, '.intermediate.rds', sep = '')), error = function(e) warning(sprintf('Error writing to file %s', out.file)))
                            }               
                    }

            }
        
        ovdt = gr2dt(ov)
        
        
        cmd = 'list(coverage = sum(width),';
        if (!is.null(events))
            cmd = paste(cmd, 'count = sum(count)', sep = '')
        else
            cmd = paste(cmd, 'count = NA', sep = '')        
        
        cov.nm = setdiff(names(values(ov)), c('coverage', 'count', 'query.id', 'subject.id'))

        if (length(ov)>0)
            {
                if (length(cov.nm)>0)
                    cmd = paste(cmd,  ',', paste(cov.nm, '= mean(', cov.nm, ')', sep = '', collapse = ', '), ')',  sep = '')
                else
                    cmd = paste(cmd, ')',  sep = '')
                
                ovdta =  ovdt[, eval(parse(text = cmd)), keyby = query.id]
                values(targets) = as(as.data.frame(ovdta[list(1:length(targets)), ]), 'DataFrame')
              
            }
        else
            targets$coverage = 0                
        
        targets$query.id = 1:length(targets)                
        
        ix = is.na(targets$coverage)
        if (any(ix))
            targets$coverage[ix] = 0

        if (!is.null(out.path))
            {
                if (file.exists(paste(out.path, '.intermediate.rds', sep = '')))
                    system(paste('rm',  paste(out.path, '.intermediate.rds', sep = '')))
                tryCatch(saveRDS(targets, out.path), error = function(e) warning(sprintf('Error writing to file %s', out.file)))             
            }
        
        return(targets)        
    }



#' aggregate.targets
#'
#' Gathers annotated targets across a vector "by" into meta-intervals returned as a GRangesList, and returns the
#' aggregated statistics for these meta intervals by summing coverage and counts, and performing a weighted average of all other meta data fields
#' (except query.id)
#'
#' If rolling = TRUE, will return a rolling collapse of the sorted input where "rolling" specifies the number of adjacent intervals that are aggregated in a rolling manner.
#' (only makes sense for tiled target sets)
#'
#' If by = NULL and targets is a vector of path names, then aggregation will be done "sample wise" on the files, ie each .rds input will be assumed to comprise the same 
#' intervals in teh same order and aggregation will be computed coverage-weighted mean of covariates, a sum of coverage and counts, and (if present) a Fisher combined
#' of $p values.  Covariates are inferred from the first file in the list.  
#' 
#' @param targets annotated GRanges of targets with fields $coverage, optional field, $count and additional numeric covariates, or path to .rds file of the same
#' @param by  character vector with which to split into meta-territories
#' @param fields by default all meta data fields of targets EXCEPT reserved field names $coverage, $counts, $query.id
#' @param rolling if specified, positive integer specifying how many (genome coordinate) adjacent to aggregate in a rolling fashion
#' @return GRangesList of input targets annotated with new aggregate covariate statistics OR GRanges if rolling is specified
#' @author Marcin Imielinski
#' @import zoo
#' @importFrom data.table setkey := data.table as.data.table
#' @importFrom S4Vectors values values<-
#' @importFrom GenomeInfoDb seqnames
#' @export
aggregate.targets = function(targets, ## path to bed or rds containing genomic target regions with optional target name 
    by = NULL,
    fields = NULL,
    rolling = NULL, ## positive integer with which to performa rolling sum / weighted average WITHIN chromosomes of "rolling" ranges" --> return a granges
    disjoint = TRUE, ## only take disjoint bins of inpu
    na.rm = FALSE, ## only applicable for sample wise aggregation (i.e. if by = NULL)
    FUN = list(), ## only applies (for now) if by = NULL, this is a named list of functions, where each item named "nm" corresponds to an optional function of how to alternatively aggregate field "nm" per samples, for alternative aggregation of coverage and count.  This function is applied at every iteration of loading a new sample and adding to the existing set.   It is normally sum [for coverage and count] and coverage weighted mean [for all other covariates].  Alternative coverage / count aggregation functions should have two arguments (val1, val2) and all other alt covariate aggregation functions should have four arguments (val1, cov1, val2, cov2) where val1 is the accumulating vector and val2 is the new vector of values. 
    verbose = TRUE
    )
    {

  V1 = sn = st = en = keep = count = width = NULL ## NOTE fix
        if (is.null(by) & is.character(targets))
            {
                cat('Applying sample wise merging\n')
            }        
        else if (is.null(by) & is.null(rolling))
            stop('by must be specified and same length as targets or rolling must be non NULL')

        if (is.null(by) & is.character(targets))
            {
                if (!all(ix <- (file.exists(targets)) & grepl('\\.rds$', targets)))
                    {
                        warning(sprintf('%s of the  %s input files for sample wise merging either do not exist or are not .rds files.  Sample wise merging (i.e. when by is null) requires .rds files of equal dimension GRanges (same intervals, same meta data column names)', sum(!ix), length(ix)))
                        if (sum(ix)==0)
                            stop('No files to process')
                        targets = targets[ix]                                                        
                    }
                
                out = readRDS(targets[1])
                gr = out
                if (is.null(out$coverage))
                    stop('Coverage missing for input targets')

                core.fields = c('coverage', 'count', 'p', 'query.id')
                
                cfields = setdiff(names(values(out)), core.fields)

                if (!is.null(fields))
                    cfields = intersect(fields, cfields)

                values(out) = values(out)[, intersect(c(core.fields, cfields), names(values(out)))]

                if (!is.null(out$p))
                    {                        
                        psum = 0
                        psum.df = rep(0, length(out))
                    }

                ### initialize everything to 0
                for (cf in cfields)
                    values(out)[, cf] = 0

                if (!is.null(out$count))
                    out$count = 0

                out$coverage = 0
                out$numcases = length(targets)

                ## for (nm in cfields) ## de normalize coverage
                ##     values(out)[, nm] = as.numeric(values(out)[, nm])*out$coverage

                if (length(targets)>1)
                    for (i in 1:length(targets))
                        {
                            if (verbose)
                                cat('Processing target file', targets[i], '\n')

                            if (i > 1)
                                gr = readRDS(targets[i])
                                                            
                            if (!is.null(out$count))
                                if (!is.null(FUN[['count']]))
                                    out$count = do.call(FUN[['count']], list(out$count, gr$count))
                                else
                                    out$count = as.numeric(out$count) + gr$count
                            
                            for (cf in cfields)
                                {
                                    if (cf %in% names(values(gr)))                                        
                                        val = as.numeric(values(gr)[, cf])
                                    else
                                        {
                                            warning(paste(targets[i], 'missing column', cf))
                                            val = NA
                                        }
                                                                        
                                    if (na.rm)
                                        if (!is.null(FUN[[cf]]))
                                            values(out)[, cf] = ifelse(!is.na(val),
                                                           do.call(FUN[[cf]], list(values(out)[, cf], out$coverage, + val, gr$coverage)),
                                                           values(out)[, cf])
                                        else                                            
                                            values(out)[, cf] = ifelse(!is.na(val),
                                                           (values(out)[, cf]*out$coverage + val*gr$coverage)/(out$coverage + gr$coverage),
                                                           values(out)[, cf])
                                    else
                                        if (!is.null(FUN[[cf]]))
                                            values(out)[, cf] = do.call(FUN[[cf]], list(values(out)[, cf], out$coverage, val, gr$coverage))
                                        else
                                            values(out)[, cf] = (values(out)[, cf]*out$coverage + val*gr$coverage)/(out$coverage + gr$coverage)

                                }
                            
                            if (!is.null(out$p))
                                {
                                    if (!is.null(gr$p))
                                        {
                                            has.val = is.na(gr$p)
                                            psum = ifelse(has.val, psum - 2*log(gr$p), psum)
                                            psum.df = ifelse(has.val, psum.df + 1, psum.df)
                                        }                                            
                                     else
                                        warning(paste(targets[i], 'missing p value column, ignoring for fisher combined computation'))                                                                            
                                }

                            if (is.null(FUN[['coverage']]))
                                out$coverage = as.numeric(out$coverage) + gr$coverage
                            else
                                out$coverage = do.call(FUN[['coverage']], list(as.numeric(out$coverage), gr$coverage))                               
                        }

                if (!is.null(out$p))
                    out$p = pchisq(psum, psum.df, lower.tail = FALSE)                
                return(out)
            }

        if (is.null(fields))
            fields = names(values(targets))
        
        if (any(nnum <- !(sapply(setdiff(fields, 'query.id'), function(x) class(values(targets)[, x])) %in% 'numeric')))
            {
                warning(sprintf('%s meta data fields (%s) fit were found to be non-numeric and not aggregated', sum(nnum), paste(fields[nnum], collapse = ',')))
                fields = fields[!nnum]
            }
        
        cfields = intersect(names(values(targets)), c('coverage', 'count'))
        
        if (is.null(rolling))
            {
                by = as.character(cbind(1:length(targets), by)[,2])
                                                
                if (disjoint)
                    {                
                        tmp.sn = paste(by, seqnames(targets), sep = '_')
                        tmp.dt = data.table(sn = paste(by, seqnames(targets), sep = '_'), st = start(targets), en = end(targets), ix = 1:length(targets))
                        setkey(tmp.dt, sn, st)
                        tmp.dt[, keep := c(TRUE, st[-1]>en[-length(st)]), by = sn]
                        setkey(tmp.dt, ix)
                        targets = targets[tmp.dt$keep, ]
                        if (verbose)
                            cat(sprintf('Removing %s non-disjoint within group intervals, keeping %s\n', prettyNum(sum(!tmp.dt[, keep]), big.mark = ','), prettyNum(sum(tmp.dt[, keep]), big.mark = ',')))
                        by = by[tmp.dt$keep]
                    }
                
                if (verbose)
                    cat('Splitting into GRangesList\n')
                
                out = split(targets, by)

                values(out)[, 'name'] = names(out)
                values(out)[, 'numintervals'] = table(by)[names(out)]
                
                tadt = gr2dt(targets)
                
                if (verbose)
                    cat('Aggregating columns \n')
                
                for (f in cfields)
                    {
                        if (verbose)
                            cat(f, '\n')
                        values(out)[, f] = tadt[, sum(eval(parse(text=f)), na.rm = TRUE), keyby = list(by = by)][names(out), V1]
                    }
            }
        else
            {
                if (is.na(rolling <- as.integer(rolling)))
                    stop('rolling must be a positive integer')

                if (is.na(rolling<=1))
                    stop('rolling must be a positive integer')

                if (verbose)
                    cat('Rolling using window of', rolling, '(output will be coordinate sorted)\n')
                
                tadt = gr2dt(sort(targets))

                tadt[, width := as.numeric(width)]

                tadt <- tadt[seqnames %in% c(seq(22), "X")]
                
                if ('count' %in% cfields ) {
                  print("rolling count")
                    out = tadt[, list(
                        count = rollapply(count, rolling, sum, na.rm = TRUE, fill = NA),
                        start = rollapply(start, rolling, min, fill = NA),
                        end = rollapply(end, rolling, max, fill = NA),
                        coverage = rollapply(coverage, rolling, sum, fill = NA)
                    ), by = seqnames]
                } else {
                    out = tadt[, list(
                        start = rollapply(start, rolling, min, fill = NA),
                        end = rollapply(end, rolling, max, fill = NA),
                        coverage = rollapply(coverage, rolling, sum, fill = NA)
                    ), by = seqnames]
                  }
                nna.ix = !is.na(out$start)
                
                if (!any(nna.ix))
                    stop('Malformed input, only NA ranges produced.  Reduce value of running')
                
                out = seg2gr(out[nna.ix])
                
                ## rolling weighted average, used below
                .rwa = function(v, w)
                    rollapply(v*w, rolling, sum, na.rm = TRUE, fill = NA)/rollapply(w*as.numeric(!is.na(v)), rolling, sum, na.rm = TRUE, fill = NA)
                
            }        
        
        fields = setdiff(fields, c('coverage', 'count', 'query.id'))
        for (f in fields)
            {
                if (verbose)
                    cat(f, '\n')

                if (is.null(rolling))
                    values(out)[, f] = tadt[, sum(width*eval(parse(text=f)), na.rm = TRUE)/sum(width[!is.na(eval(parse(text=f)))]), keyby = list(by = by)][names(out), V1]
                else  ## rolling weighted average
                    values(out)[, f] = tadt[, .rwa(eval(parse(text=f)), width), by = seqnames][, V1][nna.ix]

            }
        return(out)
    }


#' score.targets
#'
#' Scores targets based on covariates using gamma-poisson model with coverage as constant
#' 
#' @param targets annotated targets with fields $coverage, optional field, $count and additional numeric covariates
#' @return GRanges of scored results
#' @author Marcin Imielinski
#' @import GenomicRanges
#' @export
score.targets = function(targets, covariates = names(values(targets)),    
    model = NULL, ## fit existing model --> covariates must be present
    return.model = FALSE,
    nb = TRUE, ## negative binomial, if false then use poisson
    verbose = TRUE,
    iter = 200,
    subsample = 1e5,
    seed = NULL,
    p.randomized = TRUE,
    classReturn = FALSE)
    {
        
        require(MASS)
        require(data.table)
        covariates = setdiff(covariates, c('count', 'coverage', 'query.id'))        
        
        if (any(nnin <- !(covariates %in% names(values(targets)))))
            stop(sprintf('%s covariates (%s) missing from input data', sum(nnin), paste(covariates[nnin], collapse = ',')))
            
        if ( any(nnum <- !(sapply(covariates, function(x) class(values(targets)[, x])) %in% c('factor', 'numeric'))))
            {
                warning(sprintf('%s covariates (%s) fit are non-numeric or factor, removing from model', sum(nnum), paste(covariates[nnum], collapse = ',')))
                covariates = covariates[!nnum]
            }
        
        if (!all(c('count', 'coverage') %in% names(values(targets))))
            stop('Targets must have count, coverage, and query.id fields populated')

        if (verbose)
            cat('Setting up problem\n')

        values(targets)$count = round(values(targets)$count)

        if (length(unique(values(targets)$count))<=1)
            stop('score.targets input malformed --> count does not vary!')

        set.seed(42) ## to ensure reproducibility
        
        if (is.null(model))
            {
                tdt = as.data.table(as.data.frame(values(targets)[, c('count', 'coverage', covariates)]))
                tdt$coverage = log(tdt$coverage)

                if (subsample>nrow(tdt))
                    subsample = NULL

                tdt = tdt[rowSums(is.na(tdt[, c('count', 'coverage', covariates), with = FALSE]))==0,]

                if (nrow(tdt)==0)
                    stop('No rows with non NA counts, coverage, and covariates')

                if (!is.null(subsample))
                    {
                        if (subsample<1)
                            subsample = ceiling(pmax(0, subsample)*nrow(tdt))

                        if (verbose)
                            cat(sprintf('Subsampling ..\n'))

                        tdt = tdt[sample(1:nrow(tdt), subsample), ]
                    }
                
                if (verbose)
#                    cat(sprintf('Fitting model with %s data points and 10 covariates\n', prettyNum(nrow(tdt), big.mark = ','), length(covariates)))
                    cat(sprintf('Fitting model with %s data points and %s covariates\n', prettyNum(nrow(tdt), big.mark = ','), length(covariates)))
                formula = eval(parse(text = paste('count', " ~ ", paste(c('offset(1*coverage)', covariates), collapse = "+")))) ## make the formula with covariateso
                if (nb)
                    g = glm.nb(formula, data = as.data.frame(tdt), maxit = iter)
                else
                    {
                        g = glm(formula, data = as.data.frame(tdt), family = poisson)
                        g$theta = 1
                    }
            }
        else
            {
                g = model                
                #                covariates = setdiff(names(coefficients(g)), '(Intercept)')                
            }
                
        if(!(classReturn)){
            if (return.model)
                return(g)            
        }
        if (is(targets, 'GRanges'))
            res = as.data.frame(targets)
        else
            res = as.data.frame(values(targets))

        if ( any(is.fact <- (sapply(covariates, function(x) class(res[, x])) %in% c('factor'))))
            {
                ix = which(is.fact)
                new.col = lapply(ix, function(i)
                    {
                        val = res[, covariates[i]]
                        if (verbose)
                            cat('Factorizing column', covariates[i], 'with', length(val),
                                'across', length(levels(val)), 'levels\n')
                        tmp.mat = matrix(as.numeric(rep(val, each = length(levels(val))) == levels(val)),
                            ncol = length(levels(val)), byrow = TRUE)
                        colnames(tmp.mat) = paste(covariates[i], levels(val), sep = '')
                        return(tmp.mat)                        
                    })

                res = cbind(res[, -match(covariates[ix], names(res))], as.data.frame(do.call('cbind', new.col)))
                covariates = c(covariates[-ix], do.call('c', lapply(new.col, function(x) colnames(x))))
            }

        if (any(nnin <- !(covariates %in% names(res))))
            stop(sprintf('%s covariates (%s) missing from input data', sum(nnin), paste(covariates[nnin], collapse = ',')))
        
        coef = coefficients(g)
        na.cov = is.na(coef)
        
        if (any(na.cov))
            {
                warning(sprintf('%s covariates (%s) fit with an NA value, consider removing', sum(na.cov), paste(names(coef[na.cov]), collapse = ',')))
                covariates = setdiff(covariates, names(coef)[na.cov])
                coef = coef[-which(na.cov)]
            }
            
        ## if (verbose)
        ##     if (nb)
        ##         cat(sprintf('final model fit: \n  count ~ gamma-poisson(e^(%s), %s) \n',
        ##                     paste(signif(coef,2), names(coef), sep = '*', collapse = ' + '), signif(g$theta, 2)))
        ##     else
        ##         cat(sprintf('final model fit: \n  count ~ poisson(e^(%s)) \n',
        ##                     paste(signif(coef,2), names(coef), sep = '*', collapse = ' + ')))

        if (verbose)
            cat('Scoring results\n')

       
        M = cbind(1, as.matrix(res[, c('coverage', names(coef[-1])), drop = FALSE]))
            
        M[, 'coverage'] = log(M[, 'coverage'])
        res$count.pred = exp(M %*% c(coef[('(Intercept)')], 1, coef[colnames(M)[-c(1:2)]]))
        res$count.pred.density = res$count.pred / res$coverage
        res$count.density = res$count / res$coverage

        ## compute "randomized" p values (since dealing with counts data / discrete distributions
        if (nb)
            {
#                p = pnbinom(ifelse(res$count==0, -1e-15, res$count), mu = res$count.pred, size = g$theta, lower.tail = F)
                pval = pnbinom(res$count-1, mu = res$count.pred, size = g$theta, lower.tail = F)
                if (p.randomized)
                    {
                        pval.right = pnbinom(res$count, mu = res$count.pred, size = g$theta, lower.tail = F)
                        pval.right = ifelse(is.na(pval.right), 1, pval.right)
                        pval = ifelse(is.na(pval), 1, pval)
                        pval = runif(nrow(res), min = pval.right, max = pval)
                    }
                res$p = signif(pval, 2)
            }           
        else            
            {
                pval = ppois(res$count-1, lambda = res$count.pred, lower.tail = F)                
                if (p.randomized)
                    {
                        pval.right = ppois(res$count, lambda = res$count.pred, lower.tail = F)
                        pval.right = ifelse(is.na(pval.right), 1, pval.right)
                        pval = ifelse(is.na(pval), 1, pval)
                        pval = runif(nrow(res), min = pval.right, max = pval)
                    }                
                res$p = signif(pval, 2)
            }
        
        res$q = signif(p.adjust(res$p, 'BH'), 2)
        if (nb)
            res$p.neg = signif(pnbinom(res$count, mu = res$count.pred, size = g$theta, lower.tail = T), 2)
        else
            res$p.neg = signif(ppois(res$count, lambda = res$count.pred, lower.tail = T), 2)
        res$q.neg = signif(p.adjust(res$p.neg, 'BH'), 2)
        res$effectsize = log2(res$count / res$count.pred)

        if(!(classReturn)){
            return(as.data.table(res))
        }
        return(list(as.data.table(res),g))
    }

#' Cov
#'
#' Stores Covariate for passing to FishHook object. To be packaged in the Cov_Array Class by calling c(Cov1,Cov2,Cov3)
#' 
#' @param Covariate object of type, GRanges, ffTrack, RleList or character. Note that character objects must be paths to files containing one of the other types as a .rds file
#' @param type a string indicating the type of Covariate, valid options are: numeric, sequence, interval. See Annotate Targets for more information on Covariate types
#' @param signature In the case where a ffTrack object is of type sequence, a signature field is required, see fftab in ffTrack for more information.
#' @return Cov object that can be passed to FishHook object constructor
#' @author Zoran Z. Gajic
#' @importFrom R6 R6Class
#' @export
Cov <- R6Class("Cov",
                      public = list(
                          
                          initialize = function(Covariate = NULL, type = NULL, signature = NULL,
                                                name = "", pad = NULL, na.rm = NULL, field = NULL,
                                                grep = NULL, chr.sub = FALSE){

                              # Checks to see if covariates and type were supplied
                              if(is.null(Covariate) | is.null(type)){
                                  stop("Both Covariate and track must be non-null.")
                              }

                              # Checks to see if covariate type is one of the specified types in private$COV.TYPES
                              if(!(type %in% private$COV.TYPES)){
                                  stop('"type" must be "numeric", or "sequence", or "interval"')
                              }

                              # Checks to see if the class of the covariates is one of the specified classes in private$COV.CLASSES
                              if(!(class(Covariate) %in% private$COV.CLASSES)){
                                  stop(sprintf('Malformed covariate input: Covariate  must be of class GRanges, ffTrack, Rle, object, or character path to rds object with the latter or .bw, .bed file, $type must have value %s', paste(private$COV.TYPES, collapse = ',')))
                              }


                              ## Requires any name provided to be a character
                              if(!(is.character(name))){
                                  stop("Name must be of type 'character'")
                              }

                              
                              ## Sequence Covariates
                              if (type == 'sequence'){

                                  ## If the covariate is a path to a file
                                  if (is.character(Covariate)) {
                                      cov$track = tryCatch(readRDS(cov$track), error = function(e) 'not ffTrack')
                                  }

                                  ## Requires the covariate that was read above or provided in the initial arguements to be  an ffTrack.
                                  if (class(Covariate) != 'ffTrack'){
                                      stop('sequence tracks must have ffTrack object as $track field or $track must be a path to an ffTrack object rds file')
                                  }                                  
                              }

                              ## Checks to see that the signature was provided if using a sequence covariate
                              if (class(Covariate) == 'ffTrack' & type == 'sequence'){
                                  if (is.null(signature))
                                      stop('sequence tracks must be ffTracks and have a $signature field specified (see fftab in ffTrack)')                    
                              }

                              ## Assigns and initialized the Cov if all of the above is satisfied

                              if(class(Covariate) == "GRanges" & chr.sub){
                                  seqlevels(Covariate) = gsub('chr','',seqlevels(Covariate))
                              }
                              private$Covariate = Covariate
                              private$type = type
                              private$signature = signature
                              private$name = name
                              private$pad = pad
                              private$na.rm = na.rm
                              private$field = field
                              private$grep = grep
                          },

                          seqlevels = function(...){
                              if(class(private$Covariate) == "GRanges"){
                                  return(seqlevels(private$Covariate))
                              }
                              return (NULL)                              
                          },
                          
                          ## prints covariate to output
                          print = function(...){
                              cat(c("type: ",private$type,"\tsignature: ", private$signature,
                                    "\nfield: ",private$field, "\tpad: ", private$pad,
                                    "\nna.rm: ", private$na.rm, "\tgrep: ", private$grep,
                                    "\ntrack:\n", private$track),collapse = "")
                              print(private$Covariate)
                          },
                          ## Returns Covariate name (private$name)
                          getName = function(...){return(private$name)},
                          getCovariate = function(...){
                              return (private$Covariate)
                          },
                          ## Returns Covariate type (private$type)
                          getType = function(...){
                              return (private$Type)
                          },
                          ## Returns Covariate signature (private$ signature)
                          getSig = function(...){
                              return (private$signature)
                          },
                          ## Returns all of the valid covariate classes
                          getCovariateClasses = function(...){return (private$COV.CLASSES)},
                          ## Returns all of the valid covaraite types
                          getCovariateTypes = function(...){return (private$COV.TYPES)},

                          ## Checks to see if this covariate has seqlevels that begin with chr. must be a GRanges to do this
                          getChr = function(...){
                              if(class(private$Covariate) == "GRanges"){
                                  return(any(grepl("chr",seqlevels(private$Covariate))))
                              }
                              return (NULL)

                          },
                          
                          ## Converts a Cov object to a list that can be passed as input for annotate.targets
                          toList = function(...){
                              if(!(is.null(private$signature)) & class(private$Covariate) == "ffTrack"){
                                  return (list(track = private$Covariate, type = private$type,
                                               signature = private$signature,pad = private$pad,
                                               na.rm = private$na.rm,field = private$field,
                                               grep = private$grep))
                              }
                              else{                                 
                                  return (list(track = private$Covariate, type = private$type,
                                               signature = private$signature,pad = private$pad,
                                               na.rm = private$na.rm,field = private$field,
                                               grep = private$grep))

                              }                              
                          }                          
                          
                      ),## hcc1143

               
                      private = list(
                          
                          ## Active Covariates
                          Covariate = NULL,

                          ## Active Track
                          type = NULL,

                          ## signature for use with ffTrack sequence covariates
                          signature = NULL,

                          ## Pad for use with annotate targets
                          pad = NULL,

                          ## na.rm for use with annotate targets
                          na.rm = NULL,

                          ## field specifies the column of the covariate  to use for numeric covariates
                          field = NULL,

                          ## grep for use with sequence covariates
                          grep = NULL,
                          
                          ## Covariate name
                          name = NULL,
                          
                          ## Valid Covariate Types
                          COV.TYPES = c('numeric', 'sequence', 'interval'),
                          
                          ## Valid Covariate Classes
                          COV.CLASSES = c('GRanges', 'RleList', 'ffTrack', 'character')

                         
                      ),
                      

               
                      active = list()
                      )


#' c.Cov
#'
#' Override the c operator for covariates so that when you type: c(Cov1,Cov2,Cov3) it returns a Cov_Arr object that support vector like operation.
#' 
#' @param ... A series of Covariates, note all objects must be of type Cov
#' @return Cov_Arr object that can be passed directly into the FishHook object constructor
#' @author Zoran Z. Gajic
#' @export
'c.Cov' <- function(...){
    Covs = list(...)
    isc = sapply(Covs, function(x) class(x)[1] == "Cov")
    if(any(!isc))
        stop("All inputs must be of class Cov.")
    return(Cov_Arr$new(...))
}


#' Cov_Arr
#'
#' Stores Covariates for passing to FishHook object constructor.Standard initialization involves calling c(Cov1,Cov2,Cov3). Cov_Arr serves to mask the underlieing list implemenations of Covariates in the FishHook Object. This class attempts to mimic a vector in terms of subsetting and in the future will add more vector like operations.
#' @param ... several Cov objects for packaging.
#' @return Cov_Arr object that can be passed directly to the FishHook object constructor
#' @author Zoran Z. Gajic
#' @importFrom R6 R6Class
#' @export
Cov_Arr <- R6Class("Cov_Arr",

                   public = list(

                       ## Creates a Cov_Arr object. This function is noramlly called from using the c operator on Cov objects
                       ## e.g. class(c(Cov_1,Cov_2,Cov_3))[1] == Cov_Arr
                       initialize = function(...){
                           Covs = list(...)
                           ## Checks to make sure that all items within the c() operator are all of class Cov
                           isc = sapply(Covs, function(x) class(x)[1] == "Cov")
                           if(any(!isc))
                               stop("All inputs must be of class Cov.")
                           ## Extracts names from all Covs
                           names = sapply(Covs,function(x) x$getName())
                           names(Covs) = names
                           ## Assigning and intializing.
                           private$Covs = Covs
                           private$names = names

                       },

                       ## Returns the Covariates as a lsit
                       getArr = function(...){return (private$Covs)},


                       ## Returns a vector where TRUE indicates a chr based seqlevels
                       ## Note that non-GRanges Covariates will not return anything
                       getChr = function(...){
                           return(unlist(lapply(private$Covs, function(x) x$getChr())))
                       },
                       seqlevels = function(...){
                           return(lapply(private$Covs, function(x) x$seqlevels()))
                       },
                                                  
                       ## Returns a subset of the Covariates as a list
                       subset = function(range,...){
                           private$Covs = private$Covs[range]
                           private$names = private$names[range]
                       },

                       ## Creates a list of lists for passing to annotate.targets.
                       ## The inner list constitutes a Cov$toList()
                       ## The outer list serves to hold all of the Cov lists
                       toList = function(...){
                           if(length(private$Covs) == 0){
                               return(list())
                           }
                           out = lapply(private$Covs, function(x) x$toList())
                           names(out) = private$names
                           return(out)

                       },
                       print = function(...){print(self$getArr())}
                   ),
                   private = list(
                       Covs = list(),
                       names = c()
                   ),
                   active = list()
                   )


#' [.Cov_Arr
#'
#' Overrides the subset operator x[] for use with Cov_Arr to allow for vector like subsetting
#'
#' @param obj This is the Cov_Arr to be subset
#' @param range This is the range of Covariates to return, like subsetting a vector. e.g. c(1,2,3,4,5)[3:4] == c(3,4)
#' @return A new Cov_Arr object that contains only the Covs within the given range
#' @author Zoran Z. Gajic
#' @export
'[.Cov_Arr' <- function(obj,range){
    #return(a$getTargets()[b])
    # a = obj b = range
    ret = obj$clone()
    ret$subset(range)
    return (ret)
}



#' FishHook
#'
#' Stores Events, Targets, Eligible, Covariates. 
#'
#' @param targets Examples of targets are genes, enhancers, 1kb tiles of the genome that we can then convert into a rolling window. This param must be of class "GRanges".
#' @param events Events are the given mutational regions and must be of class "GRanges". Examples of events are mutational data (e.g. C->G) copy number variations and fusion events. Targets are the given regions of the genome to annotate and must be of class "GRanges". 
#' @param eligible Eligible are the regions of the genome that we feel are fit to score. For example in the case of exome sequencing where not all regions are equally represented, eligible can be a set of regions that meet an arbitrary coverage threshold. Another example of when to use eligibility is in the case of whole genomes, where your targets are 1kb tiles. Regions of the genome you would want to exclude in this case are highly repetative regions such as centromeres, telomeres, and satelite repeates. This param must be of class "GRanges".
#' @param covariates Covariates are genomic covariates that you belive will cause your given type of event (mutations, CNVs, fusions) that are not linked to the process you are investigating (e.g. cancer biology). In the case of cancer biology we are looking for regions that are mutated as part of cancer progression, and regions that are more suceptable to random mutagenesis such as late replicating or non-expressed region (transcription coupled repair) are potential false positives. Includinig covariates for these will reduce thier prominence in the final data. This param must be of type "Cov_Arr" which can be created by wrapping Cov objects in c(). e.g. c(Cov1,Cov2,Cov3).
#' @return FishHook object that can be annotated.
#' @author Zoran Z. Gajic
#' @importFrom R6 R6Class
#' @export
FishHook <- R6Class("FishHook",
                    public = list(
                        
                        initialize = function(targets = NULL, out.path = NULL, eligible = NULL, ... ,events = NULL, covariates = NULL){ 

                            ## This next portion checks to make sure that the seqlevels are in the same format
                            
                            ## Gets whther seqlevels of covariates are chr or not chr
                            seqLevelsStatus_Covariates = covariates$getChr()
                            ## Warns if there is a heterogenetiy of seqlevels (chr or not)
                            if(length(unique(seqLevelsStatus_Covariates)) > 1){
                                warning("Covariates appears to have mismatched seqlevels, make sure all Covariates have seqlevels that start with chr or don't", call.=TRUE)
                            }
                            ## gets the seqlevels and looks for chr to indicate USCS format
                            seqLevelsStatus_Targets = any(grepl("chr",seqlevels(targets)))
                            seqLevelsStatus_Eligible = any(grepl("chr",seqlevels(eligible)))
                            seqLevelsStatus_Events = any(grepl("chr",seqlevels(events)))

                            if(seqLevelsStatus_Targets != seqLevelsStatus_Covariates){
                                warning("seqlevels of Targets and Covariates appear to be in different formats")
                            }
                            if(seqLevelsStatus_Targets != seqLevelsStatus_Eligible){
                                warning("seqlevels of Targets and Eligible appear to be in different formats")
                            }   
                            if(seqLevelsStatus_Targets != seqLevelsStatus_Events){
                                warning("seqlevels of Targets and Events appear to be in different formats")
                            }
                            

                            ## This next portion checks to make sure there is atleast some overlap of seqlevels i.e. some mapability
                            if(!any(seqlevels(targets) %in% seqlevels(events))){
                                stop("Error, there are no seqlevels of events that match targets")
                            }
                            if(!any(seqlevels(targets) %in% seqlevels(eligible))){
                                stop("Error, there are no seqlevels of eligible that match targets")
                            }

                            
                            if(any(!(unlist(lapply(covariates$seqlevels(),function(x) any(x%in%seqlevels(targets))))))){                                
                                warning("Warning, atleast one of the covariates has no seqlevels in common with targets")
                            }
                            
                            

                            
                            ## Initializes and Validates targets                            
                            self$replaceTargets(targets)                            
                            
                            ## Initializes and Validates out.path
                            self$replaceOut.Path(out.path)
                            
                            ## Initializes and Validates covariates 
                            if(is.null(covariates)){
                                covariates = Cov_Arr$new()
                            }
                            if(class(covariates)[1] != "Cov_Arr"){
                                stop("Covariates must be a vector of covariates or an object of class 'Cov_Arr'")
                            }
                            private$covariates = covariates
                                                                
                            ## Initializes and Validates events 
                            self$replaceEvents(events)                            

                            ## Initializes and Validates eligible
                            if(!(is.null(eligible))){
                                self$replaceEligible(eligible)
                            }                            
                        },

                        
                        toList = function(x){
                            return (x$toList())
                        },

                        print = function(){
                            targ = paste( "Contains" , length(private$targets), "hypothesize." ,collapse = "")
                            eve = paste("Contains", length(private$events), "events to map to hypothesize.", collapse = "")
                            if(is.null(private$eligible)){
                                elig = "All regions are elgible."
                            }
                            else{
                                elig = "Will map only eliglble regions."
                            }
                            if(is.null(private$covariates)){
                                covs = "No covariates will be used."
                            }
                            else{                                
                                cov.names = names(private$covariates$getArr())
                                covs = cov.names
                                ## cov.class = sapply(private$covariates,
                                ##                    function(x) if (!is.null(x$class)) x$class else NA)
                                ## cov.types = sapply(private$covariates,
                                ##                    function(x) if (!is.null(x$type)) x$type else NA)
                                ## covs = paste(cov.names,cov.class,cov.types, sep = " ", collapse = "\n")

                            }
                            meta = paste("Targets contains", ncol(values(private$targets)), "metadata columns")
                            cat(targ,eve,elig,"Covariates:",covs,meta,sep = "\n",collapse = "\n")
                            
                            
                        },
                        
                        ## Access Functions 
                        getCovariates= function(...) {return (private$covariates)},
                        getTargets = function(...) {return (private$targets[,"name"])},
                        getEvents = function(...) {return (private$events)},
                        getEligible = function(...) {return (private$eligible)},
                        getOut.Path = function(...) {return (private$out.path)},
                        getMeta = function(...) { return (private$targets)},


                        ## Remove Functions Note that Events and Targets Cannot be Removed
#                        removeCovariates = function(...){private$covariates = NULL},
                        removeEligible = function(...){private$eligible = NULL},
                        removeOut.Path = function(...){private$out.path = NULL},
                        
                        
                        ## Replace Functions
                       

                        
                       
                       replaceTargets = function(targets = private$targets){                       
                            
                            ## checks if targets is NULL
                            if (is.null(targets))
                                stop('Targets cannot be "NULL".')
                            
                            ## checks to see if targets is a path & import if so
                            if (is.character(targets))
                                if (grepl('\\.rds$', targets[1]))
                                    targets = readRDS(targets[1])
                                else if (grepl('(\\.bed$)', targets[1]))
                                    targets = import.ucsc(targets[1])
                            
                            
                            ## Forces targets to be a GRanges Objects
                            if(class(targets) != "GRanges"){
                                stop('Loaded or provided class of targets must be "GRanges"')
                            }
                            
                            ## checks to see if target contains any data
                            if (length(targets)==0)
                                stop('Must provide non-empty targets')
                            
                            ## Looks for a "name" field to index/Identify targets by name
                            ## If no such field is found creates a set of indexes
                            if(is.null(targets$name)){
                                targets$name = 1:length(targets)
                            }
                            
                            private$targets = targets                            
                            
                        },

                       
                       replaceEvents = function(events = private$events){

                           ## Forces Events to Exist
                            if(is.null(events)){
                                stop('Events must exist and cannot be "NULL"')
                            }

                           ## Forces Events to be a GRanges  Object
                            if (class(events) != "GRanges" )
                                stop('events must be of class "GRanges"')

                           private$events = events
                            
                        },

                       
                       replaceEligible = function(eligible = private$eligible){
                           ## Forces Eligible to either not exist or to be a GRanges
                            if (class(eligible) != "GRanges" & !is.null(eligible))
                                stop('eligible must be of class "GRanges"')

                           private$eligible = eligible
                        },

                       replaceOut.Path = function(out.path = private$out.path){                                                        
                            ## checks to see if out.path is able to be written to
                            ## throws a warning if unable to write
                            if (!is.null(out.path))
                                tryCatch(saveRDS(private$targets, paste(gsub('.rds', '', out.path),
                                                                '.source.rds', sep = '')),
                                         error = function(e)
                                             warning(sprintf('Error writing to file %s', out.path)))

                            private$out.path = out.path
                       },

                       ## Takes a series of run params and passes them as well as
                       ## the internal data to the initialize function for the annotate object
                       ## Returns the created annotate object
                       annotateTargets = function(mc.cores = 1, na.rm = TRUE, pad = 0,
                                                  verbose = TRUE,max.slice = 1e3, ff.chunk = 1e6,
                                                  max.chunk = 1e11){
                           
                          res = Annotated$new(targets = private$targets,
                                               covered = private$eligible,
                                               events = private$events,
                                               mc.core = mc.cores,
                                               na.rm = na.rm,
                                               pad = pad,
                                               verbose = verbose,
                                               max.slice = max.slice,
                                               ff.chunk = ff.chunk,
                                               max.chunk = max.chunk,
                                               out.path = private$out.path,
                                               meta = private$targets,
                                               private$covariates$toList())
                           return(res)
                           
                       }

                       

                       
                       

                    ),

                    private = list(
                        ## Genomic Ranges Object that Indicates Hypotheses
                        targets = NULL,
                        
                        ## Eligible Regions for Targets
                        eligible = NULL,

                        ## Events to Count
                        events = NULL,

                        ## Covariates list for passing to fish.hook
                        covariates = NULL,

                        ## Potentially allow access to covariates via indexing

                        ## Valid Covariate Types
                        COV.TYPES = c('numeric', 'sequence', 'interval'),

                        ## Valid Covariate Classes
                        COV.CLASSES = c('GRanges', 'RleList', 'ffTrack', 'character'),

                        ## Path to write the annotated data to
                        ## Useful if working with long proccessing times due to
                        ## Many Covariates
                        out.path = NULL
                        
                    )
)


#' Annotated
#'
#' Stores the annotated data from a FishHook object. and allows users to aggregate,manipulate and score that data. This object should be generated by calling FishHook$annotateTargets(). Note that this is where the meat of the computational burden lies. For example, in our test cases, running 8k pts worth of exome seq on 20k genes took 20seconds without covariates and 20sec + ~5min per covariate added. 
#'
#'@param targets Examples of targets are genes, enhancers, 1kb tiles of the genome that we can then convert into a rolling window. This param must be of class "GRanges".
#' @param events Events are the given mutational regions and must be of class "GRanges". Examples of events are mutational data (e.g. C->G) copy number variations and fusion events. Targets are the given regions of the genome to annotate and must be of class "GRanges". 
#' @param covered This is equivalent to Eligible in the FishHook class. Eligible are the regions of the genome that we feel are fit to score. For example in the case of exome sequencing where not all regions are equally represented, eligible can be a set of regions that meet an arbitrary coverage threshold. Another example of when to use eligibility is in the case of whole genomes, where your targets are 1kb tiles. Regions of the genome you would want to exclude in this case are highly repetative regions such as centromeres, telomeres, and satelite repeates. This param must be of class "GRanges".
#' @param covariates Covariates are genomic covariates that you belive will cause your given type of event (mutations, CNVs, fusions) that are not linked to the process you are investigating (e.g. cancer biology). In the case of cancer biology we are looking for regions that are mutated as part of cancer progression, and regions that are more suceptable to random mutagenesis such as late replicating or non-expressed region (transcription coupled repair) are potential false positives. Includinig covariates for these will reduce thier prominence in the final data. This param must be of type "Cov_Arr" which can be created by wrapping Cov objects in c(). e.g. c(Cov1,Cov2,Cov3).
#' @return Annotate Obeject that can be scored & manipulated and aggregated.
#' @author Zoran Z. Gajic
#' @importFrom R6 R6Class
#' @export
Annotated <- R6Class("Annotate",
                     public = list(
                         ## Initialize takes in all of the same params as annotate.targets
                         ## Including a meta param that remembers the meta data associated with
                         ## the FishHook object.
                         ## Passes all params except meta to annotate.targets()
                         ## Makes a meta variable for rembering metadata
                         initialize = function(targets = NULL,
                                               covered = NULL,
                                               events = NULL,
                                               mc.cores = 1,
                                               na.rm = TRUE,
                                               pad = 0,
                                               verbose = TRUE,
                                               max.slice = 10000,
                                               ff.chunk = 1e6,
                                               max.chunk = 1e11,
                                               out.path = NULL,
                                               covariates ,
                                               meta = NULL){
                             private$annotated_targets = annotate.targets(targets = targets,
                                                                          covered = covered,
                                                                          events = events,
                                                                          mc.cores = mc.cores,
                                                                          na.rm = na.rm,
                                                                          pad = pad,
                                                                          verbose = verbose,
                                                                          max.slice = max.slice,
                                                                          ff.chunk = ff.chunk,
                                                                          max.chunk = max.chunk,
                                                                          out.path = out.path,
                                                                          covariates = covariates)

                             private$meta = meta
                         },
                         print = function(...){

                             anno = paste("Annotated has", length(private$annotated_targets), "hypothesize.")
                             if(class(private$annotated_targets) == "GRangesList"){
                                 agg = "The data has been aggregated."
                             }
                             else{
                                 agg = "The data has not been aggregated."
                             }
                             meta = paste("Annotated has", ncol(values(private$meta)), "metadata columns")
                             cat(anno,agg,meta,sep = "\n", collapse = "\n")




                         },

                         setAnno = function(anno_new){private$annotated_targets = anno_new},
                         ## Access Functions
                         getTargets = function(...) {
                             ## Checks if data has been aggregated, if so just return the GrangesList
                             if(class(private$annotated_targets) == "GRangesList"){
                                 return (private$annotated_targets)
                             }
                             
                             ## Adds "name" to targets for convienience if not already Aggregated
                             tmp = private$annotated_targets
                             tmp$name = private$meta$name
                             return (tmp)
                         },
                         getMeta =  function(...) {return (private$meta)},

                         ## Takes in a series of run parameters for the fish.hook function
                         ## score.targets and allows the user to change the meta data at
                         ## this point if they so desire, however this is not reccomended for
                         ## standard users.
                         ## Initializes a "Score" object using the provided params
                         scoreTargets = function(annotated = private$annotated_targets,
                                                 model = NULL,
                                                 return.model = FALSE,
                                                 nb = TRUE,
                                                 verbose = TRUE,
                                                 iter = 200,
                                                 subsample = 1e5,
                                                 seed = 42,
                                                 p.randomize = TRUE,
                                                 classReturn = TRUE,
                                                 meta = private$meta){

                             ## Assume that the data has not been aggregated
                             grl = FALSE
                             ## Checks if the data has been aggregated, if so:
                             ## get rid of the old meta as it will not pertain to the
                             ## output and you can get it from accessing the current object
                             ## set the control variable grl to TRUE
                             if(class(private$annotated_targets) == "GRangesList"){
                                 meta = NULL
                                 grl = TRUE
                             }

                             res = Score$new(annotated = annotated, model = model,
                                             return.model = return.model, nb = nb,
                                             verbose = verbose, iter = iter,
                                             subsample = subsample, seed = seed,
                                             p.randomize = p.randomize, meta = meta,grl = grl,
                                             classReturn = classReturn)
                             return(res)

                         },

                         ## Passes a series of parameteres to the fish.hook function "aggregate.targets"
                         ## The only required param here is by, which is a vector upon which to
                         ## group the targets in the GRanges object
                         ## See fish.hook documentation for more info on rolling and disjoint

                         aggregateTargets = function(by = NULL, fields = NULL, rolling = NULL,
                                                     disjoint = TRUE, na.rm = FALSE, FUN = list(),
                                                     verbose = TRUE){
                             ret = self$clone()
                             ret$aggregateTargets_sys(by = by, fields = fields, rolling = rolling,
                                                      disjoint = disjoint, na.rm = na.rm, FUN = FUN,
                                                      verbose = verbose)
                             ret$deleteMeta()
                             return(ret)
                             
                         },

                         aggregateTargets_sys = function(by = NULL, fields = NULL, rolling = NULL,
                                                     disjoint = TRUE, na.rm = FALSE, FUN = list(),
                                                     verbose = TRUE){
                             
                             ## Checks to make sure the data is not aggregated as aggreagting a grl
                             ## is much more complicated and can be done much more simply by breaking
                             ## up the grls into GRanges and aggregating each one, and then merging
                             ## the resultant grls into one large grl
                             if(class(private$annotated_targets) == "GRangesList"){
                                 stop("Data already aggregated")
                             }
                             
                             private$annotated_targets= aggregate.targets(private$annotated_targets,
                                                                          by = by, fields = fields,
                                                                          rolling = rolling,
                                                                          disjoint = disjoint,
                                                                          na.rm = na.rm,
                                                                          FUN = FUN,
                                                                          verbose = verbose)
                                                      
                             

                         },

                         subset = function(range){
                             size = length(private$annotated_targets)
                             private$annotated_targets = private$annotated_targets[range]
                             if (length(private$meta) != size){
                                 warning("Meta and Targets are of different length, removing Meta.")
                                 private$meta = NULL
                                 return()
                             }
                             else{
                                 private$meta = private$meta[range]
                                 return()
                             }


                         },

                         deleteMeta = function(){private$meta = NULL}
                         

               
                         
                         
                     ),
                     
                     private = list(
                         ## the direct output of the annotate.targets function run at initialization
                         ## Or the results of aggregate
                         ## Discriminator: if output of annotate.targets -> class(annotated_targets) = "GRanges
                         ## if result of aggregate -> class(annotated_targets) = "GRangesList
                         annotated_targets = NULL,

                         ## An extra param that remembers meta data for the targets
                         meta = NULL
                     ),
                     
                     active = list()
                     )

#' [.Annotate
#'
#' Overrides the "[" operator for the Annotated object. This allows subsetting of the annotated data in Annotated Objects.
#'
#' @param obj This is the Annotated object to be subset
#' @param range This is the range of targets to return, like subsetting a vector. e.g. c(1,2,3,4,5)[3:4] == c(3,4)
#' @return Annotated object that can manipulated and scored, but cannot be aggregated again.
#' @author Zoran Z. Gajic
#' @export
'[.Annotate' <- function(obj,range){
    #return(a$getTargets()[b])
    # a = obj b = range
    ret = obj$clone()
    ret$subset(range)
    return (ret)
}



#' Score
#'
#' Stores the scored targets. Note that this constructors should be called from Annotated$scoreTargets(). Scores can also be plotted on qqplots using included functions. For other params see score.targets()
#'
#' @param annotated The annotated targets as an output from Annotated$scoreTargets() or the standard score.targets().
#' @return Score object that can be plotted/analyzed
#' @author Zoran Z. Gajic
#' @importFrom R6 R6Class
#' @export
Score <- R6Class("Score",
                 public = list(

                     ## Initialize takes in all of the params for the fish.hook function
                     ## "score.targets" and an optional meta param that will be used to remember
                     ## meta data.
                     ## meta data is automatically initialized to the grl "numinterval" and "name"
                     ## if the passed annotated var is a grl.
                     ## This will generally be called by the Annotated object method "score()"
                     ## Passes all params except meta to score.targets from the fish.hook package
                     ## the meta param is stored internally to be accessed later
                     ## grl is a control variable that indicates if the passed "annotated" object is a grl
                     ## if true -> assume it is a grl, if false -> assume it is a GRanges
                     initialize = function(annotated = NULL,
                                           model = NULL,
                                           return.model = FALSE,
                                           nb = TRUE,
                                           verbose = TRUE,
                                           iter = 200,
                                           subsample = 1e5,
                                           seed = 42,
                                           p.randomize = TRUE,
                                           meta = NULL,
                                           classReturn = TRUE,
                                           grl = FALSE){

                         score_Output = score.targets(annotated,
                                                      covariates = names(values(annotated)),
                                                      return.model = return.model,
                                                      nb = nb,
                                                      verbose = verbose,
                                                      iter = iter,
                                                      subsample = subsample,
                                                      seed = seed,
                                                      classReturn = classReturn,
                                                      p.randomize = p.randomize)

                         private$score = score_Output[[1]]
                         private$model = score_Output[[2]]
                         ## Check if we are in grl mode, if so overwrite meta and set to "numintervals" and "name"
                         if(grl){
                             meta = as.data.table(cbind(private$score$numintervals,private$score$name))
                             colnames(meta) = c("numintervals","name")
                         }
                         ## Removes the name from the internal score variable as it is already in meta,
                         ## and when accessed via the getAll()
                         ## function, we pull the name from the meta.
                         ## This is to support grl functionality
                         private$score$name = NULL
                         private$meta = meta
                    
                     
                     },
                     print = function(...){
                         
                         anno = paste("Scored has", nrow(private$score), "hypothesize.")
                         if(class(private$meta) == "GRanges"){
                             meta = paste("Scored has", ncol(values(private$meta)), "metadata columns")
                         }
                         else{
                             meta = paste("Scored has", ncol(private$meta), "metadata columns")
                         }
                         cat(anno,meta,sep = "\n", collapse = "\n")

                         


                     },
                     ## Produces a plotly html output of the scored targets
                     ## plotly = FALSE will produce a standard R graph
                     ## You can assign metadata to annotations to plot
                     ## If no annotation is provided will use defaults
                     ## use annotation = list() to have no annotations
                     qq_plot = function(plotly = TRUE, annotations = NULL, ...){
                         
                        
                         if(is.null(annotations)){                      
                             res = self$getAll()
                             if(is.null(res$name)){
                                 names = c(1:nrow(res))
                             }
                             else{
                                 names = res$name
                             }
                             annotations = list(Hypothesis_ID = names,
                                                Count = res$count,
                                                Effectsize = round(res$effectsize,2),
                                                q = res$q)
                         }
                         else{
                             res = self$getAll()
                         }
                         return(qq_pval(res$p ,annotations = annotations,
                                        gradient = list(Count = res$count),
                                        titleText = "" ,  plotly = plotly))
                         
                     },

                     ## produces a html widjet and write to file ~/public_html/plot.html
                     ## render = function(annotations = NULL){
                     ##    wij(self$qq_plot(annotations = annotations))

                     ## },

                     ## Returns the scored targets as well as target names
                     getScore = function(...){
                         if(is.null(private$meta)){
                             return(private$score)
                         }                         
                                                    
                         tmp = private$score
                         tmp$name = private$meta$name
                         return(tmp)
               
                     },

                     Model = function(...){return (private$model)},

                     ## Returns the meta data with target names
                     getMeta = function(...){return(private$meta)},
                     
                     ## Returns both the scored data and the meta data with target names
                     getAll = function(...){
                         if(is.null(private$meta)){
                             return(private$score)
                         }
                         if(class(private$meta) == "GRanges"){
                             return(seg2gr(cbind(gr2dt(private$meta),values(seg2gr(private$score)))))
                         }
                         return(cbind(private$score,private$meta))
                     },

                     ## Alows the user to replace meta data.
                     ## This is mainly intended for use with GRangesLists as the internals
                     ## do not handle the GRangesLists metadata as some might want.
                     ## "meta" should be of type "data.table" or "data.frame" but data.table is prefereable
                     replaceMeta = function(meta = NULL){
                         ## if no meta is specified assume the user wants to delete it
                         if(is.null(meta)){
                             private$meta = NULL
                         }
                         ## enforce the data.table or data.frame typing
                         if(class(meta) != "data.frame"){
                             stop("replaceMeta takes a data.frame or a data.table as args")
                         }
                         ## make sure that the meta aligns with score
                         if(nrow(meta) != nrow(private$score)){
                             stop("metadata must be of same length as targets")
                         }
                         ## enforce a column "name" on meta that will identify the scored items
                         if(is.null(meta$name)){
                             stop('metadata must have a "name" column')
                         }
                         private$meta = meta
                         
                     }

                     
                 ),
                 private = list(
                     ## The output of the score.targets() function called during initialization
                     score = NULL,

                     ## The output model of score.targets()
                     model = NULL,
                     
                     ## Target Metadata
                     meta = NULL
                     
                 ),
                 active = list()
                 )




#' @name qq_pval
#' @title qq plot given input p values
#' @param obs vector of pvalues to plot, names of obs can be intepreted as labels
#' @param highlight optional arg specifying indices of data points to highlight (ie color red)
#' @param samp integer, optional specifying how many samples to draw from input data (default NULL)
#' @param lwd integer, optional, specifying thickness of line fit to data
#' @param pch integer dot type for scatter plot
#' @param cex integer dot size for scatter plot
#' @param conf.lines logical, optional, whether to draw 95 percent confidence interval lines around x-y line
#' @param max numeric, optional, threshold to max the input p values
#' @param label character vector, optional specifying which data points to label (obs vector has to be named, for this to work)
#' @param plotly toggles between creating a pdf (FALSE) or an interactive html widget (TRUE)
#' @param annotations named list of vectors containing information to present as hover text (html widget), must be in same order as obs input 
#' @param gradient named list that contains one vector that color codes points based on value, must bein same order as obs input 
#' @param titleText title for plotly (html) graph only
#' @author Marcin Imielinski, Eran Hodis, Zoran Z. Gajic
#' @import plotly
#' @export
qq_pval = function(obs, highlight = c(), exp = NULL, lwd = 1, bestfit=T, col = NULL, col.bg='black', pch=18, cex=1, conf.lines=T, max=NULL, max.x = NULL, max.y = NULL, qvalues=NULL, label = NULL, plotly = FALSE, annotations = list(), gradient = list(), titleText = "", subsample = NA, ...)
{
    if(!(plotly)){
        is.exp.null = is.null(exp)

        if (is.null(col))
            col = rep('black', length(obs))

        ix1 = !is.na(obs)
        if (!is.null(exp))
            if (length(exp) != length(obs))
                stop('length of exp must be = length(obs)')
            else
                ix1 = ix1 & !is.na(exp)

        if (is.null(highlight))
            highlight = rep(FALSE, length(obs))
        else if (is.logical(highlight))
            {
                if (length(highlight) != length(obs))
                    stop('highlight must be either logical vector of same length as obs or a vector of indices')
            }
        else
            highlight = 1:length(obs) %in% highlight

        obs = -log10(obs[ix1])
        col = col[ix1]

        highlight = highlight[ix1]
        if (!is.null(exp))
            exp = -log10(exp[ix1])

        ix2 = !is.infinite(obs)
        if (!is.null(exp))
            ix2 = ix2 &  !is.infinite(exp)

        obs = obs[ix2]
        col = col[ix2]

        highlight = highlight[ix2]
        if (!is.null(exp))
            exp = exp[ix2]

        N <- length(obs)
        ## create the null distribution
        ## (-log10 of the uniform)

        if (is.null(exp))
            exp <- -log(1:N/N,10)
        else
            exp = sort(exp)

        if (is.null(max))
            max = max(obs,exp) + 0.5
        
        if (!is.null(max) & is.null(max.x))
            max.x = max
        
        if (!is.null(max) & is.null(max.y))
            max.y  = max
        
        if (is.null(max.x))
            max.x <- max(obs,exp) + 0.5

        if (is.null(max.y))
            max.y <- max(obs,exp) + 0.5

        if (is.exp.null)
            {
                tmp.exp = rev(seq(0, 7, 0.01))
                ix = 10^(-tmp.exp)*N
                c95 <-  qbeta(0.975,ix,N-ix+1)
                c05 <-  qbeta(0.025,ix,N-ix+1)

                if (conf.lines){
                    ## plot the two confidence lines
                    plot(tmp.exp, -log(c95,10), ylim=c(0,max.y), xlim=c(0,max.x), type="l", axes=FALSE, xlab="", ylab="")
                    par(new=T)
                    plot(tmp.exp, -log(c05,10), ylim=c(0,max.y), xlim=c(0,max.x), type="l", axes=FALSE, xlab="", ylab="")
                    par(new=T)

                    p1 <- rep(tmp.exp[1], 2)
                    p2 <- c(-log(c95,10)[1], -log(c05,10)[1])

                    lines(x=p1, y=p2)
                    x.coords <- c(tmp.exp,rev(tmp.exp))
                    y.coords <- c(-log(c95,10),rev(-log(c05,10)))
                    polygon(x.coords, y.coords, col='light gray', border=NA)
                    par(new=T)
                }
            }

        ord = order(obs)

                                        #colors = vector(mode = "character", length = length(obs)); colors[] = "black";

        colors = col
        colors[highlight] = "red";

        dat = data.table(x = sort(exp), y = obs[ord], colors = colors[ord], pch = pch, cex = cex)
        if (!is.null(names(obs)))
            {
                names = names(obs[ord])
                setkey(dat, names)
            }

        if (nrow(dat)>1e5) ## rough guide to subsmapling the lower p value part of the plot
            subsample = 5e4/nrow(dat)

        if (is.na(subsample[1]))
            dat[, plot(x, y, xlab = expression(Expected -log[10](italic(P))), ylab = expression(Observed -log[10](italic(P))), xlim = c(0, max.x), col = colors, ylim = c(0, max.y), pch=pch, cex=cex, bg=col.bg, ...)]
        else
            {
                subsample = pmin(pmax(0, subsample[1]), 1)
                dat[ifelse(x<=2, ifelse(runif(length(x))<subsample, TRUE, FALSE), TRUE), plot(x, y, xlab = expression(Expected -log[10](italic(P))), ylab = expression(Observed -log[10](italic(P))), xlim = c(0, max.y), col = colors, ylim = c(0, max.y), pch=pch, cex=cex, bg=col.bg, ...)]
            }

        if (!is.null(label))
            {
                if (length(label)>0)
                    if (is.null(key(dat)))
                        warning('Need to provide names to input vector to draw labels')
                    else
                        dat[list(label), text(x, y, labels=label, pos=3)];
            }

        lines(x=c(0, max(max.y, max.x)), y = c(0, max(max.x, max.y)), col = "black", lwd = lwd)

        if (!is.na(subsample))
            dat = dat[sample(nrow(dat), subsample*nrow(dat)), ]

        lambda = lm(y ~ x-1, dat)$coefficients;

        lines(x=c(0, max.x), y = c(0, lambda*max.y), col = "red", lty = 2, lwd = lwd);
        legend('bottomright',sprintf('lambda=\n %.2f', lambda), text.col='red', bty='n')
    }

    else{
        
        if(length(annotations) < 1){
            hover <- do.call(cbind.data.frame, list(p = obs))
        }
        else{
            hover <- do.call(cbind.data.frame, list(annotations, p = obs))
        }    
        hover <- as.data.table(hover)

        
        is.exp.null = is.null(exp)
        if (is.null(col)) 
            col = "black"
        ix1 = !is.na(hover$p)
        if (!is.null(exp)) 
            if (length(exp) != length(hover$p)) 
                stop("length of exp must be = length(hover$obs)")
            else ix1 = ix1 & !is.na(exp)
        if (is.null(highlight)) 
            highlight = rep(FALSE, length(hover$p))
        else if (is.logical(highlight)) {
            if (length(highlight) != length(hover$p)) 
                stop("highlight must be either logical vector of same length as obs or a vector of indices")
        }
        else highlight = 1:length(hover$p) %in% highlight
        hover$obs = -log10(hover$p[ix1])
        hover = hover[ix1]
        highlight = highlight[ix1]
        if (!is.null(exp)) 
            exp = -log10(exp[ix1])
        ix2 = !is.infinite(hover$obs)
        if (!is.null(exp)) 
            ix2 = ix2 & !is.infinite(exp)
        hover = hover[ix2]
        highlight = highlight[ix2]
        if (!is.null(exp)) 
            exp = exp[ix2]
        N <- length(hover$obs)
        if (is.null(exp)) 
            exp <- -log(1:N/N, 10)
        else exp = sort(exp)
        if (is.null(max)) 
            max <- max(hover$obs, exp) + 0.5
        else max <- max
        if (is.exp.null) {
            tmp.exp = rev(seq(0, 7, 0.01))
            ix = 10^(-tmp.exp) * N
            c95 <- qbeta(0.975, ix, N - ix + 1)
            c05 <- qbeta(0.025, ix, N - ix + 1)
            if (FALSE) {   ##Don't need if not using conf.line (might put this in the future)
                plot(tmp.exp, -log(c95, 10), ylim = c(0, max), xlim = c(0, max),
                type = "l", axes = FALSE, xlab = "", ylab = "")

                par(new = T)
                plot(tmp.exp, -log(c05, 10), ylim = c(0, max), xlim = c(0, max),
                type = "l", axes = FALSE, xlab = "", ylab = "")

                par(new = T)
                p1 <- rep(tmp.exp[1], 2)
                p2 <- c(-log(c95, 10)[1], -log(c05, 10)[1])
                lines(x = p1, y = p2)
                x.coords <- c(tmp.exp, rev(tmp.exp))
                y.coords <- c(-log(c95, 10), rev(-log(c05, 10)))
                polygon(x.coords, y.coords, col = "light gray", border = NA)
                par(new = T)
            }
        }
                                  
        #creating the ploting data.table (dat) and organizing the annotations to create hover text
        ord = order(hover$obs)
        hover = hover[ord]
        dat = hover
        hover$obs = NULL

        #Creating the hover text
        if(length(colnames(hover)) > 1){                                   
            annotation_names  = sapply(colnames(hover), paste0, " : ")
            annotation_names_wLineBreak  = paste("<br>", annotation_names[2:length(annotation_names)],
            sep = "")
            annotation_names = c(annotation_names[1], annotation_names_wLineBreak)
        }
        else{
            annotation_names  = sapply(colnames(hover), paste0, " : ")
        }

        #Checking if there is a gradient and if so adding it to the plotting data.table (dat)
        gradient_control = FALSE
        if(length(gradient )!= 0){    
            dat$grad = gradient[[1]][ord]
            gradient_control = TRUE
        }
        else {   
            dat$grad = c()
        }
        
        
        dat$x = sort(exp)
        dat$y = dat$obs    
        
        #declare so we can use in If statement
        p <- NULL    

        #hacky subsampling but works really well, just maxing out the number of points at 8k
        #and removing the extra from the non-sig
        #(looks to be -logp of 2.6 here can make this more dynamic later )

        if (nrow(dat) <=  8000){

            dat4 = dat
            dat4$obs = NULL
            dat4$x = NULL
            dat4$y = NULL
            dat4$grad = NULL
            
            trans = t(dat4)
            hover_text = c()
            for (i in 1:dim(trans)[2]){
                outstr = paste(c(rbind(annotation_names, trans[,i])), sep = "", collapse = "")
                hover_text = c(hover_text,outstr)
            }

            if(gradient_control){
                p <- dat[, plot_ly(data = dat, x=x, y=y, hoverinfo = "text",text = hover_text, color = grad,
                                   colors = c("blue2","gold"),marker = list(colorbar = list(title = names(gradient[1]))),
                                   mode = "markers",type = 'scatter')
                    %>% layout(xaxis = list(title = "<i>Expected -log<sub>10</sub>(P)</i>"),
                               yaxis = list(title = "<i>Observed -log<sub>10</sub>(P)</i>")) ]
            }
            else{
                p <- dat[, plot_ly(data = dat, x=x, y=y, hoverinfo = "text",text = hover_text,
                                   mode = "markers",type = 'scatter')
                    %>% layout(xaxis = list(title = "<i>Expected -log<sub>10</sub>(P)</i>"),
                               yaxis = list(title = "<i>Observed -log<sub>10</sub>(P)</i>")) ]
            }
        }

        
        else {
                                    
            dat$ID = c(1:nrow(dat))
            dat2 = dat[ y < 2.6,]
            dat3 = as.data.frame(dat2)
            dat3 = as.data.table(dat3[ sample(nrow(dat3), min(4000,nrow(dat3))), ])
            dat2 = rbind(dat3,dat[!(ID%in%dat2$ID),])
            dat2$ID = NULL
            
            dat4 = dat2
            dat4$obs = NULL
            dat4$x = NULL
            dat4$y = NULL
            dat4$grad = NULL

            trans = t(dat4)
            hover_text = c()
            for (i in 1:dim(trans)[2]){
                outstr = paste(c(rbind(annotation_names, trans[,i])), sep = "", collapse = "")
                hover_text = c(hover_text,outstr)
            }
            
            if(gradient_control){
                p <- dat2[, plot_ly(data = dat2, x=x, y=y,hoverinfo = "text", text = hover_text, color = grad,
                                    colors = c("blue2","gold"),marker = list(colorbar = list(title = names(gradient[1]))),
                                    mode = "markers",type = 'scatter')
                     %>% layout(xaxis = list(title = "<i>Expected -log<sub>10</sub>(P)</i>"),
                                yaxis = list(title = "<i>Observed -log<sub>10</sub>(P)</i>")) ]
            }
            else{
                p <- dat2[,  plot_ly(data = dat2, x=x, y=y,hoverinfo = "text", text = hover_text,
                                    mode = "markers",type = 'scatter')
                     %>% layout(xaxis = list(title = "<i>Expected -log<sub>10</sub>(P)</i>"),
                                yaxis = list(title = "<i>Observed -log<sub>10</sub>(P)</i>")) ]
            }
            
        }

        
        #Calculating lambda, Note that this is using the whole data set not the subsampled one
        lambda = lm(y ~ x - 1, dat)$coefficients
        lambda_max = max*as.numeric(lambda)

        
        ##adding shapes (lines) + title  note that html <b></b> style is used for mods and plotting lines
        ##is done by specifying two points on the line (x0/y0 and x1/y1) 
        p <- layout(p,title = sprintf("<b>%s</b>" ,titleText),titlefont = list(size = 24),
                    shapes = list(list(type = "line",line = list(color = 'black'),
                    x0 = 0, x1  = max, xref = "x", y0 = 0, y1 = max,yref ="y"),
                    list( type = "line", line = list(color = "red"),
                    x0 = 0, x1 = max, xref = "x", y0 = 0, y1 = lambda_max, yref = "y")),
                    annotations = list(
                        x = (0.9 * max),
                        y = (0.03 * max),
                        text = paste("lambda =",sprintf("%.2f", signif(lambda,3)), collapse = " "),
                        font = list(
                            color = "red",
                            size = 20
                            ),
                        showarrow = FALSE,
                        xref = "x",
                        yref = "y"
                        ),
                    margin = list(
                        t = 100
                        
                        ),
                    hovermode = "compare")
    }
}

