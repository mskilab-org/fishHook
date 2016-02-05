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
#' @export
annotate.targets = function(targets, ## path to bed or rds containing genomic target regions with optional target name
    covered = NULL,
    events = NULL, ## events are in the 
    ...,
    mc.cores = 1,
    na.rm = TRUE,
    pad = 0,
    verbose = TRUE,
    max.slice = 1e3, ## max slice of intervals to evaluate with  gr.val
    ff.chunk = 1e6, ## max chunk to evaluate with fftab
    max.chunk = 1e11, ## gr.findoverlaps parameter
    out.path = NULL)
    {
        if (is.character(targets))
            if (grepl('\\.rds$', targets[1]))
                targets = readRDS(targets[1])
            else if (grepl('(\\.bed$)', targets[1]))
                targets = import.ucsc(targets[1])
        
        if (!is.null(out.path))
            tryCatch(saveRDS(targets, paste(gsub('.rds', '', out.path), '.source.rds', sep = '')), error = function(e) warning(sprintf('Error writing to file %s', out.file)))
        
        COV.TYPES = c('numeric', 'sequence', 'interval')
        COV.CLASSES = c('GRanges', 'RleList', 'ffTrack', 'character')
        covariates = list(...)
        
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
                        else ## assume it is an Rle
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
        
        ovdt = grdt(ov)
        
        
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
                
                tadt = grdt(targets)
                
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
                require(zoo)
                
                if (is.na(rolling <- as.integer(rolling)))
                    stop('rolling must be a positive integer')

                if (is.na(rolling<=1))
                    stop('rolling must be a positive integer')

                if (verbose)
                    cat('Rolling using window of', rolling, '(output will be coordinate sorted)\n')
                
                tadt = grdt(sort(targets))

                tadt[, width := as.numeric(width)]

                if ('count' %in% cfields )
                    out = tadt[, list(
                        count = rollapply(count, rolling, sum, na.rm = TRUE, fill = NA),
                        start = rollapply(start, rolling, min, fill = NA),
                        end = rollapply(end, rolling, max, fill = NA),
                        coverage = rollapply(coverage, rolling, sum, fill = NA)
                    ), by = seqnames]
                else
                    out = tadt[, list(
                        start = rollapply(start, rolling, min, fill = NA),
                        end = rollapply(end, rolling, max, fill = NA),
                        coverage = rollapply(coverage, rolling, sum, fill = NA)
                    ), by = seqnames]

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
#' @export
score.targets = function(targets, covariates = names(values(targets)),    
    model = NULL, ## fit existing model --> covariates must be present
    return.model = FALSE,
    nb = TRUE, ## negative binomial, if false then use poisson
    verbose = TRUE,
    iter = 200,
    subsample = 1e5,
    seed = NULL,
    p.randomized = TRUE)                         
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

        targets$count = round(targets$count)

        if (length(unique(targets$count))<=1)
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
                

        if (return.model)
            return(g)            

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
            
        if (verbose)
            if (nb)
                cat(sprintf('final model fit: \n  count ~ gamma-poisson(e^(%s), %s) \n',
                            paste(signif(coef,2), names(coef), sep = '*', collapse = ' + '), signif(g$theta, 2)))
            else
                cat(sprintf('final model fit: \n  count ~ poisson(e^(%s)) \n',
                            paste(signif(coef,2), names(coef), sep = '*', collapse = ' + ')))

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
            
        return(as.data.table(res))
    }
