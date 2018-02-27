#############################################################################
## Marcin Imielinski
## Weill Cornell Medicine  mai9037@med.cornell.edu
## New York Genome Center mimielinski@nygenome.org

## Zoran Gajic
## New York Genome Center zgajic@NYGENOME.ORG

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

#' @name annotate.targets
#' @title title
#' @description
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
#'
#' @param targets path to bed or rds containing genomic target regions with optional target name
#' @param covered  optional path to bed or rds containing  granges object containing "covered" genomic regions (default = NULL)
#' @param events  optional path to bed or rds containing ranges corresponding to events (ie mutations etc) (default = NULL)
#' @param mc.cores integer info (default = 1)
#' @param na.rm info (default = TRUE)
#' @param pad  info (default = 0)
#' @param verbose boolean verbose flag (default = FALSE)
#' @param max.slice integer Max slice of intervals to evaluate with  gr.val (default = 1e3)
#' @param ff.chunk integer Max chunk to evaluate with fftab (default = 1e6)
#' @param max.chunk integer gr.findoverlaps parameter (default = 1e11)
#' @param out.path  out.path to save variable to (default = NULL)
#' @param covariates list of lists where each internal list represents a covariate, the internal list can have elements: track, type,signature,name,pad,na.rm = na.rm,field,grep. See Cov_Arr class for descriptions of what
#' each of these elements do. Note that track is equivalent to the 'Covariate' parameter in Cov_Arr
#' @param maxpatientpergene Sets the maximum number of events a patient can contribute per target (default = Inf)
#' @param ptidcol string Column where patient ID is stored
#' @param weightEvetns boolean If TRUE, will weight events by their overlap with targets. e.g. if 10% of an event overlaps with a target
#' region, that target region will get assigned a score of 0.1 for that event. If false, any overlap will be given a weight of 1.
#' @param ... paths to sequence covariates whose output names will be their argument names, and each consists of a list with (default = FALSE)
#' $track field corresponding to a GRanges, RleList, ffTrack object (or path to rds containing that object), $type which can
#' have one of three values "numeric", "sequence", "interval".
#' Numeric tracks must have $score field if they are GRanges), and can have a $na.rm logical field describing how to treat NA values
#' (set to na.rm argument by default)
#' Sequence covariates must be ffTrack objects (or paths to ffTrack rds) and require an additional variables $signatures, which
#' will be used as input to fftab, and can have optional logical argument $grep to specify inexact matches (see fftab)
#' fftab signature: signatures is a named list that specify what is to be tallied.  Each signature (ie list element)
#' consist of an arbitrary length character vector specifying strings to %in% (grep = FALSE)
#' or length 1 character vector to grepl (if grep = TRUE)
#' or a length 1 or 2 numeric vector specifying exact value or interval to match (for numeric data)
#' Every list element of signature will become a metadata column in the output GRanges
#' specifying how many positions in the given interval match the given query
#'
#' Interval covariates must be Granges (or paths to GRanges rds) or paths to bed files
#' @return GRanges of input targets annotated with covariate statistics (+/- constrained to the subranges in optional argument covered)
#' @author Marcin Imielinski
annotate.targets = function(targets, covered = NULL, events = NULL,  mc.cores = 1, na.rm = TRUE, pad = 0, verbose = TRUE, max.slice = 1e3,
    ff.chunk = 1e6, max.chunk = 1e11, out.path = NULL, covariates = list(), maxpatientpergene = Inf, ptidcol = NULL, weightEvents = FALSE, ...)
{
    if(weightEvents){
        maxpatientpergene = NULL
    }

    if (is.character(targets)){
        if (grepl('\\.rds$', targets[1])){
            targets = readRDS(targets[1])
        } else if (grepl('(\\.bed$)', targets[1])){
            targets = import.ucsc(targets[1])
        }
    }

    if (length(targets)==0){
        stop('Error: Must provide non-empty targets')
    }

    if (!is.null(out.path)){
        tryCatch(saveRDS(targets, paste(gsub('.rds', '', out.path), '.source.rds', sep = '')), error = function(e) warning(sprintf('Error writing to file %s', out.file)))
    }

    COV.TYPES = c('numeric', 'sequence', 'interval')
    COV.CLASSES = c('GRanges', 'RleList', 'ffTrack', 'character')


    cov.types = sapply(covariates, function(x) if (!is.null(x$type)) x$type else NA)
    cov.classes = lapply(covariates, function(x) if (!is.null(x$track)) class(x$track) else NA)

    if (is.list(cov.types)){
        cov.type = ''
    }

    if (is.list(cov.classes)){
        cov.classes = ''
    }

    cov.types[is.na(cov.types)] = ''

    if (!all(cov.types %in% COV.TYPES) & !(all(cov.classes %in% COV.CLASSES))){
        stop(sprintf('Error: Malformed covariate input: each covariate must be a list with fields $tracks and $type, $track must be of class GRanges, ffTrack, Rle, object, or character path to rds object with the latter or .bw, .bed file, $type must have value %s', paste(COV.TYPES, collapse = ',')))
    }

    if (any(ix = (cov.types == 'sequence'))){
        for (cov in covariates[ix]){
            if (is.character(cov$track)){
                cov$track = tryCatch(readRDS(cov$track), error = function(e) 'Error: not ffTrack')
            }
            if (class(cov$track) != 'ffTrack'){
                stop('Error: Sequence tracks must have ffTrack object as $track field or $track must be a path to an ffTrack object rds file')
            }
        }
    }

    if (any(ix = (cov.classes == 'ffTrack' & cov.types == 'sequence'))){
        if (!all(sapply(covariates, function(x) !is.null(x$signature)))){
            stop('Error: Sequence tracks must be ffTracks and have a $signature field specified (see fftab in ffTrack)')
        }
    }

    if (verbose){
        cat('Overlapping with covered intervals\n')
    }

    if (!is.null(covered)){
        ov = gr.findoverlaps(targets, covered, verbose = verbose, max.chunk = max.chunk, mc.cores = mc.cores)
    } else {
        ov = targets[, c()]
        ov$query.id = ov$subject.id = 1:length(targets)
    }

    if (verbose){
        cat('Finished overlapping with covered intervals\n')
    }

    counts.unique = NULL

    if (length(ov) > 0){

        if (!is.null(events)){

            if (is(events, 'GRanges')){

                ev = gr.fix(events[gr.in(events, ov)])

                ## weighing each event by width means that each event will get one
                ## total count, and if an event is split between two tile windows
                ## then it will contribute a fraction of event proprotional to the number
                ## oof bases overlapping
                counts = coverage(ev, weight = 1/width(ev))
                oix = which(gr.in(ov, events))

                if(!is.null(maxpatientpergene)){

                    if(!is.numeric(maxpatientpergene)){
                        stop('Error: maxpatientpergene must be of type numeric')
                    }

                    if(!("ID" %in% colnames(values(events))) & is.null(ptidcol)){
                        events$ID = c(1:length(events))
                    }

                    ev2 = gr.findoverlaps(events,ov, max.chunk = max.chunk, mc.cores = mc.cores)

                    if(is.null(ptidcol)){
                        ev2$ID = events$ID[ev2$query.id]
                    } else{
                        ev2$ID = mcols(events)[,ptidcol][ev2$query.id]
                    }

                    ev2$target.id = ov$query.id[ev2$subject.id]
                    tab = as.data.table(cbind(ev2$ID, ev2$target.id))
                    counts.unique = tab[, dummy :=1][, .(count = sum(dummy)), keyby =.(V1, V2)][, count := pmin(maxpatientpergene, count)][, .(final_count = sum(count)), keyby = V2]
                }

            } else{
                ## assume it is an Rle of event counts along the genome
                counts = events
                oix = 1:length(ov)
            }

            if (verbose){
                cat('Computing event counts\n')
            }

            ov$count = 0

            if (length(oix)>0 & is.null(maxpatientpergene)){
                ov$count[oix] = fftab(counts, ov[oix], chunksize = ff.chunk, na.rm = TRUE, mc.cores = mc.cores, verbose = verbose)$score
            }

            if (!is.null(out.path)){
                tryCatch(saveRDS(ov, paste(out.path, '.intermediate.rds', sep = '')), error = function(e) warning(sprintf('Error writing to file %s', out.file)))
            }

            if (verbose){
                cat('Finished counting events\n')
            }
        }
    }






    for (nm in names(covariates)){

        cov = covariates[[nm]]

        if (verbose){
            cat('Starting track', nm, '\n')
        }

        if (cov$type == 'sequence'){

            if (is.null(cov$grep)){
                cov$grep = FALSE
            }

            if (verbose){
                cat('Starting fftab for track', nm, '\n')
            }

            if (!is.list(cov$signature)){
                cov$signature = list(cov$signature)
            }

            if (is.na(names(cov$signature))){
                if (length(cov$signature) > 1){
                    names(cov$signature) = 1:length(cov$signature)
                }

            }

            if (!is.na(names(cov$signature))){
                names(cov$signature) = paste(nm, names(cov$signature), sep = '.')
            } else{
                names(cov$signature) = nm
            }

            if (is.na(cov$pad)){
                cov$pad = pad
            }

            val = fftab(cov$track, ov + cov$pad, cov$signature, chunksize = ff.chunk, verbose = verbose, FUN = mean, na.rm = TRUE, grep = cov$grep, mc.cores = mc.cores)
            values(ov) = values(val)

            if (verbose){
                cat('Finished fftab for track', nm, '\n')
            }

            if (!is.null(out.path)){
                tryCatch(saveRDS(ov, paste(out.path, '.intermediate.rds', sep = '')), error = function(e) warning(sprintf('Error writing to file %s', out.file)))
            }
        } else if (cov$type == 'numeric'){
            if (is.character(cov$track)){
                if (grepl('.rds$', cov$track)){
                    cov$track = readRDS(cov$track)
                } else{
                    ## assume it is a UCSC format
                    require(rtracklayer)
                    cov$track = import(cov$track)
                }
            }

            if (is.na(cov$pad)){
                cov$pad = pad
            }
            if (is(cov$track, 'ffTrack') | is(cov$track, 'RleList')){
                val = fftab(cov$track, ov + cov$pad, signature = cov$signature, FUN = sum, verbose = verbose, chunksize = ff.chunk, grep = cov$grep, mc.cores = mc.cores)
                values(ov) = values(val)
            } else{
                if (is.na(cov$field)){
                    ## then must be GRanges
                    cov$field = 'score'
                }
                if (is.na(cov$na.rm)){

                }
                new.col = data.frame(val = values(gr.val(ov + cov$pad, cov$track, cov$field, mc.cores = mc.cores, verbose = verbose,  max.slice = max.slice, max.chunk = max.chunk, mean = TRUE, na.rm = cov$na.rm))[, cov$field])
                names(new.col) = nm
                values(ov) = cbind(values(ov), new.col)
            }

            if (!is.null(out.path)){
                tryCatch(saveRDS(ov, paste(out.path, '.intermediate.rds', sep = '')), error = function(e) warning(sprintf('Error writing to file %s', out.file)))
            }
        } else if (cov$type == 'interval'){

            if (is.character(cov$track)){

                if (grepl('.rds$', cov$track)){
                    cov$track = readRDS(cov$track)
                } else{
                    ## assume it is a UCSC format
                    require(rtracklayer)
                    cov$track = import(cov$track)
                }

            }

            if (is(cov, 'GRanges')){
                stop('Error: Interval tracks must be GRanges')
            }

            if (is.null(cov$pad)){
                cov$pad = pad
            }

            if (is.null(cov$na.rm)){
                cov$na.rm = na.rm
            }

            cov$track = reduce(cov$track)

            new.col = data.frame(val = gr.val(ov + cov$pad, cov$track[, c()], mean = FALSE, weighted = TRUE,  mc.cores = mc.cores, max.slice = max.slice, max.chunk = max.chunk, na.rm = TRUE)$value/(width(ov)+2*cov$pad))
            new.col$val = ifelse(is.na(new.col$val), 0, new.col$val)
            names(new.col) = nm
            values(ov) = cbind(values(ov), new.col)

            if (!is.null(out.path)){
                tryCatch(saveRDS(ov, paste(out.path, '.intermediate.rds', sep = '')), error = function(e) warning(sprintf('Error writing to file %s', out.file)))
            }
        }
    }
    
    ovdt = gr2dt(ov)
    
    
    cmd = 'list(coverage = sum(width), ';

    if (!is.null(events)){
        cmd = paste(cmd, 'count = sum(count)', sep = '')
    } else{
        cmd = paste(cmd, 'count = NA', sep = '')
    }


    cov.nm = setdiff(names(values(ov)), c('coverage', 'count', 'query.id', 'subject.id'))

    if (length(ov) > 0){

        if (length(cov.nm) > 0){
            cmd = paste(cmd,  ',', paste(cov.nm, '= mean(', cov.nm, ')', sep = '', collapse = ', '), ')',  sep = '')
        } else{
            cmd = paste(cmd, ')',  sep = '')
        }

        ovdta =  ovdt[, eval(parse(text = cmd)), keyby = query.id]
        values(targets) = as(as.data.frame(ovdta[list(1:length(targets)), ]), 'DataFrame')

        if(!is.null(maxpatientpergene)){
            targets$count = 0
            targets$count[as.numeric(counts.unique$V2)] = counts.unique$final_count
        }

    } else{
        targets$coverage = 0
    }

    targets$query.id = 1:length(targets)

    ix = is.na(targets$coverage)
    if (any(ix)){
        targets$coverage[ix] = 0
        if(!is.null(maxpatientpergene)){
            targets$count[targets$coverage == 0] = NA
        }
    }
    if (!is.null(out.path)){
        if (file.exists(paste(out.path, '.intermediate.rds', sep = ''))){
            system(paste('rm',  paste(out.path, '.intermediate.rds', sep = '')))  ## error catch above
        }
        tryCatch(saveRDS(targets, out.path), error = function(e) warning(sprintf('Error writing to file %s', out.file)))
    }

    return(targets)

}




#' @name aggregate.targets
#' @title title
#' @description
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
#' @param targets annotated GRanges of targets with fields $coverage, optional field, $count and additional numeric covariates, or path to .rds file of the same; path to bed or rds containing genomic target regions with optional target name
#' @param by character vector with which to split into meta-territories (default = NULL)
#' @param fields by default all meta data fields of targets EXCEPT reserved field names $coverage, $counts, $query.id (default = NULL)
#' @param rolling if specified, positive integer specifying how many (genome coordinate) adjacent to aggregate in a rolling fashion; positive integer with which to performa rolling sum / weighted average WITHIN chromosomes of "rolling" ranges" --> return a granges (default = NULL)
#' @param disjoint boolean only take disjoint bins of input (default = TRUE)
#' @param na.rm boolean only applicable for sample wise aggregation (i.e. if by = NULL) (default = FALSE)
#' @param FUN list only applies (for now) if by = NULL, this is a named list of functions, where each item named "nm" corresponds to an optional function of how to alternatively aggregate field "nm" per samples, for alternative aggregation of coverage and count.  This function is applied at every iteration of loading a new sample and adding to the existing set.   It is normally sum [for coverage and count] and coverage weighted mean [for all other covariates].  Alternative coverage / count aggregation functions should have two arguments (val1, val2) and all other alt covariate aggregation functions should have four arguments (val1, cov1, val2, cov2) where val1 is the accumulating vector and val2 is the new vector of values.
#' @param verbose boolean verbose flag (default = TRUE)
#' @return GRangesList of input targets annotated with new aggregate covariate statistics OR GRanges if rolling is specified
#' @author Marcin Imielinski
#' @import zoo
#' @importFrom data.table setkey := data.table as.data.table
#' @importFrom S4Vectors values values<-
#' @importFrom GenomeInfoDb seqnames
aggregate.targets = function(targets, by = NULL, fields = NULL, rolling = NULL, disjoint = TRUE, na.rm = FALSE, FUN = list(), verbose = TRUE)
{

    V1 = sn = st = en = keep = count = width = NULL ## NOTE fix
    if (is.null(by) & is.character(targets)){
        cat('Applying sample wise merging\n')
    } else if (is.null(by) & is.null(rolling)){
        stop('Error: argument "by" must be specified and same length as targets or "rolling" must be non NULL')
    }

    if (is.null(by) & is.character(targets)){

        if (!all(ix <- (file.exists(targets)) & grepl('\\.rds$', targets))){

            warning(sprintf('Warning: %s of the  %s input files for sample wise merging either do not exist or are not .rds files. Sample wise merging (i.e. when by is null) requires .rds files of equal dimension GRanges (same intervals, same meta data column names)', sum(!ix), length(ix)))
            if (sum(ix)==0){
                stop('No files to process')
            }
            targets = targets[ix]
        }

        out = readRDS(targets[1])
        gr = out
        if (is.null(out$coverage)){
            stop('Coverage missing for input targets')
        }


        core.fields = c('coverage', 'count', 'p', 'query.id')

        cfields = setdiff(names(values(out)), core.fields)

        if (!is.null(fields)){
            cfields = intersect(fields, cfields)
        }


        values(out) = values(out)[, intersect(c(core.fields, cfields), names(values(out)))]

        if (!is.null(out$p)){
            psum = 0
            psum.df = rep(0, length(out))
        }

        ### initialize everything to 0
        for (cf in cfields){
            values(out)[, cf] = 0
        }

        if (!is.null(out$count)){
            out$count = 0
        }

        out$coverage = 0
        out$numcases = length(targets)

        if (length(targets)>1)
            for (i in 1:length(targets)){

                if (verbose){
                    cat('Processing target file', targets[i], '\n')
                }

                if (i > 1){
                    gr = readRDS(targets[i])
                }

                if (!is.null(out$count)){
                    if (!is.null(FUN[['count']])){
                        out$count = do.call(FUN[['count']], list(out$count, gr$count))
                    } else{
                        out$count = as.numeric(out$count) + gr$count
                    }
                }


                for (cf in cfields){

                    if (cf %in% names(values(gr))){
                        val = as.numeric(values(gr)[, cf])
                    } else{
                        warning(paste(targets[i], 'missing column', cf))
                        val = NA
                    }

                    if (na.rm){
                        if (!is.null(FUN[[cf]])){
                            values(out)[, cf] = ifelse(!is.na(val), do.call(FUN[[cf]], list(values(out)[, cf], out$coverage, + val, gr$coverage)), values(out)[, cf])
                        } else{
                            values(out)[, cf] = ifelse(!is.na(val), (values(out)[, cf]*out$coverage + val*gr$coverage)/(out$coverage + gr$coverage), values(out)[, cf])
                        }
                    } else{
                        if (!is.null(FUN[[cf]])){
                            values(out)[, cf] = do.call(FUN[[cf]], list(values(out)[, cf], out$coverage, val, gr$coverage))
                        } else{
                            values(out)[, cf] = (values(out)[, cf]*out$coverage + val*gr$coverage)/(out$coverage + gr$coverage)
                        }
                    }
                }

                if (!is.null(out$p)){

                    if (!is.null(gr$p)){
                        has.val = is.na(gr$p)
                        psum = ifelse(has.val, psum - 2*log(gr$p), psum)
                        psum.df = ifelse(has.val, psum.df + 1, psum.df)
                    } else{
                        warning(paste(targets[i], 'missing p value column, ignoring for fisher combined computation'))
                    }
                }

                if (is.null(FUN[['coverage']])){
                    out$coverage = as.numeric(out$coverage) + gr$coverage
                } else{
                    out$coverage = do.call(FUN[['coverage']], list(as.numeric(out$coverage), gr$coverage))
                }


                if (!is.null(out$p)){
                    out$p = pchisq(psum, psum.df, lower.tail = FALSE)
                }
                return(out)
            }
        }

        if (is.null(fields)){
            fields = names(values(targets))
        }

        if (any(nnum <- !(sapply(setdiff(fields, 'query.id'), function(x) class(values(targets)[, x])) %in% 'numeric'))){
            warning(sprintf('%s meta data fields (%s) fit were found to be non-numeric and not aggregated', sum(nnum), paste(fields[nnum], collapse = ',')))
            fields = fields[!nnum]
        }


        cfields = intersect(names(values(targets)), c('coverage', 'count'))

        if (is.null(rolling)){

            by = as.character(cbind(1:length(targets), by)[,2])

            if (disjoint){
                tmp.sn = paste(by, seqnames(targets), sep = '_')
                tmp.dt = data.table(sn = paste(by, seqnames(targets), sep = '_'), st = start(targets), en = end(targets), ix = 1:length(targets))
                setkey(tmp.dt, sn, st)
                tmp.dt[, keep := c(TRUE, st[-1]>en[-length(st)]), by = sn]
                setkey(tmp.dt, ix)
                targets = targets[tmp.dt$keep, ]
                if (verbose){
                    cat(sprintf('Removing %s non-disjoint within group intervals, keeping %s\n', prettyNum(sum(!tmp.dt[, keep]), big.mark = ','), prettyNum(sum(tmp.dt[, keep]), big.mark = ',')))
                }
                by = by[tmp.dt$keep]
            }

            if (verbose){
                cat('Splitting into GRangesList\n')
            }

            out = split(targets, by)

            values(out)[, 'name'] = names(out)
            values(out)[, 'numintervals'] = table(by)[names(out)]

            tadt = gr2dt(targets)

            if (verbose){
                cat('Aggregating columns \n')
            }

            for (f in cfields){
                if (verbose){
                    cat(f, '\n')
                }
                values(out)[, f] = tadt[, sum(eval(parse(text=f)), na.rm = TRUE), keyby = list(by = by)][names(out), V1]
            }
        } else{
            ## check not NA
            if (is.na(as.integer(rolling))){
                stop('Error: rolling must be a positive integer')
            }

            if (rolling <= 1){
                stop('Error: rolling must be a positive integer greater than one')
            }

            if (verbose){
                cat('Rolling using window of', rolling, '(output will be coordinate sorted)\n')
            }

            tadt = gr2dt(sort(targets))

            tadt[, width := as.numeric(width)]

            tadt = tadt[seqnames %in% c(seq(22), "X")]

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

            if (!any(nna.ix)){
                stop('Error: Malformed input, only NA ranges produced. Reduce value of running')
            }

             out = seg2gr(out[nna.ix])

            ## rolling weighted average, used below
            .rwa = function(v, w){
                rollapply(v*w, rolling, sum, na.rm = TRUE, fill = NA)/rollapply(w*as.numeric(!is.na(v)), rolling, sum, na.rm = TRUE, fill = NA)
            }

        }

        fields = setdiff(fields, c('coverage', 'count', 'query.id'))

        for (f in fields){

            if (verbose){
                cat(f, '\n')
            }

            if (is.null(rolling)){
                values(out)[, f] = tadt[, sum(width*eval(parse(text=f)), na.rm = TRUE)/sum(width[!is.na(eval(parse(text=f)))]), keyby = list(by = by)][names(out), V1]
            } else{
                ## rolling weighted average
                values(out)[, f] = tadt[, .rwa(eval(parse(text=f)), width), by = seqnames][, V1][nna.ix]
            }

        }

    return(out)

}




#' @name score.targets
#' @title title
#' @description
#'
#' Scores targets based on covariates using Gamma-Poisson model with coverage as constant
#'
#' @param targets annotated targets with fields $coverage, optional field, $count and additional numeric covariates
#' @param covariates chracter vector, indicates which columns of targets contain the covariates
#' @param model fit existing model --> covariates must be present (default = NULL)
#' @param return.model boolean info (default = FALSE)
#' @param nb boolean If TRUE, uses negative binomial; if FALSE then use Poisson
#' @param verbose boolean verbose flag (default = TRUE)
#' @param iter integer info (default = 200)
#' @param subsample interger info (default = 1e5)
#' @param seed integer (default = 42)
#' @param p.randomized boolean Flag info (default = TRUE)
#' @param classReturn boolean Flag info (default = FALSE)
#' @return GRanges of scored results
#' @author Marcin Imielinski
#' @import GenomicRanges
score.targets = function(targets, covariates = names(values(targets)), model = NULL, return.model = FALSE, nb = TRUE,
    verbose = TRUE, iter = 200, subsample = 1e5, seed = 42, p.randomized = TRUE, classReturn = FALSE)
{
    require(MASS)
    covariates = setdiff(covariates, c('count', 'coverage', 'query.id'))

    if (any(nnin = !(covariates %in% names(values(targets))))){
        stop(sprintf('Error: %s covariates (%s) missing from input data', sum(nnin), paste(covariates[nnin], collapse = ',')))
    }

    if (any(nnum = !(sapply(covariates, function(x) class(values(targets)[, x])) %in% c('factor', 'numeric')))){
        warning(sprintf('%s covariates (%s) fit are non-numeric or factor, removing from model', sum(nnum), paste(covariates[nnum], collapse = ',')))
        covariates = covariates[!nnum]
    }

    if (!all(c('count', 'coverage') %in% names(values(targets)))){
        stop('Error: Targets must have count, coverage, and query.id fields populated')
    }

    if (verbose){
        cat('Setting up problem\n')
    }

    values(targets)$count = round(values(targets)$count)

    if (length(unique(values(targets)$count)) <= 1){
        stop('Error: "score.targets" input malformed --> count does not vary!')
    }

    set.seed(seed) ## to ensure reproducibility

    if (is.null(model)){

        tdt = as.data.table(as.data.frame(values(targets)[, c('count', 'coverage', covariates)]))
        tdt$coverage = log(tdt$coverage)

        if (subsample > nrow(tdt)){
            subsample = NULL
        }

        tdt = tdt[rowSums(is.na(tdt[, c('count', 'coverage', covariates), with = FALSE]))==0,]

        if (nrow(tdt)==0){
            stop('Error: No rows with non NA counts, coverage, and covariates')
        }

        if (!is.null(subsample)){

            if (subsample < 1){
                subsample = ceiling(pmax(0, subsample)*nrow(tdt))
            }

            if (verbose){
                cat(sprintf('Subsampling ..\n'))
            }

            tdt = tdt[sample(1:nrow(tdt), subsample), ]
        }

        if (verbose){
            cat(sprintf('Fitting model with %s data points and %s covariates\n', prettyNum(nrow(tdt), big.mark = ','), length(covariates)))
        }

        formula = eval(parse(text = paste('count', " ~ ", paste(c('offset(1*coverage)', covariates), collapse = "+")))) ## make the formula with covariates

        if (nb){
            g = glm.nb(formula, data = as.data.frame(tdt), maxit = iter)
        } else{
            g = glm(formula, data = as.data.frame(tdt), family = poisson)
            g$theta = 1
        }
    } else{
        g = model
    }

    if(!(classReturn)){
        if (return.model){
            return(g)
        }
    }

    if (is(targets, 'GRanges')){
        res = as.data.frame(targets)
    } else{
        res = as.data.frame(values(targets))
    }

    if (any(is.fact = (sapply(covariates, function(x) class(res[, x])) %in% c('factor')))){
        ix = which(is.fact)
        new.col = lapply(ix, function(i){
            val = res[, covariates[i]]

            if (verbose){
                cat('Factorizing column', covariates[i], 'with', length(val), 'across', length(levels(val)), 'levels\n')
            }

            tmp.mat = matrix(as.numeric(rep(val, each = length(levels(val))) == levels(val)), ncol = length(levels(val)), byrow = TRUE)
            colnames(tmp.mat) = paste(covariates[i], levels(val), sep = '')
            return(tmp.mat)
        })

        res = cbind(res[, -match(covariates[ix], names(res))], as.data.frame(do.call('cbind', new.col)))
        covariates = c(covariates[-ix], do.call('c', lapply(new.col, function(x) colnames(x))))

    }

    if (any(nnin = !(covariates %in% names(res)))){
        stop(sprintf('Error: %s covariates (%s) missing from input data', sum(nnin), paste(covariates[nnin], collapse = ',')))
    }

    coef = coefficients(g)
    na.cov = is.na(coef)

    if (any(na.cov)){
        warning(sprintf('Warning: %s covariates (%s) fit with an NA value, consider removing', sum(na.cov), paste(names(coef[na.cov]), collapse = ',')))
        covariates = setdiff(covariates, names(coef)[na.cov])
        coef = coef[-which(na.cov)]
    }

    if (verbose){
        cat('Scoring results\n')
    }

    M = cbind(1, as.matrix(res[, c('coverage', names(coef[-1])), drop = FALSE]))

    M[, 'coverage'] = log(M[, 'coverage'])
    res$count.pred = exp(M %*% c(coef[('(Intercept)')], 1, coef[colnames(M)[-c(1:2)]]))
    res$count.pred.density = res$count.pred / res$coverage
    res$count.density = res$count / res$coverage

    ## compute "randomized" p values (since dealing with counts data / discrete distributions
    if (nb){
        pval = pnbinom(res$count-1, mu = res$count.pred, size = g$theta, lower.tail = F)
        if (p.randomized){
            pval.right = pnbinom(res$count, mu = res$count.pred, size = g$theta, lower.tail = F)
            pval.right = ifelse(is.na(pval.right), 1, pval.right)
            pval = ifelse(is.na(pval), 1, pval)
            pval = runif(nrow(res), min = pval.right, max = pval)
        }
        res$p = signif(pval, 2)
    } else{
        pval = ppois(res$count-1, lambda = res$count.pred, lower.tail = F)
        if (p.randomized){
            pval.right = ppois(res$count, lambda = res$count.pred, lower.tail = F)
            pval.right = ifelse(is.na(pval.right), 1, pval.right)
            pval = ifelse(is.na(pval), 1, pval)
            pval = runif(nrow(res), min = pval.right, max = pval)
        }
        res$p = signif(pval, 2)
    }

    res$q = signif(p.adjust(res$p, 'BH'), 2)
    if (nb){
       res$p.neg = signif(pnbinom(res$count, mu = res$count.pred, size = g$theta, lower.tail = T), 2)
    } else{
        res$p.neg = signif(ppois(res$count, lambda = res$count.pred, lower.tail = T), 2)
    }
    res$q.neg = signif(p.adjust(res$p.neg, 'BH'), 2)
    res$effectsize = log2(res$count / res$count.pred)

    if(!(classReturn)){
        return(as.data.table(res))
    }

    return(list(as.data.table(res),g))

}







#' @name Cov_Arr
#' @title title
#' @description
#'
#' Stores Covariates for passing to FishHook object constructor.
#'
#' Can also be initiated by passing a vector of multiple vectors of equal length, each representing one of the internal variable names
#' You must also include a list containg all of the covariates (Granges, chracters, RLELists, ffTracks)
#'
#' Cov_Arr serves to mask the underlieing list implemenations of Covariates in the FishHook Object.
#' This class attempts to mimic a vector in terms of subsetting and in the future will add more vector like operations.
#'
#'
#' @param name character vector Contains names of the covariates to be created, this should not include the names of any Cov objects passed
#' @param pad numeric vector Indicates the width to extend each item in the covarite. e.g. if you have a GRanges covariate with two ranges (5:10) and (20:30) with a pad of 5,
#' These ranges wil become (0:15) and (15:35)
#' @param type character vector Contains the types of each covariate (numeric, interval, sequencing)
#' @param signature, see ffTrack, a vector of signatures for use with ffTrack sequence covariates
#' fftab signature: signatures is a named list that specify what is to be tallied.  Each signature (ie list element)
#' consist of an arbitrary length character vector specifying strings to %in% (grep = FALSE)
#' or length 1 character vector to grepl (if grep = TRUE)
#' or a length 1 or 2 numeric vector specifying exact value or interval to match (for numeric data)
#' Every list element of signature will become a metadata column in the output GRanges
#' specifying how many positions in the given interval match the given query
#' @param field, a chracter vector for use with numeric covariates (NA otherwise) the indicates the column containing the values of that covarites.
#' For example, if you have a covariate for replication timing and the timings are in the column 'value', the parameter field should be set to the character 'Value'
#' @param na.rm, logical vector that indicates whether or not to remove nas in the covariates
#' @param grep, a chracter vector of  grep for use with sequence covariates of class ffTrack
#' The function fftab is called during the processing of ffTrack sequence covariates grep is used to specify inexact matches (see fftab)
#' @param cvs, a list of covariates that can include any of the covariate classes (GRanges, ffTrack, RleList, character)
#' @return Cov_Arr object that can be passed directly to the FishHook object constructor
#' @author Zoran Z. Gajic
#' @import R6
#' @export
Cov_Arr = R6::R6Class('Cov_Arr',
    public = list(

    ## See the class documentation
    initialize = function(..., name = NA, cvs = NULL, pad = 0, type = NA, signature = NA, field = NA, na.rm = NA, grep = NA){
        
        ##If cvs are valid and are a list of tracks concatenate with any premade covs
        if(!is.null(cvs)){
            if(class(cvs) != 'list'){
                cvs = list(cvs)
            }            
            self$cvs = c(cvs)
        }
        ##Otherwise assume that no cvs are given
        else{
            return(self)
        }
        
        self$names = name
        self$type = type
        self$signature = signature
        self$field = field
        self$pad = pad
        self$na.rm = na.rm
        self$grep = grep
            

    },

    ## Params:
    ## ... Other Cov_Arrs to be merged into this array, note that it can be any number of Cov_Arrs
    ## Return:
    ## A single Cov_Arr object that contains the contents of self and all passed Cov_Arrs
    ## UI:
    ## None
    ## Notes:
    ## This is linked to the c.Cov_Arr override and the c.Cov_Arr override should be used preferentially over this
    merge = function(...){
        return (c(self,...))
    },

    ## Params:
    ## No params required, included arguments will be ignored.
    ## Return:
    ## a logical vector where each element corresponds to a covariate and where TRUE indicates a chr based seqlevels e.g. chr14, False -> 14
    ## UI:
    ## None
    ## Note:
    ## non-GRanges Covariates will not return NA
    chr = function(...){
        if(length(private$pCovs) == 0){
            return(NULL)
        }
        chrs = lapply(c(1:length(private$pCovs)), function(x){
            if(class(private$pCovs[[x]]) == 'GRanges'){
                return(any(grepl('chr',  GenomeInfoDb::seqlevels(private$pCovs[[x]]))))
            } else{
                return(NA)
            }
        })
        return(unlist(chrs))
    },


    ## Params:
    ## No params required, included arguments will be ignored.
    ## Return:
    ## returns a list of character vectors. If the respective covariate is of class GRanges, the vector will contain all of the chromosome names,
    ## if it is not of class GRanges, will return NA
    ## UI:
    ## None
    seqlevels = function(...){
        if(length(private$pCovs) == 0){
            return(NULL)
        }
        seqs = lapply(c(1:length(private$pCovs)), function(x){
            cov = private$pCovs[[x]]
            if(class(cov) == 'GRanges'){
                return(GenomeInfoDb::seqlevels(cov))
            } else{
                return(NA)
            }
        })
        return(seqs)
    },

    ## Params:
    ## range, a numeric vector of the covariates to include. e.g. if the Cov_Arr contains the covariates (A,B,C) and the range is c(2:3),
    ## this indicates you wish to get a Cov_Arr containing (B,C). NOTE THAT THIS DOES NOT RETURN A NEW COV_ARR, IT MODIFIES THE CURRENT.
    ## Return:
    ## None, this modifies the Cov_Arr on which it was called
    ## UI:
    ## None
    ## Notes:
    ## If you want to create a new Cov_Arr containing certain covariates, use the '[' operator, e.g. Cov_Arr[2:3]
    subset = function(range, ...){
        
        private$pCovs = private$pCovs[range]
        private$pnames = private$pnames[range]
        private$ptype = private$ptype[range]
        private$psignature = private$psignature[range]
        private$pfield = private$pfield[range]
        private$ppad = private$ppad[range]
        private$pna.rm = private$pna.rm[range]
        private$pgrep = private$pgrep[range]
    },

    ## Params:
    ## No params required, included arguments will be ignored.
    ## Return:
    ## A list of lists where each internal list corresponds to the covariate and is for use internally in the annotate.targets function
    ## The list representation of the covariate will contain the following variables: type, signature, pad, na.rm, field, grep
    ## UI:
    ## None
    toList = function(...){
        if(length(private$pCovs) == 0){
            return(list())
        }
        out = lapply(c(1:length(private$pCovs)), function(x){
            if(!(is.na(private$psignature[x])) & class(private$pCovs[[x]]) == 'ffTrack'){
                return (list(track = private$pCovs[[x]], type = private$ptype[x],
                    signature = private$psignature[x],
                    pad = private$ppad[x],
                    na.rm = private$pna.rm[x],
                    field = private$pfield[x],
                    grep = private$pgrep[x]))
            } else{
                return (list(track = private$pCovs[[x]],
                    type = private$ptype[x],
                    signature = private$psignature[x],
                    pad = private$ppad[x],
                    na.rm = private$pna.rm[x],
                    field = private$pfield[x],
                    grep = private$pgrep[x]))
            }
        })
        names(out) = private$pnames
        return(out)

        },

    ## Params:
    ## No params required, included arguments will be ignored.
    ## Return:
    ## Nothing
    ## UI:
    ## Prints information about the Cov_Arr to the console with all of covariates printed in order with variables printed alongside each covariate
    print = function(...){
        if(length(private$pCovs) == 0){
            cat('Empty Cov_Arr Object\n')
            return(NULL)
        }
        
        out= sapply(c(1:length(private$pCovs)),
            function(x){
                cat(c('Covariate Number: ' , x, '\nName: ', private$pnames[x],
                '\ntype: ',private$ptype[x], '\tsignature: ', private$psignature[x],
                '\nfield: ',private$pfield[x], '\tpad: ', private$ppad[x],
                '\nna.rm: ', private$pna.rm[x], '\tgrep: ', private$pgrep[x],
                '\nCovariate Class: ', class(private$pCovs[[x]]), '\n\n'), collapse = '', sep = '')
        })
    }

    

    ),

    ## Private variables are internal variables that cannot be accessed by the user
    ## These variables will have active representations that the user can interact with the update
    ## and view these variables, all internal manipulations will be done with these private variables
        private = list(
            ## The list of covariates, each element can be of class: 'GRanges', 'character', 'RleList', 'ffTrack'
            pCovs = list(),
            ## A string vector containing the names of the covariates, the covariate will be refered to by its name in the final table
            pnames = c(),
            ## Type is a string vector of types for each covariate, can be: 'numeric','sequence', or 'interval'
            ptype = c(),
            ## A vector of signatures for use with ffTrack, se fftab
            psignature = c(),
            ## A character vector of field names for use with numeric covariates, see the Cov_Arr class definition for more info
            pfield = c(),
            ## A numeric vector of paddings for each covariate, see the 'pad' param in Cov_Arr class definition for more info
            ppad = c(),
            ## A logical vector for each covariate, see the 'na.rm' param in Cov_Arr class definition for more info
            pna.rm = c(),
            ##  A chracter vector for each covariate, see the 'grep' param in Cov_Arr class definition for more info
            pgrep = c(),

            ##  Valid Covariate Types
            COV.TYPES = c('numeric', 'sequence', 'interval', NA),

            ##  Valid Covariate Classes
            COV.CLASSES = c('GRanges', 'RleList', 'ffTrack', 'character')

        ),

    ## The active list contains a variable for each private variable.
    ## Active variables are for user interaction,
    ## Interactions can be as such
    ## class$active will call the active variable function with the value missing
    ## class$active = value will call the active variable function with the value = value
        active = list(

            ## Covariate Names
            ## Here we check to make sure that all names are of class chracter and that they are the same length as pCovs -> the internal covariate list
            ## If the names vector is equal in length to the pCovs list length we will allow the assignment
            ## If the pCovs list length is a multiple of the names vector, will will replicate the names vector such that it matches in length
            ## to the pCovs list.
            names = function(value) {

                if(!missing(value)){

                    if(!is.character(value) && !all(is.na(value)) ){
                        stop('Error: names must be of class character')
                    }

                    if(length(value) != length(private$pCovs) & length(private$pCovs) %% length(value) != 0){
                        stop('Error: Length of names must be of length equal to the number of Covariates or a divisor of number of covariates.')
                    }

                    if(length(private$pCovs) / length(value) != 1){
                        private$pnames = rep(value, length(private$pCovs)/length(value))
                       return(private$pnames)
                    }

                    private$pnames = value
                    return(private$pnames)

                } else{
                    return(private$pnames)
                }
            },

            ## Covariate type
            ## Here we check to make sure that all types are of class chracter and that they are the same length as pCovs -> the internal covariate list
            ## We then check to make sure that each type is a valid type as defined by the COV.TYPES private parameter
            ## If the types vector is equal in length to the pCovs list length we will allow the assignment
            ## If the pCovs list length is a multiple of the types vector, will will replicate the types vector such that it matches in length
            ## to the pCovs list.
            type = function(value) {

                if(!missing(value)){
                    if(!is.character(value) && !all(is.na(value))){
                        stop('Error: type must be of class character')
                    }

                    if(length(value) != length(private$pCovs) & length(private$pCovs) %% length(value) != 0){
                        stop('Error: Length of type must be of length equal to the number of Covariates or a divisor of number of covariates.')
                    }

                    if(!all(value %in% private$COV.TYPES)){
                        stop('"type" must be "numeric", "sequence", or "interval"')
                    }

                    if(length(private$pCovs) / length(value) != 1){
                        private$ptype = rep(value, length(private$pCovs)/length(value))
                        return(private$ptype)
                    }

                    private$ptype = value
                    return(private$ptype)

                } else{
                    return(private$ptype)
                }
            },

            ##Covariate Signature
            ## Here we check to make sure that all signatures are list within lists  and that they are the same length as pCovs -> the internal covariate list
            ## If the signature vector is equal in length to the pCovs list length we will allow the assignment
            ## If the pCovs list length is a multiple of the signature vector, will will replicate the signature vector such that it matches in length
            ## to the pCovs list.
            signature = function(value) {
                if(!missing(value)){
                    if(!all(is.na(value))){
                        stop('Error: signature must be of class list')
                    }
                    if(length(value) != length(private$pCovs) & length(private$pCovs) %% length(value) != 0){
                        stop('Error: Length of signature must be of length equal to the number of Covariates or a divisor of number of covariates.')
                    }
                    if(length(private$pCovs) / length(value) != 1){
                        private$psignature = rep(value, length(private$pCovs)/length(value))
                        return(private$psignature)
                    }

                    private$psignature = value
                    return(private$psignature)

                } else{
                    return(private$psignature)
                }
            },

            ##Covariate Field
            ## Here we check to make sure that all fields are of class chracter and that they are the same length as pCovs -> the internal covariate list
            ## If the fields vector is equal in length to the pCovs list length we will allow the assignment
            ## If the pCovs list length is a multiple of the fields vector, will will replicate the fields vector such that it matches in length
            ## to the pCovs list.
            field = function(value) {
                if(!missing(value)){
                    if(!is.character(value) && !all(is.na(value))){
                       stop('Error: field must be of class character')
                    }
                    if(length(value) != length(private$pCovs) & length(private$pCovs) %% length(value) != 0){
                        stop('Error: Length of field must be of length equal to the number of Covariates or a divisor of number of covariates.')
                    }

                    if(length(private$pCovs) / length(value) != 1){
                        private$pfield = rep(value, length(private$pCovs)/length(value))
                        return(private$pfield)
                    }

                    private$pfield = value
                    return(private$pfield)

                } else{
                    return(private$pfield)
                }
            },

            ## Covariate Paddinig
            ## Here we check to make sure that all pad are of class numeric and that they are the same length as pCovs -> the internal covariate list
            ## If the pad vector is equal in length to the pCovs list length we will allow the assignment
            ## If the pCovs list length is a multiple of the pad vector, will will replicate the pad vector such that it matches in length
            ## to the pCovs list.
            pad = function(value) {
                if(!missing(value)){
                    if(!is.numeric(value) && !all(is.na(value))){
                        stop("Error: pad must be of class character")
                    }
                    if(length(value) != length(private$pCovs) & length(private$pCovs) %% length(value) != 0){
                        stop("Error: Length of pad must be of length equal to the number of Covariates or a divisor of number of covariates.")
                    }
                    if(length(private$pCovs) / length(value) != 1){
                        private$ppad = rep(value, length(private$pCovs)/length(value))
                        return(private$ppad)
                    }

                    private$ppad = value
                    return(private$ppad)

                } else{
                    return(private$ppad)
                }
            },

            ##Covariate na.rm
            ## Here we check to make sure that all na.rms are of class logical and that they are the same length as pCovs -> the internal covariate list
            ## If the na.rms vector is equal in length to the pCovs list length we will allow the assignment
            ## If the pCovs list length is a multiple of the na.rms vector, will will replicate the na.rms vector such that it matches in length
            ## to the pCovs list.
            na.rm = function(value) {
                if(!missing(value)){
                    if(!is.logical(value) && !all(is.na(value))){
                        stop("Error: na.rm must be of class character")
                    }
                    if(length(value) != length(private$pCovs) & length(private$pCovs) %% length(value) != 0){
                        stop("Error: Length of na.rm must be of length equal to the number of Covariates or a divisor of number of covariates.")
                    }

                    if(length(private$pCovs) / length(value) != 1){
                        private$pna.rm = rep(value, length(private$pCovs)/length(value))
                        return(private$pna.rm)
                    }

                    private$pna.rm = value
                    return(private$pna.rm)

                } else{
                    return(private$pna.rm)
                }
            },

            ## Covariate Grep
            ## Here we check to make sure that all greps are of class chracter and that they are the same length as pCovs -> the internal covariate list
            ## If the greps vector is equal in length to the pCovs list length we will allow the assignment
            ## If the pCovs list length is a multiple of the greps vector, will will replicate the greps vector such that it matches in length
            ## to the pCovs list.
            grep = function(value) {
                if(!missing(value)){
                    if(!is.character(value) && !all(is.na(value))){
                        stop("Error: grep must be of class character")
                    }
                    if(length(value) != length(private$pCovs) & length(private$pCovs) %% length(value) != 0){
                        stop("Error: Length of grep must be of length equal to the number of Covariates or a divisor of number of covariates.")
                    }
                    if(length(private$pCovs) / length(value) != 1){
                        private$pgrep = rep(value, length(private$pCovs)/length(value))
                        return(private$pgrep)
                    }

                    private$pgrep = value
                    return(private$pgrep)
                } else{
                    return(private$pgrep)
                }
            },

            ##Covariate Covs
            cvs = function(value) {
                if(!missing(value)){
                    private$pCovs = value
                    return(private$pCovs)
                } else{
                    return(private$pCovs)
                }
            }
    ),
)

    



#' @name c.Cov_Arr
#' @title title
#' @description
#'
#' Override the c operator for covariates so that you can merge them like a vector
#'
#' @param ... A series of Covariates, note all objects must be of type Cov_Arr 
#' @return Cov_Arr object that can be passed directly into the FishHook object constructor that contains all of the Cov_Arr covariates
#' Passed in the ... param
#' @author Zoran Z. Gajic
#' @export
'c.Cov_Arr' = function(...){

    ##Ensure that all params are of type Cov_Arr
    Cov_Arrs = list(...)
    isc = sapply(Cov_Arrs, function(x)  class(x)[1] == 'Cov_Arr')

    if(any(!isc)){
        stop('Error: All inputs must be of class Cov_Arr.')
    }

    ## Merging vars of the covariates
    names  = unlist(sapply(Cov_Arrs, function(x) x$names))
    type  = unlist(sapply(Cov_Arrs, function(x) x$type))
    signature  = unlist(sapply(Cov_Arrs, function(x) x$signature))
    field  = unlist(sapply(Cov_Arrs, function(x) x$field))
    pad  = unlist(sapply(Cov_Arrs, function(x) x$pad))
    na.rm  = unlist(sapply(Cov_Arrs, function(x) x$na.rm))
    grep  = unlist(sapply(Cov_Arrs, function(x) x$grep))

    ## Merging Covariates
    covs = lapply(Cov_Arrs, function(x) x$cvs)
    Covs = unlist(covs, recursive = F)

    ##Creating a new Cov_Arr and assigning all of the merged variables to it
    ret = Cov_Arr$new()
    ret$cvs = Covs
    ret$names = names
    ret$type = type
    ret$signature = signature
    ret$field = field
    ret$pad = pad
    ret$na.rm = na.rm
    ret$grep = grep

    return(ret)
}




#' @name [.Cov_Arr
#' @title title
#' @description
#'
#' Overrides the subset operator x[] for use with Cov_Arr to allow for vector like subsetting
#'
#' @param obj Cov_Arr This is the Cov_Arr to be subset
#' @param range vector This is the range of Covariates to return, like subsetting a vector. e.g. c(1,2,3,4,5)[3:4] == c(3,4)
#' @return A new Cov_Arr object that contains only the Covs within the given range
#' @author Zoran Z. Gajic
#' @export
'[.Cov_Arr' = function(obj, range){
    ##Clone the object so we don't mess with the original
    ret = obj$clone()
    ##Call the subset function of the Cov_Arr class that will modify the cloned Cov_Arr
    ret$subset(range)
    return (ret)
}






#' @name FishHook
#' @title title
#' @description
#'
#' Stores Events, Targets, Eligible, Covariates.
#'
#' @param targets Examples of targets are genes, enhancers, or even 1kb tiles of the genome that we can then convert into a rolling/tiled window. This param must be of class "GRanges".
#' @param events Events are the given mutational regions and must be of class "GRanges". Examples of events are SNVs (e.g. C->G) somatic copy number alterations (SCNAs), fusion events, etc.
#' @param eligible Eligible regions are the regions of the genome that have enough statistical power to score. For example, in the case of exome sequencing where all regions are not equally
#' represented, eligible can be a set of regions that meet an arbitrary exome coverage threshold. Another example of when to use eligibility is in the case of whole genomes,
#' where your targets are 1kb tiles. Regions of the genome you would want to exclude in this case are highly repetative regions such as centromeres, telomeres, and satelite repeates.
#' This param must be of class "GRanges".
#' @param covariates Covariates are genomic covariates that you belive will cause your given type of event (mutations, CNVs, fusions, case control samples) that are not linked to the process you are
#' investigating (e.g. cancer drivers). In the case of cancer drivers, we are looking for regions that are mutated as part of cancer progression. As such, regions that are more suceptable to
#' random mutagenesis such as late replicating or non-expressed region (transcription coupled repair) could become false positives. Including covariates for these biological processes will
#' reduce thier visible effect in the final data. This param must be of type "Cov_Arr".
#' @param out.path A character that will indicate a system path in which to save the results of the analysis.
#' @param use_local_mut_density A logical that when true, creates a covariate that will represent the mutational density in the genome, whose bin size will be determined by local_mut_density_bin.
#' This covariate can be used when you have no other covariates as a way to correct for variations in mutational rates along the genome under the assumption that driving mutations
#' will cluster in local regions as opposed to global regions. This is similar to saying, in the town of foo, there is a crime rate of X that we will assume to be the local crime rate
#' If a region in foo have a crime rate Y such that Y >>>>> X, we can say that region Y has a higher crime rate than we would expect.
#' @param local_mut_density_bin A numeric value that will indicate the size of the genomic bins to use if use_local_mut_density = TRUE. Note that this paramter should be a few orders of
#' magnitude greater than the size of your targets. e.g. if your targets are 1e5 bps long, you may want a local_mut_density_bin of 1e7 or even 1e8
#' @param gennome A character value indicating which build of the human genome to use, by default set to hg19
#' @param mc.cores A numeric value that indicates the amount of computing cores to use when running fishHook. This will mainly be used during the annotation step of the analysis, or during
#' initial instantiation of the object if use_local_mut_density = T
#' @param na.rm A logical indicating how you handle NAs in your data, mainly used in fftab and gr.val, see these function documentations for more information
#' @param pad A numeric indicating how far each covariate range should be extended, see Cov_Arr for more information, not that this will only be used if atleast on of the
#' Covariates have pad = NA
#' @param vebose A logical indicating whether or not to print information to the console when running FishHook
#' @param max.slice integer Max slice of intervals to evaluate with  gr.val (default = 1e3)
#' @param ff.chunk integer Max chunk to evaluate with fftab (default = 1e6)
#' @param max.chunk integer gr.findoverlaps parameter (default = 1e11)
#' @param ptidcol A character, that indicates the column name containing the patient ids, this is for use in conjunction with maxpatientpergene. If max patientpergene is specified and
#' and the column referenced by ptidcol exists, we will limit the contributions of each patient to each target to maxpatientpergene. e.g. if Patient A has 3 events in target A and Patient B
#' has 1 event in target A, and maxpatientpergene is set to 2, with thier ID column specified, target A will have a cournt of 3, 2 coming from patient A and 1 coming from patient B
#' @param maxpatientpergene a numeric that indicates the max number of events any given patient can contribute to a given target. for use in conjction with ptidcol. see ptidcol for more info.
#' @param weightEvents a logical that indicates if the events should be weighted by thier overlap with the targets. e.g. if we have a SCNA spanning 0:1000 and a target spanning 500:10000, the overlap
#' of the SCNA and target is 500:1000 which is half of the original width of the SCNA event. thus if weightEvent = T, we will credit a count of 0.5 to the target for this SCNA. This deviates from
#' the expected input for the gamma poisson as the gamma poisson measures whole event counts.
#' @param nb boolean negative binomial, if false then use poisson
#' @return FishHook object ready for annotation/scoring.
#' @author Zoran Z. Gajic
#' @importFrom R6 R6Class
#' @export
FishHook = R6::R6Class('FishHook',

    public = list(

        ##See class documentation for params
        initialize = function(targets = NULL, out.path = NULL, eligible = NULL, ... ,events = NULL, covariates = NULL,
            use_local_mut_density = FALSE, local_mut_density_bin = 1e6, genome = 'BSgenome.Hsapiens.UCSC.hg19::Hsapiens',
            mc.cores = 1, na.rm = TRUE, pad = 0, verbose = TRUE, max.slice = 1e3, ff.chunk = 1e6, max.chunk = 1e11, ptidcol = NULL,
            maxpatientpergene = Inf, weightEvents = FALSE, nb = TRUE){
             ## This next portion checks to make sure that the seqlevels are in the same format
            if(!is.null(covariates)){
                ## Gets whther seqlevels of covariates are chr or not chr
                seqLevelsStatus_Covariates = covariates$chr()
                ## Warns if there is a heterogenetiy of seqlevels (chr or not)
                if(length(unique(seqLevelsStatus_Covariates)) > 1){
                    warning('Warning:Covariates appears to have mismatched seqlevels, make sure all Covariates have seqlevels that start with chr or do not', call.=TRUE)
                }
            }

            ## gets the seqlevels and looks for chr to indicate USCS format
            seqLevelsStatus_Targets = any(grepl('chr', GenomeInfoDb::seqlevels(targets)))
            seqLevelsStatus_Events = any(grepl('chr', GenomeInfoDb::seqlevels(events)))

            if(!is.null(covariates)){
                if(any(!(seqLevelsStatus_Targets %in% seqLevelsStatus_Covariates))){
                    warning('Warning: seqlevels of Targets and Covariates appear to be in different formats')
                }
            }

            if(!is.null(eligible)){
                seqLevelsStatus_Eligible = any(grepl("chr", GenomeInfoDb::seqlevels(eligible)))
                if(seqLevelsStatus_Targets != seqLevelsStatus_Eligible){
                    warning('Warning: seqlevels of Targets and Eligible appear to be in different formats')
                }
            }
            if(seqLevelsStatus_Targets != seqLevelsStatus_Events){
                warning('Warning: seqlevels of Targets and Events appear to be in different formats')
            }

            ## This next portion checks to make sure there is atleast some overlap of seqlevels i.e. some mapability
            if(!any(GenomeInfoDb::seqlevels(targets) %in% GenomeInfoDb::seqlevels(events))){
                stop('Error: there are no seqlevels of events that match targets')
            }

            if(!is.null(eligible)){
                if(!any(GenomeInfoDb::seqlevels(targets) %in% GenomeInfoDb::seqlevels(eligible))){
                    stop('Error: there are no seqlevels of eligible that match targets')
                }
            }

            if(!is.null(covariates)){
                if(any(!(unlist(lapply(covariates$seqlevels(),function(x) any(x %in% GenomeInfoDb::seqlevels(targets))))))){
                    warning("Warning: atleast one of the covariates has no seqlevels in common with targets")
                }
            }

            ## Initializes and Validates targets
            self$targets = targets

            ## Initializes and Validates out.path
            self$out.path = out.path

            ## Initializes and Validates covariates
            if(is.null(covariates)){
                covariates = Cov_Arr$new()
            }

            if(class(covariates)[1] != "Cov_Arr"){
                stop("Error: Covariates must be a vector of covariates or an object of class 'Cov_Arr'")
            }

            self$cvs = covariates

            ## Initializes and Validates events
            self$events = events

            ## Initializes and Validates eligible
            if(!(is.null(eligible))){
                self$eligible = eligible
            }

            ##Creating the local mutational denisty track
            if(use_local_mut_density){
                Sys.setenv(DEFAULT_BSGENOME = genome)
                bins = gr.tile(hg_seqlengths(), local_mut_density_bin)
                f1 = FishHook$new(targets = bins, events = events, eligible = eligible)
                f1$annotate(mc.cores = mc.cores, na.rm = na.rm, verbose = verbose, max.slice = max.slice, ff.chunk = ff.chunk, max.chunk = max.chunk)
                f1$score()
                local_mut_density = seg2gr(f1$scores)[,'count.density']
                cd = local_mut_density$count.density
                avg_cd = mean(cd, na.rm = T)
                cd[is.na(cd) | cd == Inf] = avg_cd
                local_mut_density$count.density = cd
                if(length(private$pcovariates$toList()) == 0 ){
                    private$pcovariates = c(Cov_Arr$new(cvs = c(local_mut_density), type = c('numeric'), name = c("Local Mutation Density"), field = c("count.density")))
                } else{
                    private$pcovariates = c(Cov_Arr$new(cvs = c(local_mut_density), type = c('numeric'), name = c("Local Mutation Density"), field = c("count.density")), private$pcovariates)
                }
            }


            private$pmc.cores = mc.cores
            private$pna.rm = na.rm
            private$ppad = pad
            private$pverbose = verbose
            private$pmax.slice = max.slice
            private$pff.chunk = ff.chunk
            private$pmax.chunk = max.chunk
            private$pptidcol = ptidcol
            private$pmaxpatientpergene = maxpatientpergene
            private$pweightEvents = weightEvents
            private$pnb = nb
        },

        ## Params:
        ## x Cov calls toList on X
        ## Return:
        ## x$toList
        ## UI:
        ## None
        ## Notes:
        ## This function only exists for legacy purposes
        toList = function(x){
            return (x$toList())
        },


        ## Params:
        ## no params, any passed params will be ignored
        ## Return:
        ## none
        ## UI:
        ## prints a summary of the internal state of the FishHook object
        print = function(){
            targ = paste('Contains' , length(private$ptargets), "hypotheses." ,collapse = "")
            eve = paste('Contains', length(private$pevents), "events to map to hypotheses.", collapse = "")
            if(is.null(private$peligible)){
                elig = "All regions are elgible."
            } else{
                elig = "Will map only eliglble regions."
            }
            if(is.null(private$pcovariates$names)){
                covs = "No covariates will be used."
            } else{
                cov.names = private$pcovariates$names
                covs = cov.names
            }
            meta = paste('Targets contains', ncol(values(private$ptargets)), 'metadata columns')
            state = paste('Current State:', private$pstate)
            cat(targ, eve, elig, 'Covariates:', covs, meta, state, sep = '\n', collapse = '\n')
        },

        ## Params:
        ## mc.cores, see FishHook class documentation for more info
        ## na.rm, see FishHook class documentation for more info
        ## pad, see FishHook class documentation for more info
        ## verbose, see FishHook class documentation for more info
        ## max.slice, see FishHook class documentation for more info
        ## ff.chunk, see FishHook class documentation for more info
        ## max.chunk, see FishHook class documentation for more info
        ## ptidcol, see FishHook class documentation for more info
        ## maxpatientpergene, see FishHook class documentation for more info
        ## Return:
        ## None
        ## UI:
        ## If verbose = T, will print updates as the annotation proceeds
        ## Notes:
        ## This function changes the internal state of the fishHook object and sets the state to 'Annotated'
        annotate = function(mc.cores = private$pmc.cores, na.rm = private$pna.rm, pad = private$ppad,
            verbose = private$pverbose, max.slice = private$pmax.slice, ff.chunk = private$pff.chunk,
            max.chunk = private$pmax.chunk, ptidcol = private$pptidcol, maxpatientpergene = private$maxpatientpergene,
            weightEvents = private$pweightEvents){

            if(private$pstate == "Scored"){
                stop("Error: You have a scored object already created. If you want to re-run the analysis you can clear the scored and annotated objects using fish$clear()")
            }

            private$panno = annotate.targets(targets = private$ptargets,
                covered = private$peligible,
                events = private$pevents,
                mc.cores = mc.cores,
                na.rm = na.rm,
                pad = pad,
                verbose = verbose,
                max.slice = max.slice,
                ff.chunk = ff.chunk,
                max.chunk = max.chunk,
                out.path = private$pout.path,
                covariates = private$pcovariates$toList(),
                ptidcol = ptidcol,
                maxpatientpergene = maxpatientpergene,
                weightEvents = weightEvents)

            private$pstate = "Annotated"

        },



        ## Params:
        ## targets, a GRanges that is the output of annotate.targets. note that this is for admin degbugging and
        ## should never be touched by the user unless you absoltely know exactly what you are doing and why you are doing it.
        ## by, character vector with which to split into meta-territories (default = NULL)
        ## fields, a character vector indicating which columns to be used in aggregateion by default all meta data
        ## fields of targets EXCEPT reserved field names $coverage, $counts, $query.id (default = NULL)
        ## rolling, positive numeric (integer) specifying how many (genome coordinate) adjacent to aggregate in a rolling
        ## fashion; positive integer with which to performa rolling sum / weighted average WITHIN chromosomes of "rolling" ranges" --> return a granges (default = NULL)
        ## For example, if we cut a chromosome into 5 pieces (1,2,3,4,5) and set rolling = 3, we will get an aggregated dataset (123,234,345) as the internal value
        ## This is mainly for use with whole genome analysis in order to speed up the annotation step
        ## disjoint, boolean only take disjoint bins of input (default = TRUE)
        ## na.rm,  boolean only applicable for sample wise aggregation (i.e. if by = NULL) (default = FALSE)
        ## FUN, list only applies (for now) if by = NULL, this is a named list of functions, where each item named "nm" corresponds to an
        ## optional function of how to alternatively aggregate field "nm" per samples, for alternative aggregation of coverage and count.
        ## This function is applied at every iteration of loading a new sample and adding to the existing set.   It is normally sum [for coverage and count]
        ## and coverage weighted mean [for all other covariates].  Alternative coverage / count aggregation functions should have two
        ## arguments (val1, val2) and all other alt covariate aggregation functions should have four arguments (val1, cov1, val2, cov2)
        ## where val1 is the accumulating vector and val2 is the new vector of values. 
        ## verbose, boolean verbose flag (default = TRUE)
        ## Return:
        ## None
        ## UI:
        ## If verbose = T, will print updates as the aggregation proceeds
        ## Notes:
        ## This function changes the internal state of the fishHook object and sets the state to 'Aggregated'
        aggregate = function(targets = private$panno, by = NULL, fields = NULL, rolling = NULL, disjoint = TRUE, na.rm = FALSE, FUN = list(), verbose = private$pverbose){

            if(private$pstate == "Initialized"){
                self$annotate()
            }

            agg = aggregate.targets(targets,
                by = by,
                fields = fields,
                rolling = rolling,
                disjoint = disjoint,
                na.rm = na.rm,
                FUN = list(),
                verbose = TRUE
            )

            private$paggregated = agg

            dump = self$clear("Aggregated")

        },



        ## Params:        
        ## targets,  a GRanges that is the output of annotate.targets. note that this is for admin degbugging and
        ## should never be touched by the user unless you absoltely know exactly what you are doing and why you are doing it.
        ## annotated targets with fields $coverage, optional field, $count and additional numeric covariates
        ## model, boolean if true,  fit existing model --> covariates must be present (default = NULL)
        ## nb, boolean If TRUE, uses negative binomial; if FALSE then use Poisson
        ## verbose, boolean verbose flag (default = TRUE)
        ## iter, integer info (default = 200)
        ## subsample, interger info (default = 1e5)
        ## seed, numeric (integer) indicated the random number seed to be used.  (default = 42)
        ## p.randomized, boolean Flag info (default = TRUE)
        ## Return:
        ## None
        ## UI:
        ## If verbose = T, will print updates as the scoring proceeds
        ## Notes:
        ## This function changes the internal state of the fishHook object and sets the state to 'Scored'
        score = function(nb = private$pnb, verbose = private$pverbose,  iter = 200, subsample = 1e5, seed = 42, p.randomize = TRUE){

            if(private$pstate == "Initialized"){
                self$annotate()
            }

            ## If we are aggregated we should score that, if we are not we should score anno
            if(private$pstate == "Aggregated"){
                targ = private$paggregated
                covs = names(values(private$paggregated[[1]]))
            } else{
                targ = private$panno
                covs = names(values(private$panno))
            }

                                        #pls

            #print(targ)

            ## Scoring
            score = score.targets(targ,
                covariates = covs,
                return.model = TRUE,
                nb = nb,
                verbose = verbose,
                iter = iter,
                subsample = subsample,
                seed = seed,
                classReturn = TRUE,
                p.randomize = p.randomize)

            private$pscore = score[[1]]
            private$pmodel = score[[2]]
            private$pstate = 'Scored'
        },

        ## Params:
        ## state, a character indicating which state to revert to, e.g. if you are 'Scored' you can revert to 'Annotated', 'Initialized', an possibly 'Aggregated'
        ## however if your state is 'Initialized' you cannot revert to a 'Scored' state
        ## Return:
        ## none
        ## UI:
        ## None
        ## Notes: !!==WARNING==!! This function will delete any data pertaining to steps that come after the reversion state, so if you revert a scored object to 'Initialized'
        ## You will lose the scored data and the annotated data as well as any aggregate data
        clear = function(state = 'Initialized'){
            if(state == 'Initialized'){
                private$pstate = 'Initialized'
                private$pmodel = NULL
                private$pscore = NULL
                private$panno = NULL
                private$paggregated = NULL
                return('Clear Completed')
            }
            if(state == 'Annotated'){
                private$pstate = 'Annotated'
                private$pmodel = NULL
                private$pscore = NULL
                private$paggregated = NULL
                return('Clear Completed')
            }
            if(state == 'Aggregated'){
                private$pstate = 'Aggregated'
                private$pmodel = NULL
                private$pscore = NULL
                return('Clear Completed')
            }

            return('Valid reversion state not specified. This is not a major error, just letting you know that nothing has been chaged')
        },


        ## Params:
        ## plotly, boolean
        ## columns, character vector, that indicates the names of the columns from the fishHook$all output to use in ploting.
        ## Note that this is only used if plotly = T
        ## annotations, a named list of character vectors. Each vector must have the same number of rows  as the fishHook$score datatable
        ## table. Each vector will be used to annotate the plot, only if plotly = T
        ## key, a character that is passed to the plotly function that will link each point to a give value. For example, if key is set to gene_name
        ## The ploted points are refered to by thier gene_name. This is useful when integrating with shiny or any other tool that can
        ## integrate with plotly plots.
        ## Return:
        ## plotly object that can be plotted
        ## UI:
        ## None
        qq_plot = function(plotly = TRUE, columns = NULL, annotations = NULL, key = NULL, ...){
            res = self$all

            if(!is.null(columns)){
                columns = columns[columns %in% names(res)]
                annotation_columns = lapply(columns, function(x) as.character(unlist(res[,x,with=FALSE])))
                names(annotation_columns) = columns
            } else{
                annotation_columns = list()
            }

            if(is.null(annotations) & is.null(columns)){

                if(is.null(res$name)){
                    names = c(1:nrow(res))
                } else{
                    names = res$name
                }

                annotations = list(Hypothesis_ID = names,
                    Count = res$count,
                    Effectsize = round(res$effectsize,2),
                    q = res$q)

            } else{
                res = self$all
            }

            return(qq_pval(res$p ,annotations = c(annotations,annotation_columns), gradient = list(Count = res$count), titleText = "" ,  plotly = plotly, key = key))

        }
    ),

    
    ## Private variables are internal variables that cannot be accessed by the user
    ## These variables will have active representations that the user can interact with the update
    ## and view these variables, all internal manipulations will be done with these private variables
    private = list(
        ## Genomic Ranges Object that Indicates the Hypotheses
        ptargets = NULL,

        ## Eligible Regions, this could be things such as regions of the genome capture with
        ## whole exome sequencing or in the case of whole genome sequencing non-repetative regions
        peligible = NULL,

        ## Events to Count, these can be things such as wes single nucleotide variants, microarracy somatic copy number variations, fusions, breakpoints, etc.
        pevents = NULL,

        ## Covariates list for passing to fishHook
        pcovariates = NULL,

        ## Valid Covariate Types
        pCOV.TYPES = c('numeric', 'sequence', 'interval'),

        ## Valid Covariate Classes
        pCOV.CLASSES = c('GRanges', 'RleList', 'ffTrack', 'character'),

        ## Path to write the annotated data to
        ## Useful if working with long proccessing times due to
        ## Many Covariates
        pout.path = NULL,

        ##The number of cores to use for the analysis
        pmc.cores = 1,

        ##The internal state of the object
        pstate = "Initialized",

        ##na.rm, see fishHook$annotate/fishHook$aggregate/fishHook$score for more info
        pna.rm = TRUE,

        ##padding to use with the events, see annotate.targets for more info
        ppad = 0,

        ##global verbose paramter
        pverbose = TRUE,

        ##see annotate.targets for more info
        pmax.slice = 1e3,

        ##see annotate.targets for more info
        pff.chunk = 1e6,

        ##see annotate.targets for more info
        pmax.chunk = 1e11,

        ##see annotate.targets for more info
        pptidcol = NULL,

        ##see annotate.targets for more info
        pmaxpatientpergene = Inf,

        ##see annotate.targets for more info
        pweightEvents = FALSE,

        ##The variable containing the output of fishHook$annotate()
        panno = NULL,

        ##The variable containing the output of fishHook$score()
        pscore = NULL,

        ##see score.targets for more info
        pmodel = NULL,

        ##see score.targets for more info
        preturn.model = TRUE,

        ##see score.targets for more info
        pnb = TRUE,

        ##The variable containing the output of fishHook$aggregate()
        paggregated = NULL

        ),

    ## The active list contains a variable for each private variable.
    ## Active variables are for user interaction,
    ## Interactions can be as such
    ## class$active will call the active variable function with the value missing
    ## class$active = value will call the active variable function with the value = value
    active = list(

        ## Covariates = cvs
        ## Here we check to make sure that all cvs  are of class Cov_Arr
        ## We then reset the object to its initialized state so as to not introduce incosistencies amongst variables
        cvs = function(value) {
            if(!missing(value)){
                if(!(class(value)[1] == 'Cov_Arr')  & !is.null(value)){
                    stop('Error: covariates must be of class Cov_Arr')
                }

                self$clear()

                private$pcovariates = value

                return(private$pcovariate)

            } else{
                return(private$pcovariates)
            }
        },

        ## Eligible
        ## Here we check to make sure that eligible is of class GRanges
        ## We then reset the object to its initialized state so as to not introduce incosistencies amongst variables
        eligible = function(value) {
            if(!missing(value)){
                if((!class(value) == 'GRanges') & !is.null(value)){
                    stop('Error: eligible must be of class GRanges')
                }

                self$clear()

                private$peligible = value
                return(private$peligible)

            } else{
                return(private$peligible)
            }
        },

        ## Targets
        ## Here we check to make sure that targets is of class GRanges or a chracter path and are not NULL
        ## We then reset the object to its initialized state so as to not introduce incosistencies amongst variables
        targets = function(value) {

            if(!missing(value)){
                if(!(class(value) == 'GRanges') && !(class(value) == 'character')){
                    stop('Error: targets must be of class GRanges')
                }

            targets = value
            ## checks if targets is NULL
            if (is.null(targets)){
                stop('Targets cannot be NULL.')
            }

            ## checks to see if targets is a path & import if so
            if (is.character(targets)){
                if (grepl('\\.rds$', targets[1])){
                    targets = readRDS(targets[1])
                } else if (grepl('(\\.bed$)', targets[1])){
                    targets = import.ucsc(targets[1])
                }
            }

            ## Forces targets to be a GRanges Objects
            if(class(targets) != 'GRanges'){
                stop('Error: Loaded or provided class of targets must be "GRanges"')
            }

            ## checks to see if target contains any data
            if (length(targets)==0){
                stop('Error: Must provide non-empty targets')
            }

            ## Looks for a "name" field to index/Identify targets by name
            ## If no such field is found creates a set of indexes
            if(is.null(targets$name)){
                targets$name = 1:length(targets)
            }

            ## Change here when making the smart swaps
            self$clear()

            private$ptargets = targets

            return(private$ptargets)

            } else{
                return(private$ptargets)
            }
        },

        ## Events
        ## Here we check to make sure that events is of class GRanges is not NULL
        ## We then reset the object to its initialized state so as to not introduce incosistencies amongst variables
        events = function(value) {
            if(!missing(value)){
                if(!(class(value) == 'GRanges')){
                    stop('Error: Events must be of class GRanges')
               }

                events = value

               ## Forces Events to Exist
               if(is.null(events)){
                   stop('Error: Events must exist and cannot be NULL')
               }

               ## Forces Events to be a GRanges  Object
               if (class(events) != 'GRanges'){
                   stop('Error: Events must be of class "GRanges"')
               }

               private$pevents = events

               ## Change here when making the smart swaps
               self$clear()

               return(private$pevents)

            } else{
                return(private$pevents)
            }
        },
        
        ## out.path
        ## Here we check to make sure that out.path is of class character and that it exists
        out.path = function(value) {
            if(!missing(value)){
                if(!(class(value) == 'character')  && !is.null(value)){
                    stop('Error: out.path must be of class character')
                }

                out.path = value
                ## checks to see if out.path is able to be written to
                ## throws a warning if unable to write
                if (!is.null(out.path)){
                    tryCatch(saveRDS(private$ptargets, paste(gsub('.rds', '', out.path), '.source.rds', sep = '')),
                        error = function(e) warning(sprintf('Error: writing to file %s', out.path)))
                }

                private$pout.path = out.path

                return(private$pout.path)

            } else{
                return(private$pout.path)
            }
        },

        ## anno
        ## Here we check to make sure that anno is of class GRanges
        anno = function(value) {
            if(!missing(value)){
                if(!(class(value) == 'GRanges')  && !is.null(value)){
                    stop('Error: anno must be of class GRanges')
                } else{
                    warning('Warning: You are editing the annotated dataset generated by fishHook, if you are trying to change targets use fish$targets.')
                }

                private$panno = value

                return(private$panno)

           } else{
                return(private$panno)
           }
        },

        ## scores
        ## Here we check to make sure that scores is of class data.table
        scores = function(value) {
            if(!missing(value)){
                if(!(class(value) == 'data.table')  && !is.null(value)){
                    stop('Error: score must be of class data.table')
                } else{
                    warning('Warning: You are editing the annotated dataset generated by fishHook, if you are trying to change targets use fish$targets.')
                }

                private$pscore = value

                return(private$pscore)

            } else{
                return(private$pscore)
            }
        },

        ## model
        ## !!==WARNING==!! Do not edit this variable unless you really know what you're doing
        model = function(value) {
            if(!missing(value)){

                warning('Warning: You are editing the regression model generated by fishHook. Unless you know what you arere doing I would recomend reverting to a safe state using fish$clear()')

                private$pmodel = value

                return(private$pmodel)

            } else{
                return(private$pmodel)
            }
        },

        ## mc.cores
        ## Here we check to make sure that scores of of class numeric and is a positive value, it is then floored for safety
        mc.cores = function(value) {
            if(!missing(value)){
                if(!(class(value) == 'numeric')  && !is.null(value) && value > 0){
                    stop('Error: mc.cores must be of class numeric')
                }

                private$pmc.cores = floor(value)

                return(private$pmc.cores)

            } else{
                return(private$pmc.cores)
            }
        },

        ## na.rm
        ## boolean
        na.rm = function(value) {
            if(!missing(value)){
                if(!(class(value) == 'logical')  && !is.null(value)){
                    stop('Error: na.rm must be of class logical')
                }

                private$pna.rm = value

                return(private$pna.rm)

            } else{
                return(private$pna.rm)
            }
        },

        ## pad
        ## numeric
        pad = function(value) {
            if(!missing(value)){
                if(!(class(value) == 'numeric')  && !is.null(value)){
                    stop('Error: pad must be of class numeric')
                }

                private$ppad = value

                return(private$ppad)

            } else{
                return(private$ppad)
            }
        },

        ## verbose
        ## logical
        verbose = function(value) {
            if(!missing(value)){
                if(!(class(value) == 'logical')  && !is.null(value)){
                    stop('Error: verbose must be of class logical')
                }

                private$pverbose = value

                return(private$pverbose)

            } else{
                return(private$pverbose)
            }
        },

        ## max.slice
        ## numeric
        max.slice = function(value) {
            if(!missing(value)){
                if(!(class(value) == "numeric")  && !is.null(value)){
                    stop('Error: max.slice must be of class numeric')
                }

                private$pmax.slice = value

                return(private$pmax.slice)

            } else{
                return(private$pmax.slice)
            }
        },

        ## ff.chuck
        ## numeric
        ff.chunk = function(value) {
            if(!missing(value)){
                if(!(class(value) == "numeric")  && !is.null(value)){
                    stop('Error: ff.chunk  must be of class numeric')
                }

                private$pff.chunk = value

                return(private$pff.chunk)

            } else{
                return(private$pff.chunk)
                }
        },

        ## max.chunk
        ## numeric
        max.chunk = function(value) {
            if(!missing(value)){
                if(!(class(value) == "numeric")  && !is.null(value)){
                    stop('Error: max.chunk must be of class numeric')
                }

                private$pmax.chunk = value

                return(private$pmax.chunk)

            } else{
                return(private$pmax.chunk)
            }
        },

        ## ptidcol
        ## character
        ptidcol = function(value) {
            if(!missing(value)){
                if(!(class(value) == "character")  && !is.null(value)){
                    stop('Error: ptidcol must be of class character')
                }

                private$pptidcol = value

                return(private$pptidcol)

            } else{
                return(private$pptidcol)
            }
        },

        ## maxpatientpergene
        ## numeric, must be greater than 0, value is floored for safety
        maxpatientpergene = function(value) {
            if(!missing(value)){
                if(!(class(value) == "numeric")  && !is.null(value) && value > 0){
                    stop("Error: maxpatientpergene must be of class numeric")
                }

                private$pmaxpatientpergene = floor(value)

                return(private$pmaxpatientpergene)

            } else{
                return(private$pmaxpatientpergene)
            }
        },

        ## weightEvents
        ## logical
        weightEvents = function(value) {
            if(!missing(value)){
                if(!(class(value) == "logical")  && !is.null(value)){
                    stop("Error: weightEvents must be of class logical")
                }

                private$pweightEvents = value

                return(private$pweightEvents)

            } else{
                return(private$pweightEvents)
            }
        },

        ## nb
        ## logical
        nb = function(value) {
            if(!missing(value)){
                if(!(class(value) == "logical")  && !is.null(value)){
                    stop("Error: nb must be of class logical")
                }

                private$pnb = value

                return(private$pnb)

            } else{
                return(private$pnb)
            }
        },

        ## all
        ## cannot be used for assigning data, can only be used for accessing a data.table containing merged scores and meta data
        all = function(value) {
            if(!missing(value)){
                stop("Error: This is solely for accessing data. If you want to set data, use $targets")
            } else{
                meta = values(private$ptargets)
                scores = private$pscore
                return(as.data.table(cbind(scores,meta)))
            }
        },

        ## state
        ## for accessing the state of the fishHook object, see fishHook$clear() for more information
        ## This active variable cannot be used for assignment, if you want to change the state use fishHook$clear()
        state = function(value) {
            if(!missing(value)){
                stop("Error: Cannot change the state of the FishHook Object, if you want to rever to an earlier state use fishHook$clear('state')")
            } else{
                return(private$pstate)
            }
        },

        ## aggregated
        ## GRangesList containing aggregated targets, you probably shouldn't be messing with this unless
        ## you really know what you're doing
        aggregated = function(value) {
            if(!missing(value)){
                if(!(class(value) == "GRangesList")  && !is.null(value)){
                    stop('Error: aggregated must be of class GRangesList')
                } else{
                    warning('Warning: You are editing the aggregated dataset generated by fishHook, goodluck!')
                }

                private$paggregated = value

                return(private$paggregated)

            } else{
                return(private$paggregated)
            }
        }
    )
)




#' @name [.FishHook
#' @title title
#' @description
#'
#' Overrides the subset operator x[] for use with FishHook to allow for vector like subsetting, see fishHook demo for examples
#'
#' @param obj FishHook object This is the FishHookObject to be subset
#' @param i vector subset targets
#' @param j vector subset events
#' @param k vector subset covariats
#' @param l vector susbet eligible
#' @return A new FishHook object that contains only the data within the given ranges
#' @author Zoran Z. Gajic
#' @export
'[.FishHook' = function(obj, i, j, k, l){
    ret = obj$clone()

    ##i -> targets
    if(!missing(i)){
        ret$targets = ret$targets[i]
    }

    ##j -> events
    if(!missing(j)){
        ret$events = ret$events[j]
    }

    ##k -> cvs
    if(!missing(k)){
        ret$cvs = ret$cvs[k]
    }

    ##l -> eligible
    if(!missing(l)){
        ret$eligible = ret$eligible[l]
    }
    return(ret)
}





#' @name qq_pval
#' @title qq plot given input p values
#' @description
#'
#' Generates R or Shiny quantile-quantile (Q-Q) plot given (minimally) an observed vector of p values, plotted their -log1 )quantiles against corresponding -log10
#' quantiles of the uniform distribution.
#'
#' @param obs vector of pvalues to plot, names of obs can be intepreted as labels
#' @param highlight vector optional arg specifying indices of data points to highlight (i.e. color red) (default = c())
#' @param exp numeric vector, expected distribution. if default (NULL) will plot observed against a uniform distribution
#' Use this if you are expecting a non-uniform distribution. Must be equal in length to obs. (default = NULL)
#' @param lwd integer, optional, specifying thickness of line fit to data (default = 1)
#' @param col a vector of strings (colors) equivalent in length to obs, this is the color that will be used for plotting. This is only if plotly = T  (default = NULL)
#' @param col.bg string indicating the color of the background
#' @param pch integer dot type for scatter plot
#' @param cex integer dot size for scatter plot
#' @param conf.lines logical, optional, whether to draw 95 percent confidence interval lines around x-y line
#' @param max numeric, optional, threshold to max the input p values
#' @param max.x numeric, max value for the x axis
#' @param max.y numeric, max value for the y axis
#' @param label character vector, optional specifying which data points to label (obs vector has to be named, for this to work)
#' @param plotly toggles between creating a pdf (FALSE) or an interactive html widget (TRUE)
#' @param annotations named list of vectors containing information to present as hover text (html widget), must be in same order as obs input
#' @param gradient named list that contains one vector that color codes points based on value, must bein same order as obs input
#' @param titleText title for plotly (html) graph only
#' @param subsample numeric (positive integer), number of points to use for plotting, will be taken randomly from the set of obs -> p values
#' @param key a character that is passed to the plotly function that will link each point to a give value. For example, if key is set to gene_name
#' The ploted points are refered to by thier gene_name. This is useful when integrating with shiny or any other tool that can integrate with plotly plots.
#' @import plotly
#' @author Marcin Imielinski, Eran Hodis, Zoran Z. Gajic
#' @export
qq_pval = function(obs, highlight = c(), exp = NULL, lwd = 1, col = NULL, col.bg = 'black', pch = 18, cex = 1, conf.lines = TRUE, max = NULL, max.x = NULL,
    max.y = NULL,  label = NULL, plotly = FALSE, annotations = list(), gradient = list(), titleText = "", subsample = NA, key = NULL,  ...)
{
    if(!(plotly)){
        is.exp.null = is.null(exp)

        if (is.null(col)){
            col = rep('black', length(obs))
        }

        ix1 = !is.na(obs)

        if (!is.null(exp)){
            if (length(exp) != length(obs)){
                stop('Error: length of exp must be = length(obs)')
            } else{
                ix1 = ix1 & !is.na(exp)
            }
        }

        if (is.null(highlight)){
            highlight = rep(FALSE, length(obs))
        } else if (is.logical(highlight)){
            if (length(highlight) != length(obs)){
                stop('Error: argument "highlight" must be either logical vector of same length as obs or a vector of indices')
            }
        } else{
            highlight = 1:length(obs) %in% highlight
        }

        obs = -log10(obs[ix1])
        col = col[ix1]

        highlight = highlight[ix1]

        if (!is.null(exp)){
            exp = -log10(exp[ix1])
        }

        ix2 = !is.infinite(obs)

        if (!is.null(exp)){
            ix2 = ix2 &  !is.infinite(exp)
        }

        obs = obs[ix2]
        col = col[ix2]

        highlight = highlight[ix2]
        if (!is.null(exp)){
            exp = exp[ix2]
        }

        N = length(obs)

        ## create the null distribution
        ## (-log10 of the uniform)

        if (is.null(exp)){
            exp = -log(1:N/N,10)
        } else{
            exp = sort(exp)
        }

        if (is.null(max)){
            max = max(obs,exp) + 0.5
        }

        if (!is.null(max) & is.null(max.x)){
            max.x = max
        }

        if (!is.null(max) & is.null(max.y)){
            max.y  = max
        }

        if (is.null(max.x)){
            max.x <- max(obs,exp) + 0.5
        }

        if (is.null(max.y)){
            max.y <- max(obs,exp) + 0.5
        }

        if (is.exp.null){
            tmp.exp = rev(seq(0, 7, 0.01))
            ix = 10^(-tmp.exp)*N
            c95 <-  qbeta(0.975,ix,N-ix+1)
            c05 <-  qbeta(0.025,ix,N-ix+1)

            if (conf.lines){
                ## plot the two confidence lines
                plot(tmp.exp, -log(c95,10), ylim=c(0,max.y), xlim=c(0,max.x), type = 'l', axes = FALSE, xlab = '', ylab = '')
                par(new=T)
                plot(tmp.exp, -log(c05,10), ylim=c(0,max.y), xlim=c(0,max.x), type = 'l', axes = FALSE, xlab = '', ylab = '')
                par(new=T)

                p1 = rep(tmp.exp[1], 2)
                p2 = c(-log(c95,10)[1], -log(c05,10)[1])

                lines(x=p1, y=p2)
                x.coords <- c(tmp.exp,rev(tmp.exp))
                y.coords <- c(-log(c95,10),rev(-log(c05,10)))
                polygon(x.coords, y.coords, col='light gray', border=NA)
                par(new=T)
            }
        }

        ord = order(obs)

        colors = col
        colors[highlight] = 'red';

        dat = data.table(x = sort(exp), y = obs[ord], colors = colors[ord], pch = pch, cex = cex)
        if (!is.null(names(obs))){
            names = names(obs[ord])
            setkey(dat, names)
        }

        ## rough guide to subsmapling the lower p value part of the plot
        if (nrow(dat)>1e5){
            subsample = 5e4/nrow(dat)
        }

        if (is.na(subsample[1])){
            dat[, plot(x, y, xlab = expression(Expected -log[10](italic(P))), ylab = expression(Observed -log[10](italic(P))), xlim = c(0, max.x), col = colors, ylim = c(0, max.y), pch=pch, cex=cex, bg=col.bg, ...)]
        } else{
            subsample = pmin(pmax(0, subsample[1]), 1)
            dat[ifelse(x<=2, ifelse(runif(length(x))<subsample, TRUE, FALSE), TRUE), plot(x, y, xlab = expression(Expected -log[10](italic(P))), ylab = expression(Observed -log[10](italic(P))), xlim = c(0, max.y), col = colors, ylim = c(0, max.y), pch=pch, cex=cex, bg=col.bg, ...)]
        }

        if (!is.null(label)){
            if (length(label)>0){
                if (is.null(key(dat))){
                    warning('Warning: Need to provide names to input vector to draw labels')
                } else{
                    dat[list(label), text(x, y, labels=label, pos=3)];
                }
            }
        }

        lines(x=c(0, max(max.y, max.x)), y = c(0, max(max.x, max.y)), col = 'black', lwd = lwd)

        if (!is.na(subsample)){
            dat = dat[sample(nrow(dat), subsample*nrow(dat)), ]
        }

        lambda = lm(y ~ x-1, dat)$coefficients;

        lines(x=c(0, max.x), y = c(0, lambda*max.y), col = 'red', lty = 2, lwd = lwd);
        legend('bottomright', sprintf('lambda=\n %.2f', lambda), text.col = 'red', bty = 'n')
    } else{

        if(length(annotations) < 1){
            hover = do.call(cbind.data.frame, list(p = obs))
        } else{
            hover = do.call(cbind.data.frame, list(annotations, p = obs))
        }

        hover = as.data.table(hover)

        hover$key = key

        is.exp.null = is.null(exp)
        if (is.null(col)){
            col = 'black'
        }
        ix1 = !is.na(hover$p)
        if (!is.null(exp)){
            if (length(exp) != length(hover$p)){
                stop('Error: length of exp must be = length(hover$obs)')
            } else{
                ix1 = ix1 & !is.na(exp)
            }
        }
        if (is.null(highlight)){
            highlight = rep(FALSE, length(hover$p))
        } else if (is.logical(highlight)) {
            if (length(highlight) != length(hover$p)){
                stop('Error: argument "highlight" must be either logical vector of same length as obs or a vector of indices')
            }
        } else{
            highlight = 1:length(hover$p) %in% highlight
        }
        hover$obs = -log10(hover$p[ix1])
        hover = hover[ix1]
        highlight = highlight[ix1]
        if (!is.null(exp)){
            exp = -log10(exp[ix1])
        }
        ix2 = !is.infinite(hover$obs)
        if (!is.null(exp)){
            ix2 = ix2 & !is.infinite(exp)
        }
        hover = hover[ix2]
        highlight = highlight[ix2]
        if (!is.null(exp)){
            exp = exp[ix2]
        }
        N <- length(hover$obs)
        if (is.null(exp)){
            exp = -log(1:N/N, 10)
        } else{
            exp = sort(exp)
        }
        if (is.null(max)){
            max = max(hover$obs, exp) + 0.5
        } else{
            max = max
        }
        if (is.exp.null) {
            tmp.exp = rev(seq(0, 7, 0.01))
            ix = 10^(-tmp.exp) * N
            c95 = qbeta(0.975, ix, N - ix + 1)
            c05 = qbeta(0.025, ix, N - ix + 1)
            if (FALSE){
                ## Don't need if not using conf.line (might put this in the future)
                plot(tmp.exp, -log(c95, 10), ylim = c(0, max), xlim = c(0, max),
                type = "l", axes = FALSE, xlab = '', ylab = '')

                par(new = T)
                plot(tmp.exp, -log(c05, 10), ylim = c(0, max), xlim = c(0, max),
                type = "l", axes = FALSE, xlab = '', ylab = '')

                par(new = T)
                p1 = rep(tmp.exp[1], 2)
                p2 = c(-log(c95, 10)[1], -log(c05, 10)[1])
                lines(x = p1, y = p2)
                x.coords = c(tmp.exp, rev(tmp.exp))
                y.coords = c(-log(c95, 10), rev(-log(c05, 10)))
                polygon(x.coords, y.coords, col = 'light gray', border = NA)
                par(new = T)
            }
        }

        ## creating the ploting data.table (dat) and organizing the annotations to create hover text
        ord = order(hover$obs)
        hover = hover[ord]
        dat = hover
        hover$obs = NULL

        ## Creating the hover text
        if(length(colnames(hover)) > 1){
            annotation_names  = sapply(colnames(hover), paste0, ' : ')
            annotation_names_wLineBreak  = paste('<br>', annotation_names[2:length(annotation_names)],
            sep = '')
            annotation_names = c(annotation_names[1], annotation_names_wLineBreak)
        } else{
            annotation_names  = sapply(colnames(hover), paste0, ' : ')
        }

        ## Checking if there is a gradient and if so adding it to the plotting data.table (dat)
        gradient_control = FALSE
        if(length(gradient )!= 0){
            dat$grad = gradient[[1]][ord]
            gradient_control = TRUE
        } else {
            dat$grad = c()
        }

        dat$x = sort(exp)
        dat$y = dat$obs

        ## declare so we can use in If statement
        p = NULL

        ## hacky subsampling but works really well, just maxing out the number of points at 8k
        ## and removing the extra from the non-sig
        ## (looks to be -logp of 2.6 here can make this more dynamic later )

        if (nrow(dat) <= 8000){

            dat4 = dat
            dat4$obs = NULL
            dat4$x = NULL
            dat4$y = NULL
            dat4$grad = NULL

            trans = t(dat4)
            hover_text = c()
            for (i in 1:dim(trans)[2]){
                outstr = paste(c(rbind(annotation_names, trans[,i])), sep = '', collapse = '')
                hover_text = c(hover_text,outstr)
            }

            if(gradient_control){
                p = dat[, plot_ly(data = dat, x=x, y=y, key = dat$key, hoverinfo = 'text', text = hover_text, color = grad,
                                   colors = c('blue2', 'gold'), marker = list(colorbar = list(title = names(gradient[1]), len = 1)),
                                   mode = 'markers', type = 'scatter')
                    %>% layout(xaxis = list(title = '<i>Expected -log<sub>10</sub>(P)</i>'),
                               yaxis = list(title = '<i>Observed -log<sub>10</sub>(P)</i>')) ]
            } else{
                p = dat[, plot_ly(data = dat, x=x, y=y, key = dat$key, hoverinfo = 'text', text = hover_text,
                                   mode = 'markers', type = 'scatter')
                    %>% layout(xaxis = list(title = '<i>Expected -log<sub>10</sub>(P)</i>'),
                               yaxis = list(title = '<i>Observed -log<sub>10</sub>(P)</i>')) ]
            }
        } else{

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

            test3 <<-dat2

            for (i in 1:dim(trans)[2]){
                outstr = paste(c(rbind(annotation_names, trans[,i])), sep = "", collapse = "")
                hover_text = c(hover_text,outstr)
            }


            test2 <<- dat2
            if(gradient_control){
                p = dat2[, plot_ly(data = dat2, x=x, y=y, hoverinfo = 'text', key = dat2$key,  text = hover_text, color = grad,
                                    colors = c('blue2', 'gold'), marker = list(colorbar = list(title = names(gradient[1]), len = 1, lenmode = 'fraction')),
                                    mode = 'markers', type = 'scatter')
                     %>% layout(xaxis = list(title = '<i>Expected -log<sub>10</sub>(P)</i>'),
                                yaxis = list(title = '<i>Observed -log<sub>10</sub>(P)</i>')) ]
            } else{
                p = dat2[,  plot_ly(data = dat2, x=x, y=y,hoverinfo = "text", text = hover_text, mode = 'markers', type = 'scatter')
                     %>% layout(xaxis = list(title = '<i>Expected -log<sub>10</sub>(P)</i>'),
                                yaxis = list(title = '<i>Observed -log<sub>10</sub>(P)</i>')) ]
            }

        }

        ## Calculating lambda, Note that this is using the whole data set not the subsampled one
        lambda = lm(y ~ x - 1, dat)$coefficients
        lambda_max = max*as.numeric(lambda)

        ## adding shapes (lines) + title  note that html <b></b> style is used for mods and plotting lines
        ## is done by specifying two points on the line (x0/y0 and x1/y1)
        p = layout(p,
            title = sprintf('<b>%s</b>' ,titleText),
            titlefont = list(size = 24),
            shapes = list(
                list(type = 'line', line = list(color = 'black'), x0 = 0, x1  = max, xref = 'x', y0 = 0, y1 = max, yref ='y'),
                list(type = 'line', line = list(color = 'red'), x0 = 0, x1 = max, xref = 'x', y0 = 0, y1 = lambda_max, yref = 'y')),
            annotations = list(
                x = (0.9 * max),
                y = (0.03 * max),
                text = paste('lambda =', sprintf('%.2f', signif(lambda,3)), collapse = ' '),
                font = list(color = 'red', size = 20),
                showarrow = FALSE,
                xref = 'x',
                yref = 'y'
            ),
            margin = list(t = 100),
            hovermode = 'compare')
    }
}











