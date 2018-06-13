############################################################################
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




## eligible.rda  events.rda  replication_timing_cov.rda  hypotheses.rda

#' Sample events
#'
#' An object of type 'GRanges' that contains a set of events dervided from the TCGA whole exome sequencing data.
#'
#' Metadata columns:
#' id, inidcates to which sample (patient) the mutational event belongs to. There are a total
#' of 8475 patients and 1985704 total events
#' 
#' @name events
#' @docType data
#' @keywords data
#' @format \code{GRanges}
NULL


#' Sample hypotheses
#'
#' An object of type 'GRanges' that contains 19,688 human genes
#'
#' Metadata columns:
#' gene_name, inidcates the name by which this gene is refered to as. e.g. TP53
#'
#' @name hypotheses
#' @docType data
#' @keywords data
#' @format \code{GRanges}
NULL


#' Sample replication_timing, GC-content score
#'
#' An object of type 'GRanges' that contains information regarind how long each genomic region takes to replicate.
#' This will be used as a covariate in the fishHook model
#'
#' Metadata columns:
#' score, indicates the relative rate of replication timing in this region
#' 
#' @name replication_timing
#' @docType data
#' @keywords data
#' @format \code{GRanges}
NULL


#' Sample eligible
#'
#' An object of type 'GRanges' that contains all of the eligible regions of whole exome sequencing.
#' Whole exome sequencing only sequences exonic sequences and thus most of the genome should be
#' disregarded when conducting the analysis. In addition, many exonic regions are not even captured in
#' whole exome sequencing. We define an eligible (covered) region here as a region where 80% of samples
#' have mapping reads. i.e. if we sequence 10 people and only 6 (60%)  have reads in that region then we
#' would consider that region uneigible.
#'
#' Metadata columns:
#' score, indicates the percent of samples that have reads mapping to that region.
#' 
#' @name hypotheses
#' @docType data
#' @keywords data
#' @format \code{GRanges}
NULL


#' @name annotate.hypotheses
#' @title title
#' @description
#'
#' Takes input of GRanges hypotheses, an optional set of "covered" intervals, and an indefinite list of covariates which can be R objects
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
#' @param hypotheses path to bed or rds containing genomic target regions with optional target name
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
#' @param covariates list of lists where each internal list represents a covariate, the internal list can have elements: track, type,signature,name,pad,na.rm = na.rm,field,grep. See Covariate class for descriptions of what
#' each of these elements do. Note that track is equivalent to the 'Covariate' parameter in Covariate
#' @param idcap Sets the maximum number of events a patient can contribute per target (default = Inf)
#' @param idcol string Column where patient ID is stored
#' @param weightEvetns boolean If TRUE, will weight events by their overlap with hypotheses. e.g. if 10% of an event overlaps with a target
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
#' @return GRanges of input hypotheses annotated with covariate statistics (+/- constrained to the subranges in optional argument covered)
#' @author Marcin Imielinski
annotate.hypotheses = function(hypotheses, covered = NULL, events = NULL,  mc.cores = 1, na.rm = TRUE, pad = 0, verbose = TRUE, max.slice = 1e4,
    ff.chunk = 1e6, max.chunk = 1e11, out.path = NULL, covariates = list(), idcap = Inf, idcol = NULL, weightEvents = FALSE, ...)
{
  if(weightEvents){
        idcap = NULL
    }

    if (is.character(hypotheses)){
        if (grepl('\\.rds$', hypotheses[1])){
            hypotheses = readRDS(hypotheses[1])
        } else if (grepl('(\\.bed$)', hypotheses[1])){
            require(rtracklayer)
            hypotheses = rtracklayer::import(hypotheses[1], (format = "BED"))
        }
    }

    if (length(hypotheses)==0){
        stop('Error: Must provide non-empty hypotheses')
    }

    if (!is.null(out.path)){
        tryCatch(saveRDS(hypotheses, paste(gsub('.rds', '', out.path), '.source.rds', sep = '')), error = function(e) warning(sprintf('Error writing to file %s', out.file)))
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
        fmessage('Overlapping with covered intervals')
    }

    if (!is.null(covered)){
        ov = gr.findoverlaps(hypotheses, covered, verbose = verbose, max.chunk = max.chunk, mc.cores = mc.cores)
    } else {
        ov = hypotheses[, c()]
        ov$query.id = ov$subject.id = 1:length(hypotheses)
    }

    if (verbose){
        fmessage('Finished overlapping with covered intervals')
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

                if(!is.null(idcap)){

                    if(!is.numeric(idcap)){
                        stop('Error: idcap must be of type numeric')
                    }

                    if(!("ID" %in% colnames(values(events))) & is.null(idcol)){
                        events$ID = c(1:length(events))
                    }

                    ev2 = gr.findoverlaps(events,ov, max.chunk = max.chunk, mc.cores = mc.cores)

                    if(is.null(idcol)){
                        ev2$ID = events$ID[ev2$query.id]
                    } else{
                      if (!(idcol %in% names(mcols(events))))
                        {
                          stop(paste('Column', idcol, 'not found in events'))
                        }

                      ev2$ID = mcols(events)[,idcol][ev2$query.id]
                    }

                    ev2$target.id = ov$query.id[ev2$subject.id]
                    tab = as.data.table(cbind(ev2$ID, ev2$target.id))
                    counts.unique = tab[, dummy :=1][, .(count = sum(dummy)), keyby =.(V1, V2)][, count := pmin(idcap, count)][, .(final_count = sum(count)), keyby = V2]
                }

            } else{
                ## assume it is an Rle of event counts along the genome
                counts = events
                oix = 1:length(ov)
            }

            if (verbose){
                fmessage('Computing event counts')
            }

            ov$count = 0

            if (length(oix)>0 & is.null(idcap)){
                ov$count[oix] = fftab(counts, ov[oix], chunksize = ff.chunk, na.rm = TRUE, mc.cores = mc.cores, verbose = verbose)$score
            }

            if (!is.null(out.path)){
                tryCatch(saveRDS(ov, paste(out.path, '.intermediate.rds', sep = '')), error = function(e) warning(sprintf('Error writing to file %s', out.file)))
            }

            if (verbose){
                fmessage('Finished counting events')
            }
        }
    }


    for (nm in names(covariates)){

        cov = covariates[[nm]]

        if (verbose){
            fmessage('Annotating track ', nm, '')
        }

        if (cov$type == 'sequence'){

            if (is.null(cov$grep)){
                cov$grep = FALSE
            }

            if (verbose){
                fmessage('Starting fftab for track', nm, '')
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
                fmessage('Finished fftab for track', nm, '')
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
                    cov$track = rtracklayer::import(cov$track)
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
                new.col = suppressWarnings(data.frame(val = values(gr.val(ov + cov$pad, cov$track, cov$field, mc.cores = mc.cores, verbose = verbose>1,  max.slice = max.slice, max.chunk = max.chunk, mean = TRUE, na.rm = cov$na.rm))[, cov$field]))
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
                    cov$track = rtracklayer::import(cov$track)
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

            new.col = suppressWarnings(data.frame(val = gr.val(ov + cov$pad, cov$track[, c()], mean = FALSE, weighted = TRUE,  mc.cores = mc.cores, max.slice = max.slice, max.chunk = max.chunk, na.rm = TRUE)$value/(width(ov)+2*cov$pad)))
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
        values(hypotheses) = as(as.data.frame(ovdta[list(1:length(hypotheses)), ]), 'DataFrame')

        if(!is.null(idcap)){
            hypotheses$count = 0
            hypotheses$count[as.numeric(counts.unique$V2)] = counts.unique$final_count
        }

    } else{
        hypotheses$coverage = 0
    }

    hypotheses$query.id = 1:length(hypotheses)

    ix = is.na(hypotheses$coverage)
    if (any(ix)){
        hypotheses$coverage[ix] = 0
        if(!is.null(idcap)){
            hypotheses$count[hypotheses$coverage == 0] = NA
        }
    }
    if (!is.null(out.path)){
        if (file.exists(paste(out.path, '.intermediate.rds', sep = ''))){
            system(paste('rm',  paste(out.path, '.intermediate.rds', sep = '')))  ## error catch above
        }
        tryCatch(saveRDS(hypotheses, out.path), error = function(e) warning(sprintf('Error writing to file %s', out.file)))
    }


    if (is.null(hypotheses$count))
      hypotheses$count = ifelse(hypotheses$coverage == 0, NA, ifelse(is.null(events), NA, 0))
    
    return(hypotheses)

}




#' @name aggregate.hypotheses
#' @title title
#' @description
#'
#' Gathers annotated hypotheses across a vector "by" into meta-intervals returned as a GRangesList, and returns the
#' aggregated statistics for these meta intervals by summing coverage and counts, and performing a weighted average of all other meta data fields
#' (except query.id)
#'
#' If rolling = TRUE, will return a rolling collapse of the sorted input where "rolling" specifies the number of adjacent intervals that are aggregated in a rolling manner.
#' (only makes sense for tiled target sets)
#'
#' If by = NULL and hypotheses is a vector of path names, then aggregation will be done "sample wise" on the files, ie each .rds input will be assumed to comprise the same
#' intervals in teh same order and aggregation will be computed coverage-weighted mean of covariates, a sum of coverage and counts, and (if present) a Fisher combined
#' of $p values.  Covariates are inferred from the first file in the list.
#'
#' @param hypotheses annotated GRanges of hypotheses with fields $coverage, optional field, $count and additional numeric covariates, or path to .rds file of the same; path to bed or rds containing genomic target regions with optional target name
#' @param by character vector with which to split into meta-territories (default = NULL)
#' @param fields by default all meta data fields of hypotheses EXCEPT reserved field names $coverage, $counts, $query.id (default = NULL)
#' @param rolling if specified, positive integer specifying how many (genome coordinate) adjacent to aggregate in a rolling fashion; positive integer with which to performa rolling sum / weighted average WITHIN chromosomes of "rolling" ranges" --> return a granges (default = NULL)
#' @param disjoint boolean only take disjoint bins of input (default = TRUE)
#' @param na.rm boolean only applicable for sample wise aggregation (i.e. if by = NULL) (default = FALSE)
#' @param FUN list only applies (for now) if by = NULL, this is a named list of functions, where each item named "nm" corresponds to an optional function of how to alternatively aggregate field "nm" per samples, for alternative aggregation of coverage and count.  This function is applied at every iteration of loading a new sample and adding to the existing set.   It is normally sum [for coverage and count] and coverage weighted mean [for all other covariates].  Alternative coverage / count aggregation functions should have two arguments (val1, val2) and all other alt covariate aggregation functions should have four arguments (val1, cov1, val2, cov2) where val1 is the accumulating vector and val2 is the new vector of values.
#' @param verbose boolean verbose flag (default = TRUE)
#' @return GRangesList of input hypotheses annotated with new aggregate covariate statistics OR GRanges if rolling is specified
#' @author Marcin Imielinski
#' @import zoo
#' @importFrom data.table setkey := data.table as.data.table
#' @importFrom S4Vectors values values<-
#' @importFrom GenomeInfoDb seqnames
aggregate.hypotheses = function(hypotheses, by = NULL, fields = NULL, rolling = NULL, disjoint = TRUE, na.rm = FALSE, FUN = list(), verbose = TRUE)
{
    V1 = sn = st = en = keep = count = width = NULL ## NOTE fix
    if (is.null(by) & is.character(hypotheses)){
        fmessage('Applying sample wise merging')
    } else if (is.null(by) & is.null(rolling)){
        stop('Error: argument "by" must be specified and same length as hypotheses or "rolling" must be non NULL')
    }

    if (is.null(by) & is.character(hypotheses)){

        if (!all(ix <- (file.exists(hypotheses)) & grepl('\\.rds$', hypotheses))){

            warning(sprintf('Warning: %s of the  %s input files for sample wise merging either do not exist or are not .rds files. Sample wise merging (i.e. when by is null) requires .rds files of equal dimension GRanges (same intervals, same meta data column names)', sum(!ix), length(ix)))
            if (sum(ix)==0){
                stop('No files to process')
            }
            hypotheses = hypotheses[ix]
        }

        out = readRDS(hypotheses[1])
        gr = out
        if (is.null(out$coverage)){
            stop('Coverage missing for input hypotheses')
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
        out$numcases = length(hypotheses)

        if (length(hypotheses)>1)
            for (i in 1:length(hypotheses)){

                if (verbose){
                    fmessage('Processing target file', hypotheses[i], '')
                }

                if (i > 1){
                    gr = readRDS(hypotheses[i])
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
                        warning(paste(hypotheses[i], 'missing column', cf))
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
                        warning(paste(hypotheses[i], 'missing p value column, ignoring for fisher combined computation'))
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
            fields = names(values(hypotheses))
        }

        if (any(nnum <- !(sapply(setdiff(fields, 'query.id'), function(x) class(values(hypotheses)[, x])) %in% 'numeric'))){
            warning(sprintf('%s meta data fields (%s) fit were found to be non-numeric and not aggregated', sum(nnum), paste(fields[nnum], collapse = ',')))
            fields = fields[!nnum]
        }


        cfields = intersect(names(values(hypotheses)), c('coverage', 'count'))

        if (is.null(rolling)){

            by = as.character(cbind(1:length(hypotheses), by)[,2])

            if (disjoint){
                tmp.sn = paste(by, seqnames(hypotheses), sep = '_')
                tmp.dt = data.table(sn = paste(by, seqnames(hypotheses), sep = '_'), st = start(hypotheses), en = end(hypotheses), ix = 1:length(hypotheses))
                setkey(tmp.dt, sn, st)
                tmp.dt[, keep := c(TRUE, st[-1]>en[-length(st)]), by = sn]
                setkey(tmp.dt, ix)
                hypotheses = hypotheses[tmp.dt$keep, ]
                if (verbose){
                    fmessage(sprintf('Removing %s non-disjoint within group intervals, keeping %s', prettyNum(sum(!tmp.dt[, keep]), big.mark = ','), prettyNum(sum(tmp.dt[, keep]), big.mark = ',')))
                }
                by = by[tmp.dt$keep]
            }

            if (verbose){
                fmessage('Splitting into GRangesList')
            }

            out = split(hypotheses, by)

            values(out)[, 'name'] = names(out)
            values(out)[, 'numintervals'] = table(by)[names(out)]

            tadt = gr2dt(hypotheses)

            if (verbose){
                fmessage('Aggregating columns')
            }

            for (f in cfields){
                if (verbose){
                    fmessage(f, '')
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
                fmessage('Rolling using window of ', rolling, ' (output will be coordinate sorted)')
            }

            tadt = gr2dt(sort(hypotheses))

            tadt[, width := as.numeric(width)]

            tadt = tadt[seqnames %in% c(seq(22), "X")]

          if ('count' %in% cfields ) {

            if (verbose)
              {
                fmessage("Computing rolling count")
              }
                out = tadt[, list(
                    count = zoo::rollapply(count, rolling, sum, na.rm = TRUE, fill = NA),
                    start = zoo::rollapply(start, rolling, min, fill = NA),
                    end = zoo::rollapply(end, rolling, max, fill = NA),
                    coverage = zoo::rollapply(coverage, rolling, sum, fill = NA)
                ), by = seqnames]
            } else {
                out = tadt[, list(
                    start = zoo::rollapply(start, rolling, min, fill = NA),
                    end = zoo::rollapply(end, rolling, max, fill = NA),
                    coverage = zoo::rollapply(coverage, rolling, sum, fill = NA)
                ), by = seqnames]
            }
            nna.ix = !is.na(out$start)

            if (!any(nna.ix)){
                stop('Error: Malformed input, only NA ranges produced. Reduce value of running')
            }

             out = seg2gr(out[nna.ix])

            ## rolling weighted average, used below
            .rwa = function(v, w){
                zoo::rollapply(v*w, rolling, sum, na.rm = TRUE, fill = NA)/zoo::rollapply(w*as.numeric(!is.na(v)), rolling, sum, na.rm = TRUE, fill = NA)
            }

        }

        fields = setdiff(fields, c('coverage', 'count', 'query.id'))

        for (f in fields){

            if (verbose){
                fmessage(f, '')
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


#' @name score.hypotheses
#' @title title
#' @description
#'
#' Scores hypotheses based on covariates using Gamma-Poisson model with coverage as constant
#'
#' @param hypotheses annotated hypotheses with fields $coverage, optional field, $count and additional numeric covariates
#' @param covariates chracter vector, indicates which columns of hypotheses contain the covariates
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
score.hypotheses = function(hypotheses, covariates = names(values(hypotheses)), model = NULL, return.model = FALSE, nb = TRUE,
    verbose = TRUE, iter = 200, subsample = 1e5, sets = NULL, seed = 42, mc.cores = 1, p.randomized = TRUE, classReturn = FALSE)
{

    require(MASS)
    covariates = setdiff(covariates, c('count', 'coverage', 'query.id'))

    if (any(nnin = !(covariates %in% names(values(hypotheses))))){
        stop(sprintf('Error: %s covariates (%s) missing from input data', sum(nnin), paste(covariates[nnin], collapse = ',')))
    }

    if (any(nnum = !(sapply(covariates, function(x) class(values(hypotheses)[, x])) %in% c('factor', 'numeric')))){
        warning(sprintf('%s covariates (%s) fit are non-numeric or factor, removing from model', sum(nnum), paste(covariates[nnum], collapse = ',')))
        covariates = covariates[!nnum]
    }

    if (!all(c('count', 'coverage') %in% names(values(hypotheses)))){
        stop('Error: Hypotheses must have count, coverage, and query.id fields populated')
    }

    if (verbose){
        fmessage('Setting up problem')
    }

    values(hypotheses)$count = round(values(hypotheses)$count)

    if (length(unique(values(hypotheses)$count)) <= 1){
        stop('Error: "score.hypotheses" input malformed --> count does not vary!')
    }

    set.seed(seed) ## to ensure reproducibility

    if (is.null(model)){

      tdt = as.data.table(as.data.frame(values(hypotheses)[, c('count', 'coverage', covariates)]))
      tdt$coverage = ifelse(tdt$coverage ==0, NA, log(tdt$coverage))


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
                fmessage(sprintf('Subsampling ..'))
            }

            tdt = tdt[sample(1:nrow(tdt), pmin(nrow(tdt), subsample)), ]
        }

        if (verbose){
            fmessage(sprintf('Fitting model with %s data points and %s covariates', prettyNum(nrow(tdt), big.mark = ','), length(covariates)))
        }

        formula = eval(parse(text = paste('count', " ~ ", paste(c('offset(1*coverage)', covariates), collapse = "+")))) ## make the formula with covariates

        if (nb){
            g = glm.nb(formula, data = as.data.frame(tdt), maxit = iter)
        } else{
            g = glm(formula, data = as.data.frame(tdt), family = poisson)
            g$theta = Inf ## Poisson glm is glm.nb with infinite theta 
        }
    } else{
      if (verbose)
        fmessage('Scoring hypotheses against provided model without re-fitting')
        g = model
    }

    if(!(classReturn)){
        if (return.model){
            return(g)
        }
    }

    if (is(hypotheses, 'GRanges')){
        res = as.data.frame(hypotheses)
    } else{
        res = as.data.frame(values(hypotheses))
    }

    if (any(is.fact = (sapply(covariates, function(x) class(res[, x])) %in% c('factor')))){
        is.fact = (sapply(covariates, function(x) class(res[, x])) %in% c('factor'))
        ix = which(is.fact)
        new.col = lapply(ix, function(i){
            val = res[, covariates[i]]

            if (verbose){
                fmessage('Factorizing column', covariates[i], 'with', length(val), 'across', length(levels(val)), 'levels')
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
        fmessage('Computing p values for ', nrow(res), ' hypotheses.')
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

    res$fdr = signif(p.adjust(res$p, 'BH'), 2)
    if (nb){
       res$p.neg = signif(pnbinom(res$count, mu = res$count.pred, size = g$theta, lower.tail = T), 2)
    } else{
        res$p.neg = signif(ppois(res$count, lambda = res$count.pred, lower.tail = T), 2)
    }
    res$fdr.neg = signif(p.adjust(res$p.neg, 'BH'), 2)
    res$effectsize = log2(res$count / res$count.pred)

    if(!(classReturn)){
        return(as.data.table(res))
    }

  setres = NULL
  res$query.id = 1:nrow(res)
    if (!is.null(sets))
    {
      if (is.null(names(sets)))
        {
          names(sets) = paste0('set_', 1:length(sets))
        }


      names(sets) = dedup(names(sets)) ## make sure sets are all named uniquely

      setmap = data.table(names = names(sets), tmpnames = gsub('\\W', '_', paste0('set_', names(sets))))

      if (any(duplicated(setmap$tmpnames)))
        stop('Set names contain illegal special characters that cannot be resolved, please try again using only alphanumeric characters for set names')

      names(sets) = copy(setmap$tmpnames)
      setkey(setmap, tmpnames)

      if (verbose)
        fmessage('Computing p values for ', length(sets), ' hypothesis sets.')

      .score.sets = function(sid, sets, res, nb = nb)
      {
        sets = sets[sid]
        if (all(elementNROWS(sets))==0)
          return(data.table(name = names(sets), method = as.numeric(NA), p = as.numeric(NA), p.left = as.numeric(NA), p.right = as.numeric(NA), estimate = as.numeric(NA), ci.lower = as.numeric(NA), ci.upper = as.numeric(NA), effect = as.character(NA), n = 0))

        if (verbose>1)
        {
          fmessage('Scoring set(s) ', paste(names(sets), collapse = ', '))
        }

        ij = cbind(unlist(sets), rep(1:length(sets), elementNROWS(sets)))

        
        tmpres = merge(cbind(query.id = unlist(sets), set.id = rep(names(sets), elementNROWS(sets))), res)

        if (nrow(tmpres)==0)
          {
            return(data.table(name = names(sets), method = as.numeric(NA), p = as.numeric(NA), p.left = as.numeric(NA), p.right = as.numeric(NA), estimate = as.numeric(NA), ci.lower = as.numeric(NA), ci.upper = as.numeric(NA), effect = as.character(NA), n = 0))
          }

        setdata = cbind(data.frame(count = tmpres$count, count.pred = log(tmpres$count.pred)),
                        as.data.frame(as.matrix(Matrix::sparseMatrix(1:nrow(tmpres), match(tmpres$set.id, names(sets)), x = 1,
                                                             dims = c(nrow(tmpres), length(sets)),
                                                             dimnames = list(NULL, names(sets))))))
        
        
        ## make the formula with covariates
        ## this is a model using the hypothesis specific estimate as an offset and then inferring a
        ## set specific intercept
        setformula = eval(parse(text = paste('count', " ~ ", paste(c('offset(count.pred)', names(sets)), collapse = "+"), '-1'))) 

        ## reduce setdata to only rows (hypotheses) that belong to at least one set (speed things up)
        setdata = setdata[rowSums(setdata[, -c(1:2), drop = FALSE])>0 & !is.na(setdata$count), ]
        if (nb)
          {
            gset = tryCatch(glm.nb.fh(setformula, data = setdata, maxit = iter, theta = structure(g$theta, SE = g$SE.theta)),
                            error = function(e) NULL)
          } else
          {
            gset = glm(setformula, data = as.data.frame(setdata), family = poisson)
          }

        if (!is.null(gset))
          {
            tmpres = dflm(gset)[names(sets), ]
          }
        else
        {
          tmpres = data.table(name = names(sets), method = as.numeric(NA), p = as.numeric(NA), p.left = as.numeric(NA), p.right = as.numeric(NA), estimate = as.numeric(NA), ci.lower = as.numeric(NA), ci.upper = as.numeric(NA), effect = as.character(NA))
                              
        }
        tmpres[, n := elementNROWS(sets)]
        return(tmpres)

        }

      setres = rbindlist(mclapply(1:length(sets), .score.sets, sets = sets, res = res, nb = nb, mc.cores = mc.cores))
      setres$name = setmap[.(setres$name), names] ## remap set names to original names
      setnames(setres, 'name', 'setname')
      setres$method = gsub('\\(.*$', '', setres$method)
      setres$fdr = signif(p.adjust(setres$p, 'BH'), 2) ## compute q value
    }

    return(list(res = as.data.table(res),model = g, setres = setres))

}



#' @name Cov
#' @title Cov
#' @description
#'
#' function to initialize Covariates for passing to FishHook object constructor.
#'
#' Can also be initiated by passing a vector of multiple vectors of equal length, each representing one of the internal variable names
#' You must also include a list containg all of the covariates (Granges, chracters, RLELists, ffTracks)
#'
#' Covariate serves to mask the underlieing list implemenations of Covariates in the FishHook Object.
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
#' @param data, a list of covariate data that can include any of the covariate classes (GRanges, ffTrack, RleList, character)
#' @return Covariate object that can be passed directly to the FishHook object constructor
#' @author Zoran Z. Gajic
#' @import R6
#' @export
Cov = function(data = NULL, field = as.character(NA), name = as.character(NA), pad = 0, type = as.character(NA),
               signature = as.character(NA), na.rm = NA, grep = NA){
  Covariate$new(name = name, data = data, pad = pad, type = type, signature = signature,
                field = field, na.rm = na.rm, grep = grep)
}
 
#' @name Covariate
#' @title title
#' @description
#'
#' Stores Covariates for passing to FishHook object constructor.
#'
#' Can also be initiated by passing a vector of multiple vectors of equal length, each representing one of the internal variable names
#' You must also include a list containg all of the covariates (Granges, chracters, RLELists, ffTracks)
#'
#' Covariate serves to mask the underlieing list implemenations of Covariates in the FishHook Object.
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
#' @param data, a list of covariate data that can include any of the covariate classes (GRanges, ffTrack, RleList, character)
#' @return Covariate object that can be passed directly to the FishHook object constructor
#' @author Zoran Z. Gajic
#' @import R6
#' @export
Covariate = R6::R6Class('Covariate',
    public = list(

      ## See the class documentation
      initialize = function(data = NULL, field = as.character(NA), name = as.character(NA),
                            pad = 0, type = 'numeric', signature = as.character(NA),
                            na.rm = as.logical(NA), grep = as.logical(NA)){

        ##If data are valid and are a list of tracks concatenate with any premade covs
        if(is.null(data))
        {
          self$data = NULL
          self$names = name
          self$type = type
          self$signature = signature
          self$field = field
          self$pad = pad
          self$na.rm = na.rm
          self$grep = grep
          return()
                                        #        stop('No data provided to Covariate instantiation.  Data must be provided as GRanges, filepath to supported UCSC format (.bed, .bw, .bigWig), or path to .rds file of GRanges.')
        }

        if(class(data) != 'list'){
          data = list(data)
        }            

        ## replicate params and data if necessary
        params = data.table(id = 1:length(data), field = field, name = name, pad = pad, type = type, signature = signature, na.rm = na.rm, grep = grep)

        if (length(data)==1 & nrow(params)>1)
          data = rep(data, nrow(params))

        self$data = data

        params$dclasses = sapply(self$data, class)

        if (any(ix <<- params$dclasses == 'character'))
        {
          if (any(iix <<- !file.exists(unlist(self$data[ix]))))
          {
            stop(sprintf('Some covariate files not found:\n%s', paste(unlist(self$data[ix][iix]), collapse = ',')))
          }
        }

        ## for any GRanges data that are provided where there is more than one metadata
        ## column, there should be a field given, otherwise we complain
        dmeta = NULL
        if (any(ix <<- params$dclasses != 'character'))
        {
          dmeta = lapply(self$data[ix], function(x) names(values(x)))
        }
        
        ## we require field to be specified if GRanges have more than one metadata column, otherwise
        ## ambiguous
        if (length(dmeta)>0)
        {
          ## check to make sure that fields actually exist in the provided GRanges arguments
          found = mapply(params$field[ix], dmeta, FUN = function(x,y) ifelse(is.na(x), NA, x %in% y))

          if (any(!found, na.rm = TRUE))
          {
            stop('Some provided Covariate fields not found in their corresponding GRanges metadata, please check arguments')
          }
        }                

        if (na.ix <<- any(is.na(params$name)))
        {
          params[na.ix, name := ifelse(!is.na(field), field, paste0('Covariate', id))]
          params[, name := dedup(name)]       
        }

        ## label any type = NA covariates for which a field has not been specified
        ## as NA by default
        if (any(na.ix <<- !is.na(params$field) & is.na(params$type)))
        {
          params[na.ix, type := 'numeric']
        }

        if (any(na.ix <<- is.na(params$field) & is.na(params$type)))
        {
          params[na.ix, type := 'interval']
        }

        ## check names to make sure not malformed, i.e. shouldn't start with number or contain spaces or special
        ## characters

        if (any(iix <<- grepl('\\W', params$name)))
        {
          warning('Replacing spaces and special characters with "_" in Covariate names')
          params$names[iix] = gsub('\\W+', '_', params$name[iix])
        }

        if (!is.null(params$name))
          {
            if (any(iix <<- duplicated(params$name)))
            {
              warning('Deduping covariate names')
              params$name = dedup(params$name)
            }
          }

        self$names = params$name
        self$type = params$type
        self$signature = params$signature
        self$field = params$field
        self$pad = params$pad
        self$na.rm = params$na.rm
        self$grep = params$grep            
    },

    ## Params:
    ## ... Other Covariates to be merged into this array, note that it can be any number of Covariates
    ## Return:
    ## A single Covariate object that contains the contents of self and all passed Covariates
    ## UI:
    ## None
    ## Notes:
    ## This is linked to the c.Covariate override and the c.Covariate override should be used preferentially over this
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
    ## range, a numeric vector of the covariates to include. e.g. if the Covariate contains the covariates (A,B,C) and the range is c(2:3),
    ## this indicates you wish to get a Covariate containing (B,C). NOTE THAT THIS DOES NOT RETURN A NEW COV_ARR, IT MODIFIES THE CURRENT.
    ## Return:
    ## None, this modifies the Covariate on which it was called
    ## UI:
    ## None
    ## Notes:
    ## If you want to create a new Covariate containing certain covariates, use the '[' operator, e.g. Covariate[2:3]
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
    ## A list of lists where each internal list corresponds to the covariate and is for use internally in the annotate.hypotheses function
    ## The list representation of the covariate will contain the following variables: type, signature, pad, na.rm, field, grep
    ## UI:
    ## None
    toList = function(...){
        if(length(private$pCovs) == 0){
            return(list())
        }
        out = lapply(c(1:length(private$pCovs)), function(x){
            return (list(track = private$pCovs[[x]],
                         type = private$ptype[x],
                         signature = private$psignature[x],
                         pad = private$ppad[x],
                         na.rm = private$pna.rm[x],
                         field = private$pfield[x],
                         grep = private$pgrep[x]))
        })
        names(out) = private$pnames
        return(out)

        },

    ## Params:
    ## No params required, included arguments will be ignored.
    ## Return:
    ## Nothing
    ## UI:
    ## Prints information about the Covariate to the console with all of covariates printed in order with variables printed alongside each covariate
    print = function(...){
        if(length(private$pCovs) == 0){
            fmessage('Empty Covariate Object')
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
            ## A character vector of field names for use with numeric covariates, see the Covariate class definition for more info
            pfield = c(),
            ## A numeric vector of paddings for each covariate, see the 'pad' param in Covariate class definition for more info
            ppad = c(),
            ## A logical vector for each covariate, see the 'na.rm' param in Covariate class definition for more info
            pna.rm = c(),
            ##  A chracter vector for each covariate, see the 'grep' param in Covariate class definition for more info
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

                  if (any(iix <<- grepl('\\W', value)))
                  {
                    warning('Replacing spaces and special characters with "_" in provided Covariate names')
                    value[iix] = gsub('\\W+', '_', value[iix])
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
                        stop("Error: pad must be of class numeric")
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
                        stop("Error: na.rm must be of class logical")
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
            data = function(value) {
                if(!missing(value)){
                    private$pCovs = value
                    return(private$pCovs)
                } else{
                    return(private$pCovs)
                }
            }
    ),
)

    



#' @name c.Covariate
#' @title title
#' @description
#'
#' Override the c operator for covariates so that you can merge them like a vector
#'
#' @param ... A series of Covariates, note all objects must be of type Covariate 
#' @return Covariate object that can be passed directly into the FishHook object constructor that contains all of the Covariate covariates
#' Passed in the ... param
#' @author Zoran Z. Gajic
#' @export
'c.Covariate' = function(...){

    ##Ensure that all params are of type Covariate
    Covariates = list(...)
    isc = sapply(Covariates, function(x)  class(x)[1] == 'Covariate')

    if(any(!isc)){
        stop('Error: All inputs must be of class Covariate.')
    }

    ## Merging vars of the covariates
  names  = unlist(lapply(Covariates, function(x) x$names))
  type  = unlist(lapply(Covariates, function(x) x$type))
  signature  = unlist(lapply(Covariates, function(x) x$signature))
  field  = unlist(lapply(Covariates, function(x) x$field))
  pad  = unlist(lapply(Covariates, function(x) x$pad))
  na.rm  = unlist(lapply(Covariates, function(x) x$na.rm))
  grep  = unlist(lapply(Covariates, function(x) x$grep))
  
  ## Merging Covariates
  covs = lapply(Covariates, function(x) x$data)
  Covs = unlist(covs, recursive = F)
  
  ##Creating a new Covariate and assigning all of the merged variables to it
  ret = Covariate$new(data = Covs, name = names, type = type, signature = signature, field = field, pad = pad, na.rm = na.rm, grep = grep)

    return(ret)
}




#' @name [.Covariate
#' @title title
#' @description
#'
#' Overrides the subset operator x[] for use with Covariate to allow for vector like subsetting
#'
#' @param obj Covariate This is the Covariate to be subset
#' @param range vector This is the range of Covariates to return, like subsetting a vector. e.g. c(1,2,3,4,5)[3:4] == c(3,4)
#' @return A new Covariate object that contains only the Covs within the given range
#' @author Zoran Z. Gajic
#' @export
'[.Covariate' = function(obj, range){
  if (any(is.na(range)))
    stop('NA in subscripts not allowed')

  if (any(range>length(obj)))
    stop('Subscript out of bounds')

  ##Clone the object so we don't mess with the original
  ret = obj$clone()
  ##Call the subset function of the Covariate class that will modify the cloned Covariate
  ret$subset(range)
  return (ret)
}


#' @name names.Covariate
#' @title title
#' @description
#'
#' Overrides the names function for Covariate object
#'
#' @param obj Covariate This is the Covariate whose names we are extracting' 
#' @return names of covariates
#' @author Zoran Z. Gajic
#' @export
'names.Covariate' = function(x){
    ##Call the subset function of the Covariate class that will modify the cloned Covariate
    return (x$names)
}



#' @name names.FishHook
#' @title title
#' @description
#'
#' Overrides the names function for FishHook object, i.e. the names of its covariates
#'
#' @param obj FishHook This is the FishHook whose names we are extracting
#' @return names of covariates
#' @author Zoran Z. Gajic
#' @export
'names.FishHook' = function(x){
    ##Call the subset function of the Covariate class that will modify the cloned Covariate
    return (x$covariates$names)
}



#' @name length.Covariate
#' @title title
#' @description
#'
#' Overrides the length function 'length(Covariate)' for use with Covariate
#'
#' @param obj Covariate object that is passed to the length function
#' @return number of covariates contained in the Covariate object as defined by length(Covariate$data)
#' @author Zoran Z. Gajic
#' @export
'length.Covariate' = function(obj,...){
    return(length(obj$data))
}







#' @name FishHook
#' @title title
#' @description
#'
#' Stores Events, Hypotheses, Eligible, Covariates.
#'
#' @param hypotheses Examples of hypotheses are genes, enhancers, or even 1kb tiles of the genome that we can then convert into a rolling/tiled window. This param must be of class "GRanges".
#' @param events Events are the given mutational regions and must be of class "GRanges". Examples of events are SNVs (e.g. C->G) somatic copy number alterations (SCNAs), fusion events, etc.
#' @param eligible Eligible regions are the regions of the genome that have enough statistical power to score. For example, in the case of exome sequencing where all regions are not equally
#' represented, eligible can be a set of regions that meet an arbitrary exome coverage threshold. Another example of when to use eligibility is in the case of whole genomes,
#' where your hypotheses are 1kb tiles. Regions of the genome you would want to exclude in this case are highly repetative regions such as centromeres, telomeres, and satelite repeates.
#' This param must be of class "GRanges".
#' @param covariates Covariates are genomic covariates that you belive will cause your given type of event (mutations, CNVs, fusions, case control samples) that are not linked to the process you are
#' investigating (e.g. cancer drivers). In the case of cancer drivers, we are looking for regions that are mutated as part of cancer progression. As such, regions that are more suceptable to
#' random mutagenesis such as late replicating or non-expressed region (transcription coupled repair) could become false positives. Including covariates for these biological processes will
#' reduce thier visible effect in the final data. This param must be of type "Covariate".
#' @param out.path A character that will indicate a system path in which to save the results of the analysis.
#' @param use_local_mut_density A logical that when true, creates a covariate that will represent the mutational density in the genome, whose bin size will be determined by local_mut_density_bin.
#' This covariate can be used when you have no other covariates as a way to correct for variations in mutational rates along the genome under the assumption that driving mutations
#' will cluster in local regions as opposed to global regions. This is similar to saying, in the town of foo, there is a crime rate of X that we will assume to be the local crime rate
#' If a region in foo have a crime rate Y such that Y >>>>> X, we can say that region Y has a higher crime rate than we would expect.
#' @param local_mut_density_bin A numeric value that will indicate the size of the genomic bins to use if use_local_mut_density = TRUE. Note that this paramter should be a few orders of
#' magnitude greater than the size of your targetls
#'
#' e.g. if your hypotheses are 1e5 bps long, you may want a local_mut_density_bin of 1e7 or even 1e8
#' @param genome A character value indicating which build of the human genome to use, by default set to hg19
#' @param mc.cores A numeric value that indicates the amount of computing cores to use when running fishHook. This will mainly be used during the annotation step of the analysis, or during
#' initial instantiation of the object if use_local_mut_density = T
#' @param na.rm A logical indicating how you handle NAs in your data, mainly used in fftab and gr.val, see these function documentations for more information
#' @param pad A numeric indicating how far each covariate range should be extended, see Covariate for more information, not that this will only be used if atleast on of the
#' Covariates have pad = NA
#' @param vebose A logical indicating whether or not to print information to the console when running FishHook
#' @param max.slice integer Max slice of intervals to evaluate with  gr.val (default = 1e3)
#' @param ff.chunk integer Max chunk to evaluate with fftab (default = 1e6)
#' @param max.chunk integer gr.findoverlaps parameter (default = 1e11)
#' @param idcol A character, that indicates the column name containing the patient ids, this is for use in conjunction with idcap. If max patientpergene is specified and
#' and the column referenced by idcol exists, we will limit the contributions of each patient to each target to idcap. e.g. if Patient A has 3 events in target A and Patient B
#' has 1 event in target A, and idcap is set to 2, with thier ID column specified, target A will have a cournt of 3, 2 coming from patient A and 1 coming from patient B
#' @param idcap a numeric that indicates the max number of events any given patient can contribute to a given target. for use in conjction with idcol. see idcol for more info.
#' @param weightEvents a logical that indicates if the events should be weighted by thier overlap with the hypotheses. e.g. if we have a SCNA spanning 0:1000 and a target spanning 500:10000, the overlap
#' of the SCNA and target is 500:1000 which is half of the original width of the SCNA event. thus if weightEvent = T, we will credit a count of 0.5 to the target for this SCNA. This deviates from
#' the expected input for the gamma poisson as the gamma poisson measures whole event counts.
#' @param nb boolean negative binomial, if false then use poisson
#' @return FishHook object ready for annotation/scoring.
#' @author Zoran Z. Gajic
#' @importFrom R6 R6Class
#' @export
Fish = function(hypotheses = NULL, events = NULL, covariates = NULL, eligible = NULL, out.path = NULL, 
            use_local_mut_density = FALSE, local_mut_density_bin = 1e6, genome = 'BSgenome.Hsapiens.UCSC.hg19::Hsapiens',
            mc.cores = 1, na.rm = TRUE, pad = 0, verbose = TRUE, max.slice = 1e4, ff.chunk = 1e6, max.chunk = 1e11, idcol = NULL,
            idcap = Inf, weightEvents = FALSE, nb = TRUE)
{
  FishHook$new(hypotheses = hypotheses, out.path = out.path, eligible = eligible, events = events, covariates = covariates,
               use_local_mut_density = use_local_mut_density, local_mut_density_bin = local_mut_density_bin, genome = genome,
               mc.cores = mc.cores, na.rm = na.rm, pad = pad, verbose = verbose, max.slice = max.slice, ff.chunk = ff.chunk, max.chunk = max.chunk, idcol = idcol,
               idcap = idcap, weightEvents = weightEvents, nb = nb)    
}

#' @name FishHook
#' @title title
#' @description
#'
#' Stores Events, Hypotheses, Eligible, Covariates.
#'
#' @param hypotheses Examples of hypotheses are genes, enhancers, or even 1kb tiles of the genome that we can then convert into a rolling/tiled window. This param must be of class "GRanges".
#' @param events Events are the given mutational regions and must be of class "GRanges". Examples of events are SNVs (e.g. C->G) somatic copy number alterations (SCNAs), fusion events, etc.
#' @param eligible Eligible regions are the regions of the genome that have enough statistical power to score. For example, in the case of exome sequencing where all regions are not equally
#' represented, eligible can be a set of regions that meet an arbitrary exome coverage threshold. Another example of when to use eligibility is in the case of whole genomes,
#' where your hypotheses are 1kb tiles. Regions of the genome you would want to exclude in this case are highly repetative regions such as centromeres, telomeres, and satelite repeates.
#' This param must be of class "GRanges".
#' @param covariates Covariates are genomic covariates that you belive will cause your given type of event (mutations, CNVs, fusions, case control samples) that are not linked to the process you are
#' investigating (e.g. cancer drivers). In the case of cancer drivers, we are looking for regions that are mutated as part of cancer progression. As such, regions that are more suceptable to
#' random mutagenesis such as late replicating or non-expressed region (transcription coupled repair) could become false positives. Including covariates for these biological processes will
#' reduce thier visible effect in the final data. This param must be of type "Covariate".
#' @param out.path A character that will indicate a system path in which to save the results of the analysis.
#' @param use_local_mut_density A logical that when true, creates a covariate that will represent the mutational density in the genome, whose bin size will be determined by local_mut_density_bin.
#' This covariate can be used when you have no other covariates as a way to correct for variations in mutational rates along the genome under the assumption that driving mutations
#' will cluster in local regions as opposed to global regions. This is similar to saying, in the town of foo, there is a crime rate of X that we will assume to be the local crime rate
#' If a region in foo have a crime rate Y such that Y >>>>> X, we can say that region Y has a higher crime rate than we would expect.
#' @param local_mut_density_bin A numeric value that will indicate the size of the genomic bins to use if use_local_mut_density = TRUE. Note that this paramter should be a few orders of
#' magnitude greater than the size of your targetls
#'
#' e.g. if your hypotheses are 1e5 bps long, you may want a local_mut_density_bin of 1e7 or even 1e8
#' @param genome A character value indicating which build of the human genome to use, by default set to hg19
#' @param mc.cores A numeric value that indicates the amount of computing cores to use when running fishHook. This will mainly be used during the annotation step of the analysis, or during
#' initial instantiation of the object if use_local_mut_density = T
#' @param na.rm A logical indicating how you handle NAs in your data, mainly used in fftab and gr.val, see these function documentations for more information
#' @param pad A numeric indicating how far each covariate range should be extended, see Covariate for more information, not that this will only be used if atleast on of the
#' Covariates have pad = NA
#' @param vebose A logical indicating whether or not to print information to the console when running FishHook
#' @param max.slice integer Max slice of intervals to evaluate with  gr.val (default = 1e3)
#' @param ff.chunk integer Max chunk to evaluate with fftab (default = 1e6)
#' @param max.chunk integer gr.findoverlaps parameter (default = 1e11)
#' @param idcol A character, that indicates the column name containing the patient ids, this is for use in conjunction with idcap. If max patientpergene is specified and
#' and the column referenced by idcol exists, we will limit the contributions of each patient to each target to idcap. e.g. if Patient A has 3 events in target A and Patient B
#' has 1 event in target A, and idcap is set to 2, with thier ID column specified, target A will have a cournt of 3, 2 coming from patient A and 1 coming from patient B
#' @param idcap a numeric that indicates the max number of events any given patient can contribute to a given target. for use in conjction with idcol. see idcol for more info.
#' @param weightEvents a logical that indicates if the events should be weighted by thier overlap with the hypotheses. e.g. if we have a SCNA spanning 0:1000 and a target spanning 500:10000, the overlap
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
        initialize = function(hypotheses, eligible = NULL, events = NULL, covariates = NULL, out.path = NULL, 
            use_local_mut_density = FALSE, local_mut_density_bin = 1e6, genome = 'BSgenome.Hsapiens.UCSC.hg19::Hsapiens',
            mc.cores = 1, na.rm = TRUE, pad = 0, verbose = TRUE, max.slice = 1e4, ff.chunk = 1e6, max.chunk = 1e11, idcol = NULL,
            idcap = Inf, weightEvents = FALSE, nb = TRUE){

            #Make sure the format of hypotheses, events, eligible, and covarates is correct
          if (is.null(hypotheses))
            stop('Hypotheses must be defined')

            ##Hypotheses
            if(!((class(hypotheses) == 'GRanges') || class(hypotheses) == 'character')  && !is.null(hypotheses)){
                stop('Error: hypotheses must be of class GRanges or character')
            }
            
            if(class(hypotheses) == 'character'){
                self$hypotheses = hypotheses
                hypotheses = self$hypotheses
            }
        
        
            ##Events
            if(!(class(events) == 'GRanges')  && !is.null(events)){
                stop('Error: events must be of class GRanges')
            }
            
            
            ##Eligible
            if(!(class(eligible) == 'GRanges') && !is.null(eligible)){
                stop('Error: eligible must be of class GRanges')
            }

            
            ##Covariates
            if(!(class(covariates) == 'Covariate')  && !is.null(covariates)){
                stop('Error: covariates must be of class Covariate')
            }
            
            
             ## This next portion checks to make sure that the seqlevels are in the same format
            if(!is.null(covariates)){
                ## Gets whther seqlevels of covariates are chr or not chr
                seqLevelsStatus_Covariates = covariates$chr()
                ## Warns if there is a heterogenetiy of seqlevels (chr or not)
                if(length(unique(seqLevelsStatus_Covariates)) > 1){
                    warning('Warning:Covariates appears to have mismatched seqlevels, make sure all Covariates have seqlevels that start with chr or do not', call.=TRUE)
                }
            }

            ## gets the seqlevels and looks for chr to indicate UCSC format
          seqLevelsStatus_Hypotheses = any(grepl('chr', GenomeInfoDb::seqlevels(hypotheses)))

          seqLevelsStatus_Events = seqLevelsStatus_Hypotheses  
          if (!is.null(events))
            seqLevelsStatus_Events = any(grepl('chr', GenomeInfoDb::seqlevels(events)))

          if(!is.null(covariates)){
            if(any(!(seqLevelsStatus_Hypotheses %in% seqLevelsStatus_Covariates))){
              warning('Warning: seqlevels of Hypotheses and Covariates appear to be in different formats')
            }
          }

          if(!is.null(eligible)){
            seqLevelsStatus_Eligible = any(grepl("chr", GenomeInfoDb::seqlevels(eligible)))
            if(seqLevelsStatus_Hypotheses != seqLevelsStatus_Eligible){
              warning('Warning: seqlevels of Hypotheses and Eligible appear to be in different formats')
            }
          }
          if(seqLevelsStatus_Hypotheses != seqLevelsStatus_Events){
            warning('Warning: seqlevels of Hypotheses and Events appear to be in different formats')
          }

          ## This next portion checks to make sure there is atleast some overlap of seqlevels i.e. some mapability
          if (!is.null(events))
          {
            if(!any(GenomeInfoDb::seqlevels(hypotheses) %in% GenomeInfoDb::seqlevels(events))){
              stop('Error: there are no seqlevels of events that match hypotheses')
            }
          }
          
          if(!is.null(eligible)){
            if(!any(GenomeInfoDb::seqlevels(hypotheses) %in% GenomeInfoDb::seqlevels(eligible))){
              stop('Error: there are no seqlevels of eligible that match hypotheses')
            }
          }

          if(!is.null(covariates)){
            if(any(!(unlist(lapply(covariates$seqlevels(),function(x) any(x %in% GenomeInfoDb::seqlevels(hypotheses))))))){
              warning("Warning: atleast one of the covariates has no seqlevels in common with hypotheses")
            }
          }

          ## Initializes and Validates hypotheses
          ### MARCIN: we no longer allow hypotheses to be publicly reset, we will do this manually and check validity of hypotheses using function
          ##          self$hypotheses = hypotheses
          
          if (validate_hypotheses(hypotheses))
            private$phypotheses = hypotheses
          
          # Initializes and Validates out.path
          self$out.path = out.path

          private$pmc.cores = mc.cores
          private$pna.rm = na.rm
          private$ppad = pad
          private$pverbose = verbose
          private$pmax.slice = max.slice
          private$pff.chunk = ff.chunk
          private$pmax.chunk = max.chunk
          private$pidcol = idcol
          private$pidcap = idcap
          private$pweightEvents = weightEvents
          private$pnb = nb

          ## Initializes and Validates eligible (used to be self$eligible = , but this will trigger a re-annotation so avoiding)
          if(!(is.null(eligible))){
            if (validate_eligible(eligible))
              private$peligible = eligible
          }

          ## Initializes and Validates covariates
          if(is.null(covariates)){
            covariates = Covariate$new()
          }
          
          ## Initializes and Validates events (MARCIN: will call annotate now)
          ## remember we can have fishHook object with no events
          if (!is.null(events))
            {
              self$events = events
            }

          ## Initializes and Validates events (MARCIN: will call annotate as well)
          self$covariates = covariates
                   
          ##Creating the local mutational denisty track
          if(use_local_mut_density){
            Sys.setenv(DEFAULT_BSGENOME = genome)
            bins = gr.tile(hg_seqlengths(), local_mut_density_bin)
                f1 = FishHook$new(hypotheses = bins, events = events, eligible = eligible, mc.cores = mc.cores, na.rm = na.rm, verbose = verbose, max.slice = max.slice, ff.chunk = ff.chunk, max.chunk = max.chunk)
                f1$score()
                local_mut_density = seg2gr(f1$res)[,'count.density']
                cd = local_mut_density$count.density
                avg_cd = mean(cd, na.rm = T)
                cd[is.na(cd) | cd == Inf] = avg_cd
                local_mut_density$count.density = cd
                if(length(private$pcovariates$toList()) == 0 ){
                  private$pcovariates = c(Covariate$new(data = c(local_mut_density), type = c('numeric'), name = c("Local Mutation Density"), field = c("count.density")))
                } else{
                  private$pcovariates = c(Covariate$new(data = c(local_mut_density), type = c('numeric'), name = c("Local Mutation Density"), field = c("count.density")), private$pcovariates)
                }
            }

          ## private$pdata = annotate.hypotheses(hypotheses = private$phypotheses,
          ##                                     covered = private$peligible,
          ##                                     events = private$pevents,
          ##                                     mc.cores = mc.cores,
          ##                                     na.rm = na.rm,
          ##                                     pad = pad,
          ##                                     verbose = verbose,
          ##                                     max.slice = max.slice,
          ##                                     ff.chunk = ff.chunk,
          ##                                     max.chunk = max.chunk,
          ##                                     out.path = private$pout.path,
          ##                                     covariates = private$pcovariates$toList(),
          ##                                     idcol = idcol,
          ##                                     idcap = idcap,
          ##                                     weightEvents = weightEvents)

          private$pstate = "Annotated"
        },


        ## ## Params:
        ## ## mc.cores, see FishHook class documentation for more info
        ## ## na.rm, see FishHook class documentation for more info
        ## ## pad, see FishHook class documentation for more info
        ## ## verbose, see FishHook class documentation for more info
        ## ## max.slice, see FishHook class documentation for more info
        ## ## ff.chunk, see FishHook class documentation for more info
        ## ## max.chunk, see FishHook class documentation for more info
        ## ## idcol, see FishHook class documentation for more info
        ## ## idcap, see FishHook class documentation for more info
        ## ## Return:
        ## ## None
        ## ## UI:
        ## ## If verbose = T, will print updates as the annotation proceeds
        ## ## Notes:
        ## ## This function changes the internal state of the fishHook object and sets the state to 'Annotated'
        ## annotate = function(mc.cores = private$pmc.cores, na.rm = private$pna.rm, pad = private$ppad,
        ##     verbose = private$pverbose, max.slice = private$pmax.slice, ff.chunk = private$pff.chunk,
        ##     max.chunk = private$pmax.chunk, idcol = private$pidcol, idcap = private$idcap,
        ##     weightEvents = private$pweightEvents){

        ##     if(private$pstate == "Scored"){
        ##         stop("Error: You have a scored object already created. If you want to re-run the analysis you can clear the scored and annotated objects using fish$clear()")
        ##     }

        ##     private$pdata = annotate.hypotheses(hypotheses = private$phypotheses,
        ##         covered = private$peligible,
        ##         events = private$pevents,
        ##         mc.cores = mc.cores,
        ##         na.rm = na.rm,
        ##         pad = pad,
        ##         verbose = verbose,
        ##         max.slice = max.slice,
        ##         ff.chunk = ff.chunk,
        ##         max.chunk = max.chunk,
        ##         out.path = private$pout.path,
        ##         covariates = private$pcovariates$toList(),
        ##         idcol = idcol,
        ##         idcap = idcap,
        ##         weightEvents = weightEvents)

        ##     private$pstate = "Annotated"

        ## },



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
            targ = paste('Contains' , length(private$phypotheses), "hypotheses." ,collapse = "")
            eve = paste('Contains', length(private$pevents), "events to map to hypotheses.", collapse = "")
            if(is.null(private$peligible)){
                elig = "All regions are eligible."
            } else{
                elig = sprintf("Spanning %s MB of eligible territory.", round(sum(as.numeric(width(private$peligible)))/1e6, 2))
            }
            if(is.null(private$pcovariates$names)){
                covs = "No covariates will be used."
            } else{          
                cov.names = private$pcovariates$names
                if(length(cov.names) > 10){
                    covs = paste('Will use',length(cov.names),'covariates.', collapse = "")
                }
                else{
                    covs = cov.names
                }
            }
            meta = paste('Hypotheses contains', ncol(values(private$phypotheses)), 'metadata columns.')
            state = paste('Current State:', private$pstate)
            cat(targ, eve, elig, 'Covariates:', covs, meta, state, sep = '\n', collapse = '\n')
        },

        ## Params:
        ## hypotheses, a GRanges that is the output of annotate.hypotheses. note that this is for admin degbugging and
        ## should never be touched by the user unless you absoltely know exactly what you are doing and why you are doing it.
        ## by, character vector with which to split into meta-territories (default = NULL)
        ## fields, a character vector indicating which columns to be used in aggregateion by default all meta data
        ## fields of hypotheses EXCEPT reserved field names $coverage, $counts, $query.id (default = NULL)
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
        aggregate = function(hypotheses = private$pdata, by = NULL, fields = NULL, rolling = NULL, disjoint = TRUE, na.rm = FALSE, FUN = list(), verbose = private$pverbose){

            if(private$pstate == "Initialized"){
                self$annotate()
            }

            agg = aggregate.hypotheses(hypotheses,
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


      ## this function will append novel covariates from a Covariate or a fishHook object
      ## if from a fishHook object, then it will aggregate over the $data of that object
      ## i.e. not by doing a fresh aggregation over that fishHook objects covariates
      merge = function(data){
          if (is(data, 'Covariate'))
          {

            if (private$pverbose)
            {
              fmessage('Appending new Covariates to this fishHook object')
            }

            ## take care of any duplicated names among covariates
            newnames = c(private$pcovariates$names, data$names)
            if (any(duplicated(newnames)))
            {
              newnames = dedup(newnames)
              if (length(private$pcovariates)>0)
              {
                newnames = newnames[-c(1:length(private$pcovariates$names))]
              }
              data$names = newnames
            }            

            tmp.pdata = annotate.hypotheses(hypotheses = private$phypotheses,
                                            covered = private$peligible,
                                            mc.cores = private$pmc.cores,
                                            na.rm = private$pna.rm,
                                            pad = private$ppad,
                                            verbose = private$pverbose,
                                            max.slice = private$pmax.slice,
                                            ff.chunk = private$pff.chunk,
                                            max.chunk = private$pmax.chunk,
                                            out.path = private$pout.path,
                                            covariates = data$toList(),
                                            weightEvents = private$pweightEvents)
            BASIC.COLS = c('count', 'coverage', 'query.id')

            private$pcovariates = c(private$pcovariates, data)

            values(private$pdata) = cbind(values(private$pdata),
                                          values(tmp.pdata)[, setdiff(names(values(tmp.pdata)), BASIC.COLS), drop = FALSE])

            self$clear()          
          }
          else if (is(data, 'FishHook')) ## in this case we will just merge new covariates from the provided fishHook objects' $data matrix
        {
          cov.cols = data$covariates$names
          tmpdata = data$data[, cov.cols]

          if (length(cov.cols)==0)
          {
            warning('Attempting to append from fishHook object with empty covariates')
          }
          else
            {
              newnames = dedup(c(private$pcovariates$names, cov.cols))
              cov.cols = newnames[(length(private$pcovariates)+1):length(newnames)]
              
              names(values(tmpdata)) = cov.cols
              
              newdat = gr.val(self$data[, c()], tmpdata, val = cov.cols, na.rm = private$pna.rm, mc.cores = private$pmc.cores,
                              max.slice = private$pmax.slice)
              
              newcovs = do.call('c', lapply(cov.cols, function(x) Covariate$new(data = tmpdata, field = x, name = x, na.rm = TRUE, type = 'numeric')))
              private$pcovariates = c(private$pcovariates, newcovs)
              values(private$pdata) = cbind(values(private$pdata), values(newdat))
            }
        }            
        },

        ## Params:        
        ## hypotheses,  a GRanges that is the output of annotate.hypotheses. note that this is for admin degbugging and
        ## should never be touched by the user unless you absoltely know exactly what you are doing and why you are doing it.
        ## annotated hypotheses with fields $coverage, optional field, $count and additional numeric covariates
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
        score = function(nb = private$pnb, sets = private$psets, verbose = private$pverbose,  iter = 200, subsample = 1e5, seed = 42, p.randomize = TRUE, model = NULL){

            if(private$pstate == "Initialized"){
                self$annotate()
            }

          if (is.null(private$pevents))
            {
              stop('fishHook object has not been provided with events, please provide events (e.g. fish$events = events) and rescore')
            }

            ## If we are aggregated we should score that, if we are not we should score anno
            if(private$pstate == "Aggregated"){
                ## ##Rolling yeilds a GRanges
                ## if(class(private$paggregated) == 'GRanges'){
                ##     covs = names(values(private$paggregated))                    
                ## }
                ## ##by yeilds a GRangesList
                ## else{
                ##     covs = names(values(private$paggregated[[1]]))
                ## }
                targ = private$paggregated
                covs = c()
            } else{
                targ = private$pdata
                covs = names(values(private$pdata))
            }
   
          ## Scoring
          score = score.hypotheses(targ,
                                   sets = sets,
                                   covariates = covs,
                                   return.model = TRUE,
                                   nb = nb,
                                   verbose = verbose,
                                   iter = iter,
                                   subsample = subsample,
                                   seed = seed,
                                   classReturn = TRUE,
                                   model = model,
                                   mc.cores = private$pmc.cores,
                                   p.randomize = p.randomize)


          private$pscore = score$res
          private$pmodel = score$model

          if (!is.null(sets))
          {
            private$psets = sets
            private$psetscore = score$setres
          }

          private$pstate = 'Scored'
        },

      ## subset
      ## this subsets the intervals and hypotheses in the fishHook object
      subset = function(i = NULL, j = NULL){
        if (!is.null(j))
        {
          if (any(j>length(private$pcovariates)))
            stop('Covariate index out of bounds when subsetting FishHook object')

          if (is.character(j))
          {
            j = match(j, private$pcovariates$names)
            if (any(is.na(j)))
                stop('some of the provided column indices were not found among covariate names')
          }

          BASIC.COLS = c('count', 'coverage', 'query.id')
          private$pcovariates = private$pcovariates[j]
          private$pdata =  private$pdata[, c(BASIC.COLS, private$pcovariates$names)]
        }

        if (!is.null(i))
        {
          if (any(i>length(private$phypotheses)))
            stop('Hypothesis index out of bounds when subsetting FishHook object')

          private$pdata =  private$pdata[i, ]
          private$phypotheses = private$phypotheses[i]
          if (!is.null(private$pscore))
          {
            warning('Resetting hypothesis scores since object has been subsettted, please run $score() to get updated p values')
            private$pscore = private$pscore[i]
          }

          if (!is.null(private$psets))
          {
            map = data.table(old = i, new  = 1:length(i))
            setkey(map, old)
            setsdt = data.table(sid = factor(rep(names(private$psets), elementNROWS(private$psets))), ix = unlist(private$psets))
            setsdt$newix = map[.(setsdt$ix), new]
            setsdt = setsdt[!is.na(newix), ]
            if (nrow(setsdt)>0)
            {
              newsets = split(setsdt$newix, setsdt$sid)[names(private$psets)]
              names(newsets) = names(private$psets)
              private$psets = newsets
            }
            else
            {
              newsets = rep(list(as.integer(c())), length(sets))
              names(newsets) = names(private$psets)
              private$psets = newsets
            }
      
            warning('Resetting set scores since object has been subsettted, please run $score() to get updated set level p values')
            private$psetscore = NULL
          }
        }
      },
      
        ## Params:
        ## state, a character indicating which state to revert to, e.g. if you are 'Scored' you can revert to 'Annotated', 'Initialized', an possibly 'Aggregated'
        ## however if your state is 'Initialized' you cannot revert to a 'Scored' state
        ## Return:
        ## none
        ## UI:
        ## None
      clear = function(state = 'Initialized'){
        if (private$pstate == 'Scored')
          warning('Resetting scores since covariates re-defined, please run $score() to get updated p values')
        private$pmodel = NULL
        private$pstate = 'Annotated'
        private$pscore = NULL
        private$psetscore = NULL

            ## if(state == 'Initialized'){
            ##     private$pstate = 'Initialized'
            ##     private$pmodel = NULL
            ##     private$pscore = NULL
            ##     private$pdata = NULL
            ##     private$paggregated = NULL
            ##     return('Clear Completed')
            ## }
            ## if(state == 'Annotated'){
            ##     private$pstate = 'Annotated'
            ##     private$pmodel = NULL
            ##     private$pscore = NULL
            ##     private$paggregated = NULL
            ##     return('Clear Completed')
            ## }
            ## if(state == 'Aggregated'){
            ##     private$pstate = 'Aggregated'
            ##     private$pmodel = NULL
            ##     private$pscore = NULL
            ##     return('Clear Completed')
            ## }
##            return('Valid reversion state not specified. This is not a major error, just letting you know that nothing has been chaged')
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
      qqp = function(plotly = TRUE, columns = NULL, annotations = NULL, key = NULL, ...){
        if (is.null(self$res))
          stop('FishHook object not yet scored, so no qq_plot available')

            res = self$all

            if (is.null(columns))
              columns = names(values(private$phypotheses))
            
            columns = columns[columns %in% names(res)]
            annotation_columns = lapply(columns, function(x) as.character(unlist(res[,x,with=FALSE])))
            names(annotation_columns) = columns

            if(is.null(annotations) & is.null(columns)){
              
              names = res$name

                annotations = list(Hypothesis_ID = names,
                    Count = res$count,
                    Effectsize = round(res$effectsize,2),
                    fdr = res$fdr)

            }
            
            return(qqp(res$p ,annotations = c(annotations,annotation_columns), bottomrighttext = paste0('alpha =', round(self$model$theta,2)), gradient = list(Count = res$count), titleText = "" ,  plotly = plotly, key = key))

        }
    ),

    
    ## Private variables are internal variables that cannot be accessed by the user
    ## These variables will have active representations that the user can interact with the update
    ## and view these variables, all internal manipulations will be done with these private variables
    private = list(
        ## Genomic Ranges Object that Indicates the Hypotheses
        phypotheses = NULL,

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

        ##padding to use with the events, see annotate.hypotheses for more info
        ppad = 0,

        ##global verbose paramter
        pverbose = TRUE,

        ##see annotate.hypotheses for more info
        pmax.slice = 1e3,

        ##see annotate.hypotheses for more info
        pff.chunk = 1e6,

        ##see annotate.hypotheses for more info
        pmax.chunk = 1e11,

        ##see annotate.hypotheses for more info
        pidcol = NULL,

        ##see annotate.hypotheses for more info
        pidcap = Inf,

        ##see annotate.hypotheses for more info
        pweightEvents = FALSE,

        ##The variable containing the output of fishHook$annotate()
        pdata = NULL,

        ##The variable containing the output of fishHook$score()
        pscore = NULL,

      ##The variable containing the output of fishHook$score()
      psets = NULL,

      ##The variable containing the output of fishHook$score()
      psetscore = NULL,

        ##see score.hypotheses for more info
        pmodel = NULL,

        ##see score.hypotheses for more info
        preturn.model = TRUE,

        ##see score.hypotheses for more info
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

        ## Covariates = data
        ## Here we check to make sure that all data  are of class Covariate
        ## We then reset the object to its initialized state so as to not introduce incosistencies amongst variables
        covariates = function(value) {
            if(!missing(value)){
                if (!(class(value)[1] == 'Covariate')  & !is.null(value)){
                    stop('Error: covariates must be of class Covariate')
                }

                if (private$pverbose)
                {
                  fmessage('Aggregating covariates over eligible subset of hypotheses')
                }

                private$pcovariates = value


                tmp.pdata = annotate.hypotheses(hypotheses = private$phypotheses,
                                                covered = private$peligible,
                                                mc.cores = private$pmc.cores,
                                                na.rm = private$pna.rm,
                                                pad = private$ppad,
                                                verbose = private$pverbose,
                                                max.slice = private$pmax.slice,
                                                ff.chunk = private$pff.chunk,
                                                max.chunk = private$pmax.chunk,
                                                out.path = private$pout.path,
                                                covariates = private$pcovariates$toList(),
                                                weightEvents = private$pweightEvents)

                BASIC.COLS = c('count', 'coverage', 'query.id')
                other.cols = setdiff(names(values(tmp.pdata)), BASIC.COLS)

                ix = match(other.cols, names(values(tmp.pdata)))
                names(values(tmp.pdata))[ix] = value$names

                if (is.null(private$pdata))
                {
                  private$pdata = tmp.pdata
                }
                else
              {
                private$pdata = private$pdata[, BASIC.COLS]
                values(private$pdata) = cbind(values(private$pdata), values(tmp.pdata)[, setdiff(names(values(tmp.pdata)), BASIC.COLS), drop = FALSE])
              }
                self$clear()

                return(private$pcovariates)

            } else{
                return(private$pcovariates)
            }
        },

        ## Eligible
        ## Here we check to make sure that eligible is of class GRanges
        ## We then reset the object to its initialized state so as to not introduce incosistencies amongst variables
        eligible = function(value) {
          if(!missing(value)){
                if (validate_eligible(value))
                {
                  private$peligible = value
                }

                if (private$pverbose)
                {
                  fmessage('Recomputing counts and covariates over provided eligible regions')
                }

                private$pdata = annotate.hypotheses(hypotheses = private$phypotheses,
                                                    covered = private$peligible,
                                                    events = private$pevents,
                                                    mc.cores = private$pmc.cores,
                                                    na.rm = private$pna.rm,
                                                    pad = private$ppad,
                                                    verbose = private$pverbose,
                                                    max.slice = private$pmax.slice,
                                                    ff.chunk = private$pff.chunk,
                                                    max.chunk = private$pmax.chunk,
                                                    out.path = private$pout.path,
                                                    covariates = private$pcovariates$toList(),
                                                    idcol = private$pidcol,
                                                    idcap = private$pidcap,
                                                    weightEvents = private$pweightEvents)

                self$clear()

                return(private$peligible)

            } else{
                return(private$peligible)
            }
        },

        ## Hypotheses
        ## Here we check to make sure that hypotheses is of class GRanges or a chracter path and are not NULL
        ## We then reset the object to its initialized state so as to not introduce incosistencies amongst variables
        hypotheses = function(value) {
            if(!missing(value)){
              stop('Cannot reset hypotheses in existing FishHook object, please start with new object or subset this object using subsetting [] operator')    
            } else{
                return(private$phypotheses)
            }
        },

        ## Events
        ## Here we check to make sure that events is of class GRanges is not NULL
        ## then immediately reannotate
        events = function(value) {
            if(!missing(value)){
                if(!(class(value) == 'GRanges')){
                    stop('Error: Events must be of class GRanges')
                }
                
                events = value
                
                private$pevents = events
                
                ## update pdata
                if (private$pverbose)
                  {
                    fmessage('Tallying events over eligible subsets of hypotheses')
                  }

                tmp.pdata = annotate.hypotheses(hypotheses = private$phypotheses,
                                                events = private$pevents,
                                                covered = private$peligible, 
                                                mc.cores = private$pmc.cores,
                                                verbose = private$pverbose,
                                                max.slice = private$max.slice,
                                                ff.chunk = private$pff.chunk,
                                                max.chunk = private$pmax.chunk,
                                                idcol = private$pidcol,
                                                idcap = private$pidcap,
                                                weightEvents = private$pweightEvents)

                if (is.null(private$pdata))
                {
                    private$pdata = tmp.pdata
                  }
                else
                  {
                    private$pdata$count = tmp.pdata$count
                  }

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
                    tryCatch(saveRDS(private$phypotheses, paste(gsub('.rds', '', out.path), '.source.rds', sep = '')),
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
        data = function(value) {
            if(!missing(value)){
                if(!(class(value) == 'GRanges')  && !is.null(value)){
                    stop('Error: anno must be of class GRanges')
                } else{
                    warning('You are editing the annotated dataset generated by fishHook. Only do this if you are a fishHook pro!')
                }

                BASIC.COLS = c('count', 'coverage', 'query.id')
                if (!all(BASIC.COLS %in% names(values(value))))
                  stop('Provided GRanges must contain count, coverage, and query.id field, with additional <numeric> columns specifying covariate values')
                if (!identical(private$pdata[, c()], value[, c()]))
                  stop('Provided GRanges do not exactly match the hypotheses already present inside this fishHook object')

                cov.cols = setdiff(names(values(value)), BASIC.COLS)
                cov.classes = sapply(as.list(values(value))[cov.cols], class)
                if (!all(cov.classes == 'numeric'))
                  warning('Only <numeric> columns in provided granges can specify covariates, removing any non-numeric columns from given input')             

                cov.cols = cov.cols[cov.classes=='numeric']
                private$pcovariates = do.call('c', lapply(cov.cols, function(x) Covariate$new(data = value[, x], field = x, name = x, na.rm = TRUE, type = 'numeric')))
                private$pdata = value[, c(BASIC.COLS, cov.cols)]
                
                return(private$pdata)

            } else{
                return(private$pdata)
           }
        },

        ## res = results
        ## Here we check to make sure that scores is of class data.table
      res = function(value) {
          if(!missing(value)){
            stop('Scores cannot be edited, to rescore just run $score() function on object')
                ## if(!(class(value) == 'data.table')  && !is.null(value)){
                ##     stop('Error: score must be of class data.table')
                ## } else{
                ##     warning('Warning: You are editing the scored dataset generated by fishHook, if you are trying to change hypotheses use fish$hypotheses.')
                ## }

            ## private$pscore = value
            return(cbind(as.data.table(private$phypotheses), private$pscore[, 6:ncol(private$pscore)]))
          } else{
            if (is.null(private$pscore))
              stop('Model has not yet been scored, please run $score() and then retrieve results via $res')


            nms = names(private$pscore)[6:ncol(private$pscore)]

            return(cbind(as.data.table(private$phypotheses), private$pscore[, .(p, fdr, count, effectsize = round(effectsize,2), count.pred = signif(count.pred,3), count.density = signif(count.density,3), count.pred.density = signif(count.pred.density,3), p.neg, fdr.neg)]))
            }
      },

      ## all
        ## cannot be used for assigning data, can only be used for accessing a data.table containing merged scores and meta data
      all = function(value) {
        if(!missing(value)){
          stop('Error: This is solely for accessing results, only $eligible, $hypotheses, and $events can be set inside fishHook object.')
        } else{
            nms = names(private$pscore)[6:ncol(private$pscore)]

            nms.out = intersect(c('p', 'fdr', 'effectsize', 'count', 'count.pred', 'count.density', 'count.pred.density'), nms)
            nms = union(nms.out, nms)

            return(cbind(as.data.table(private$phypotheses), private$pscore[, nms, with = FALSE]))
            }
        },


        ## set scores
        ## Here we check to make sure that scores is of class data.table
        setres = function(value) {
          if(!missing(value)){
            stop('Set scores cannot be edited, to rescore just run $score(sets = mysets) or $sets = mysets function on object')
                ## if(!(class(value) == 'data.table')  && !is.null(value)){
                ##     stop('Error: score must be of class data.table')
                ## } else{
                ##     warning('Warning: You are editing the scored dataset generated by fishHook, if you are trying to change hypotheses use fish$hypotheses.')
                ## }

                ## private$pscore = value

                return(private$psetscore)

          } else{

            if (is.null(private$psetscore))
              stop('Sets results have not yet been generated, please set sets (via $sets = ) and/or run $score() to get set results')

              return(private$psetscore[, .(setname, p = signif(p, 3), fdr, effect, estimate = signif(estimate, 3), p.left = signif(p.left, 3),
                                           p.twosided = signif(as.numeric(p.twosided),3))])
            }
        },


        ## model
        model = function(value) {
            if(!missing(value)){

              warning('Warning: You are editing the regression model generated by fishHook.')

              if (!inherits(value, 'glm'))
                stop('Provided model must be glm model i.e. from another fishHook object')
              
              modelcovs = setdiff(rownames(summary(value)$coefficients), '(Intercept)')
              
              if (!identical(sort(modelcovs), sort(private$pcovariates$names)))
                stop('Mismatch between the names of covariates in provided model and the covariates in this fishHook object')

              private$pmodel = value

              if (private$pverbose)
              {
                fmessage('Rescoring fishHook data using provided model')
              }

              self$score(model = private$pmodel)

              return(private$pmodel)
            } else{
              return(private$pmodel)
            }
        },


        ## model
        sets = function(value) {
            if(!missing(value)){

              if (!inherits(value, 'list') & !all(sapply(value, class)== 'integer'))
                stop('Provided sets must be a (named) list of indices into hypotheses, each specifying a hypothesis set to be scored')
              
              if (any(sapply(value, function(x) max(c(x, 0), na.rm = TRUE))>length(self$hypotheses)))
                stop('Indices out of bounds for at least one of the provided sets')

              
              private$psets = value
              self$score(model = private$pmodel, sets = private$psets)

              return(self)
            } else{
              return(private$psets)
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
                if(!(class(value) %in% c('numeric', 'logical', 'integer'))  && !is.null(value) && length(value)==1){
                    stop('Error: verbose must be scalar of class logical, numeric, or integer')
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

        ## ff.chunk
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

        ## idcol
        ## character
        idcol = function(value) {
            if(!missing(value)){
                if(!(class(value) == "character")  && !is.null(value)){
                    stop('Error: idcol must be of class character')
                }

                if (!(value %in% names(values(events))))
                {
                  stop('Provided idcol does not exist as metadata column in $events')
                }

                private$pidcol = value

                ## re-annotate events
                self$events = self$events

                return(private$pidcol)

            } else{
                return(private$pidcol)
            }
        },

        ## idcap
        ## numeric, must be greater than 0, value is floored for safety
        idcap = function(value) {
            if(!missing(value)){
                if(!(class(value) == "numeric" | class(value) == 'integer')  && !is.null(value) && value > 0){
                    stop("Error: idcap must be of class numeric and non-negative")
                }

                private$pidcap = floor(value)

                ## re-annotate events
                self$events = self$events

                return(private$pidcap)

            } else{
                return(private$pidcap)
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

                self$events = self$events

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
        ## GRangesList containing aggregated hypotheses, you probably shouldn't be messing with this unless
        ## you really know what you're doing
        aggregated = function(value) {
          if(!missing(value)){
            stop('Aggregated cannot be set, can only be modified through $aggregate() method')

            return(private$paggregated)
            
          } else{
                return(private$paggregated)
            }
        }
    )
)




#' @name length.FishHook
#' @title title
#' @description
#'
#' Overrides the length function 'length(FishHook)' for use with FishHook 
#'
#' @param obj FishHook object that is passed to the length function
#' @return length of the hypotheses contained in the FishHook object
#' @author Zoran Z. Gajic
#' @export
'length.FishHook' = function(obj,...){
    return(length(obj$hypotheses))
}



#' @name dim.FishHook
#' @title title
#' @description
#'
#' Overrides the dim function 'dim(FishHook)' for use with FishHook 
#'
#' @param obj FishHook object that is passed to the length function
#' @return returns a numeric vector containing the lengths of various FishHook variables in the following order:
#' i : number of hypotheses
#' j : number of covariates
#' @author Zoran Z. Gajic
#' @export
'dim.FishHook' = function(obj,...){
    length_hypotheses = length(obj$hypotheses)
    length_covariates = length(obj$covariates)    
    return(c(length_hypotheses, length_covariates))
}


#' @name [.FishHook
#' @title title
#' @description
#'
#' Overrides the subset operator x[] for use with FishHook to allow for vector like subsetting, see fishHook demo for examples
#'
#' @param obj FishHook object This is the FishHookObject to be subset
#' @param i vector subset hypotheses
#' @param j vector subset covariates
#' @return A new FishHook object that contains only the given hypotheses and/or covariates
#' @author Zoran Z. Gajic
#' @export
'[.FishHook' = function(obj, i = NULL, j = NULL){
  ret = obj$clone()
  ret$subset(i, j)  
  return(ret)
}






#' @name qqp
#' @title qq plot given input p values
#' @description
#'
#' Generates R or Shiny quantile-quantile (Q-Q) plot given (minimally) an observed vector of p values, plotted their -log1 )quantiles against corresponding -log10
#' quantiles of the uniform distribution.
#'
#' @param obs vector of pvalues to plot, names of obs can be intepreted as labels, alternatively a data.frame / data.table with column $p, in which case the other (non $p) columns of obs are interpreted "annotations" in the html plot
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
#' @param annotations data.frame, data.table, or named list of vectors containing information to present as hover text (html widget), must be in same order and length as obs input, 
#' @param gradient named list that contains one vector that color codes points based on value, must bein same order as obs input
#' @param titleText title for plotly (html) graph only
#' @param subsample numeric (positive integer), number of points to use for plotting, will be taken randomly from the set of obs -> p values
#' @param key a character that is passed to the plotly function that will link each point to a give value. For example, if key is set to gene_name
#' The ploted points are refered to by thier gene_name. This is useful when integrating with shiny or any other tool that can integrate with plotly plots.
#' @import plotly
#' @author Marcin Imielinski, Eran Hodis, Zoran Z. Gajic
#' @export
qqp = function(obs, highlight = c(), exp = NULL, lwd = 1, col = NULL, col.bg = 'black', pch = 18, cex = 1, conf.lines = TRUE, max = NULL, max.x = NULL, bottomrighttext = NULL,
    max.y = NULL,  label = NULL, plotly = TRUE, annotations = list(), gradient = list(), titleText = "", subsample = NA, key = NULL,  ...)
{

  if (inherits(obs, 'data.frame') | inherits(obs, 'data.table') | inherits(obs, 'GenomicRanges'))
  {
    if (inherits(obs, 'GenomicRanges'))
    {
      obs = values(obs)
    }

    if (!("p" %in% colnames(obs)))
      stop('if obs is a data.frame or data.table or GenomicRanges then it must have a column $p corresponding to p value')


    annotations = as.list(as.data.frame(obs)[, setdiff(colnames(obs), 'p'), drop = FALSE])
    obs = obs$p

  }

  if (!is.null(bottomrighttext))
    bottomrighttext = paste0(bottomrighttext, '\n')

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
        if (is.null(bottomrighttext))
          bottomrighttext = ''
        legend('bottomright', sprintf('%slambda=\n %.2f', bottomrighttext, lambda), text.col = 'red', bty = 'n')
    } else{

        if(length(annotations) < 1){
            hover = do.call(cbind.data.frame, list(p = obs))
        } else{
            hover = cbind(as.data.frame(do.call(cbind, annotations)), data.frame(p = obs))
        }

        hover = as.data.table(hover)

      if (!is.null(key))
      {
        hover$key = key
      }

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
        hover$obs[ix1] = -log10(hover$p[ix1])
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
        } else if (!is.null(dat$grad)) {
            dat$grad = NULL
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
          if (!is.null(dat4$grad))
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
          if (!is.null(dat4$grad))
            {
              dat4$grad = NULL
            }


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
                y = (0.09 * max),
                text = paste(bottomrighttext, 'lambda =', sprintf('%.2f', signif(lambda,3)), collapse = ' '),
                font = list(color = 'red', size = 20),
                showarrow = FALSE,
                xref = 'x',
                yref = 'y'
            ),
            margin = list(t = 100),
            hovermode = 'compare')
    }
}


validate_eligible = function(value)
{
  if((!class(value) == 'GRanges') & !is.null(value)){
    stop('Error: eligible must be of class GRanges')
  }
  return(TRUE)
}

validate_hypotheses = function(value)
{
  if(!(class(value) == 'GRanges') && !(class(value) == 'character')){
    stop('Error: hypotheses must be of class GRanges or character')
  }
  
  hypotheses = value
  
  ## checks to see if hypotheses is a path & import if so
  if (is.character(hypotheses)){
    if (grepl('\\.rds$', hypotheses[1])){
      hypotheses = readRDS(hypotheses[1])
    } else if (grepl('(\\.bed$)', hypotheses[1])){
      require(rtracklayer)
      hypotheses = rtracklayer::import(hypotheses[1], (format = "BED"))
    }
  }
  
  ## checks to see if target contains any data
  if (length(hypotheses)==0){
    stop('Error: Must provide non-empty hypotheses')
  }
  
  ## Looks for a "name" field to index/Identify hypotheses by name
  ## If no such field is found creates a set of indexes
  if(is.null(hypotheses$name)){
    hypotheses$name = 1:length(hypotheses)
  }
  
  return(TRUE)
}


fmessage = function(..., pre = 'FishHook')
  message(pre, ' ', paste0(as.character(Sys.time()), ': '), ...)


dedup = function(x, suffix = '.')
{
  dup = duplicated(x);
  udup = setdiff(unique(x[dup]), NA)
  udup.ix = lapply(udup, function(y) which(x==y))
  udup.suffices = lapply(udup.ix, function(y) c('', paste(suffix, 2:length(y), sep = '')))
  out = x;
  out[unlist(udup.ix)] = paste(out[unlist(udup.ix)], unlist(udup.suffices), sep = '');
  return(out)
}




### modified glm.nb to allow setting of a (scalar or vector) theta (i.e. instead of inferring it)
###
glm.nb.fh = function (formula, data, weights, subset, na.action, start = NULL, 
    etastart, mustart, control = glm.control(...), method = "glm.fit", 
    model = TRUE, x = FALSE, y = TRUE, contrasts = NULL, ...,  theta = NULL,
    init.theta, link = log) 
{
    loglik <- function(n, th, mu, y, w) sum(w * (lgamma(th + 
        y) - lgamma(th) - lgamma(y + 1) + th * log(th) + y * 
        log(mu + (y == 0)) - (th + y) * log(th + mu)))
    link <- substitute(link)
    fam0 <- if (missing(init.theta)) 
        do.call("poisson", list(link = link))
    else do.call("negative.binomial", list(theta = init.theta, 
        link = link))
    mf <- Call <- match.call()
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "etastart", "mustart", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval.parent(mf)
    Terms <- attr(mf, "terms")
    if (method == "model.frame") 
        return(mf)
    Y <- model.response(mf, "numeric")
    X <- if (!is.empty.model(Terms)) 
        model.matrix(Terms, mf, contrasts)
    else matrix(, NROW(Y), 0)
    w <- model.weights(mf)
    if (!length(w)) 
        w <- rep(1, nrow(mf))
    else if (any(w < 0)) 
        stop("negative weights not allowed")
    offset <- model.offset(mf)
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")
    n <- length(Y)
    if (!missing(method)) {
        if (!exists(method, mode = "function")) 
            stop(gettextf("unimplemented method: %s", sQuote(method)), 
                domain = NA)
        glm.fitter <- get(method)
    }
    else {
        method <- "glm.fit"
        glm.fitter <- stats::glm.fit
    }
    if (control$trace > 1) 
      message("Initial fit:")

    fit <- glm.fitter(x = X, y = Y, w = w, start = start, etastart = etastart, 
        mustart = mustart, offset = offset, family = fam0, control = list(maxit = control$maxit, 
            epsilon = control$epsilon, trace = control$trace > 
                1), intercept = attr(Terms, "intercept") > 0)
    class(fit) <- c("glm", "lm")
    mu <- fit$fitted.values
    if (is.null(theta))
    {
      th <- as.vector(theta.ml(Y, mu, sum(w), w, limit = control$maxit, 
                               trace = control$trace > 2))
    }
    else
  {
    th = theta

    if (control$trace > 1) 
      message(gettextf("Fixing theta value to 'theta': %f", signif(th)), 
              domain = NA)
  }
      
    if (control$trace > 1) 
        message(gettextf("Initial value for 'theta': %f", signif(th)), 
            domain = NA)
    fam <- do.call("negative.binomial", list(theta = th[1], link = link))
    iter <- 0
    d1 <- sqrt(2 * max(1, fit$df.residual))
    d2 <- del <- 1
    g <- fam$linkfun
    Lm <- loglik(n, th, mu, Y, w)
    Lm0 <- Lm + 2 * d1
    while ((iter <- iter + 1) <= control$maxit && (abs(Lm0 - 
        Lm)/d1 + abs(del)/d2) > control$epsilon) {
          eta <- g(mu)
          fit <- glm.fitter(x = X, y = Y, w = w, etastart = eta, 
            offset = offset, family = fam, control = list(maxit = control$maxit, 
                epsilon = control$epsilon, trace = control$trace > 
                  1), intercept = attr(Terms, "intercept") > 
                0)
        t0 <- th
        if (is.null(theta))
        {
          th <- theta.ml(Y, mu, sum(w), w, limit = control$maxit, 
                         trace = control$trace > 2)
        } else
        {
          th = theta
        }

        fam <- do.call("negative.binomial", list(theta = th[1],  ## we don't need all the thetas here if theta is vectorized
                                                 link = link)) 

        mu <- fit$fitted.values
        del <- t0 - th ## this is where the vectorized theta makes a difference
        Lm0 <- Lm
        Lm <- loglik(n, th, mu, Y, w) ## and here - log likelihood computation
        if (control$trace) {
            Ls <- loglik(n, th, Y, Y, w)
            Dev <- 2 * (Ls - Lm)
            message(sprintf("Theta(%d) = %f, 2(Ls - Lm) = %f", 
                iter, signif(th), signif(Dev)), domain = NA)
        }
    }
    if (!is.null(attr(th, "warn"))) 
        fit$th.warn <- attr(th, "warn")
    if (iter > control$maxit) {
        warning("alternation limit reached")
        fit$th.warn <- gettext("alternation limit reached")
    }
    if (length(offset) && attr(Terms, "intercept")) {
        null.deviance <- if (length(Terms)) 
            glm.fitter(X[, "(Intercept)", drop = FALSE], Y, w, 
                offset = offset, family = fam, control = list(maxit = control$maxit, 
                  epsilon = control$epsilon, trace = control$trace > 
                    1), intercept = TRUE)$deviance
        else fit$deviance
        fit$null.deviance <- null.deviance
    }
    class(fit) <- c("negbin", "glm", "lm")
    fit$terms <- Terms
    fit$formula <- as.vector(attr(Terms, "formula"))
    Call$init.theta <- signif(as.vector(th), 10)
    Call$link <- link
    fit$call <- Call
    if (model) 
        fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    if (x) 
        fit$x <- X
    if (!y) 
        fit$y <- NULL
    fit$theta <- as.vector(th)
    fit$SE.theta <- attr(th, "SE")
    fit$twologlik <- as.vector(2 * Lm)
    fit$aic <- -fit$twologlik + 2 * fit$rank + 2
    fit$contrasts <- attr(X, "contrasts")
    fit$xlevels <- .getXlevels(Terms, mf)
    fit$method <- method
    fit$control <- control
    fit$offset <- offset
    fit
}

#' @title dflm
#' @description
#'
#' Formats lm, glm, or fisher.test outputs into readable data.table
#'
dflm = function(x, last = FALSE, nm = '')
{
  if (is.null(x))
    out = data.frame(name = nm, method = as.character(NA), p = as.numeric(NA), estimate = as.numeric(NA), ci.lower = as.numeric(NA),  ci.upper = as.numeric(NA), effect = as.character(NA))
  else if ('lm' %in% class(x))
  {

    coef = as.data.frame(summary(x)$coefficients)
    colnames(coef) = c('estimate', 'se', 'stat', 'p')
    if (last)
      coef = coef[nrow(coef), ]
    coef$ci.lower = coef$estimate - 1.96*coef$se
    coef$ci.upper = coef$estimate + 1.96*coef$se


    if (summary(x)$family$link %in% c('log', 'logit'))
    {
      coef$estimate = exp(coef$estimate)
      coef$ci.upper= exp(coef$ci.upper)
      coef$ci.lower= exp(coef$ci.lower)
    }
    if (!last)
      nm = paste(nm, rownames(coef))
    out = data.frame(name = nm, method = summary(x)$family$family, p = signif(coef$p, 3),
                     p.right = pnorm(summary(x)$coefficients[,3], lower = FALSE),
                     p.left = pnorm(summary(x)$coefficients[,3], lower = TRUE), estimate = coef$estimate, ci.lower = coef$ci.lower, ci.upper = coef$ci.upper, effect = paste(signif(coef$estimate, 3), ' [',  signif(coef$ci.lower,3),'-', signif(coef$ci.upper, 3), ']', sep = ''))
  }
  else
  {
    out = data.frame(name = nm, method = x$method, p = signif(x$p.value, 3), estimate = x$estimate, ci.lower = x$conf.int[1], ci.upper = x$conf.int[2], effect = paste(signif(x$estimate, 3), ' [',  signif(x$conf.int[1],3),'-', signif(x$conf.int[2], 3), ']', sep = ''))
  }

  out$effect = as.character(out$effect)
  out$name = gsub('\\s+$', '', gsub('^\\s+', '', as.character(out$name)))
  out$method = as.character(out$method)
  out$p.twosided = as.numeric(out$p)
  out$p = out$p.right
  out$p.right = NULL
  rownames(out) = NULL
  out = as.data.table(out)
  setkey(out, 'name')
  return(out)
}



#' @name score
#' @title score 1 or more fishHook models
#' @description
#'
#' Scores a set of K (>1) fishHook models defined over <identical> hypothesis sets.
#' Each model k in K represents a background model over a (disjoin) collection of event sets.
#'
#' In practice, each event set k could represent a different variant types (eg indels, SV, SNVs) that each
#' have a separate fit (captured in model k) with respect to a set of covariates.
#' Each event set k could also represent
#' patient (or patient set) specific background models, e.g. colon cancer vs breast cancer,
#' or a combination of patient set and variant type (e.g. indels in lung adenocarcinoma, SVs in breast
#' cancer).
#'
#' The goal is of score() is to combine all the input models / data and derive a hypothesis
#' specific p value for the mutational enrichment (or depletion).
#'
#' Since each input model k has already computed an expected value e_ik for each hypothesis i,
#' we can integrate these models through a second glm which uses this e_ik as an offset,
#' and computes a hypothesis (or hypothesis set) specific intercept.  The value of
#' this intercept and associated p value will represent the mutational enrichment (or depletion)
#' of that hypothesis interval (or hypothesis interval set) from all the various input
#' datasets.
#' 
#' @param ...  fishHook models with <identical> hypotheses
#' @param sets  named list of integers indexing the hypotheses in the input models
#' @param mc.cores integer number of cores to use
#' @param iter integer max number of iteration of glm.nb to run
#' @param verbose  logical flag of whether to run verbose, default will inherit from underlying fishHook models
#' @param ignore.theta logical flag of whether to ignore.theta and just recompute per hypothesis (not recommended)
#' @return data.table of hypotheses
#' @export
#' @author Marcin Imielinski
score = function(..., sets = NULL, mc.cores = NULL, iter = 200, verbose = NULL, ignore.theta = FALSE)
{
  models = list(...)

  if (length(models)<1)
  {
    stop('Need at least one fishHook model to score, preferably two or more')
  }


  hypotheses = models[[1]]$hypotheses
  nb = models[[1]]$nb

  if (length(models)==1)
  {
    if (is.null(sets))
        warning('Applying alternative (slower) Wald-test test for scoring a single model, for better performance with single model scoring using the $score method')
#      models[[1]]$score
#      return(models[[1]]$res)
    }
  else
  {
    for (j in 2:length(models))
    {
      if (!identical(models[[j]]$hypotheses, hypotheses))
        stop('Model hypotheses must be identical for all input fishHook models')

      if (!identical(models[[j]]$nb, nb))
        stop('Models must all be Poisson or Negative-binomial')
    }
  }

  if (is.null(mc.cores))
    {
      mc.cores = models[[1]]$mc.cores
    }


  if (is.null(verbose))
    {
      verbose = models[[1]]$verbose
    }

  if (is.null(sets))
  {
    sets = split(1:length(models[[1]]$hypotheses), paste0('set_', 1:length(models[[1]]$hypotheses)))
    
    setinfo = cbind(data.table(id = names(sets)), as.data.table(models[[1]]$hypotheses))
    names(sets) = copy(setinfo$id)
    setkey(setinfo, id)
  }
  else
  {
    if (is.null(names(sets)))
    {
      names(sets) = paste0('set_', 1:length(sets))
    }

    names(sets) = dedup(names(sets)) ## make sure sets are all named uniquely
    setinfo = data.table(id = gsub('\\W', '_', paste0('set_', names(sets))), name = names(sets))
    if (any(duplicated(setinfo$id)))
      stop('Set names contain illegal special characters that cannot be resolved, please try again using only alphanumeric characters for set names')
    names(sets) = copy(setinfo$id)
    setkey(setinfo, id)
  }  

  ## collect all results from all models
  res = rbindlist(mclapply(1:length(models), function(i)
  {
    mod = models[[i]]
    if (is.null(mod$res))
    {
      fmessage('Rescoring model ', i)
      mod$score(sets = NULL)
    }
    res = mod$res
    res[, model.id := i]
    res[, query.id := 1:.N]
    res
  }))

  if (verbose)
  {
    fmessage('Scoring ', length(models), ' models with ', length(hypotheses), ' hypotheses and ', length(sets), ' hypothesis sets ')
  }

  if (nb)
    {
      thetas = structure(sapply(models, function(x) x$model$theta), SE.theta = sapply(models, function(x) x$model$SE.theta))
    }


  .score.sets = function(sid, sets, res, nb = nb)
  {
    sets = sets[sid]
    if (verbose>1)
    {
      fmessage('Scoring set(s) ', paste(names(sets), collapse = ', '))
    }

    tmpres = merge(cbind(query.id = unlist(sets), set.id = rep(names(sets), elementNROWS(sets))), res)
    setdata = cbind(data.frame(count = tmpres$count, count.pred = log(tmpres$count.pred)),
                          as.data.frame(as.matrix(Matrix::sparseMatrix(1:nrow(tmpres), match(tmpres$set.id, names(sets)), x = 1,
                                                               dims = c(nrow(tmpres), length(sets)),
                                                               dimnames = list(NULL, names(sets))))))
    
    ## make the formula with covariates
    ## this is a model using the hypothesis specific estimate as an offset and then inferring a
    ## set specific intercept
    setformula = eval(parse(text = paste('count', " ~ ", paste(c('offset(count.pred)', names(sets)), collapse = "+"), '-1'))) 
    if (nb)
    {
      if (!ignore.theta)
        {
          thetas = structure(thetas[tmpres$model.id], SE = attr(thetas, 'SE')[tmpres$model.id])  ## note that each model has a different theta, so we will be taking this into account, via modded glm.nb.fh
        }
      gset = tryCatch(glm.nb.fh(setformula, data = setdata, maxit = iter, theta = thetas),
                      error = function(e) NULL)
    } else
    {
      gset = glm(setformula, data = as.data.frame(setdata), family = poisson)
    }
 
   if (!is.null(gset))
    {
      tmpres = dflm(gset)[names(sets), ]
    }
    else
    {
      tmpres = data.table(name = names(sets), method = NA, p = as.numeric(NA), p.right = as.numeric(NA), p.left = as.numeric(NA), estimate = as.numeric(NA), ci.lower = as.numeric(NA), ci.upper = as.numeric(NA), effect = as.character(NA))

    }

    tmpres[, n := elementNROWS(sets)]
    
    return(tmpres)

  }

  setres = rbindlist(mclapply(1:length(sets), .score.sets, sets = sets, res = res, nb = nb, mc.cores = mc.cores))
  setnames(setres, 'name', 'id')
  setres = merge(setinfo, setres, by = 'id')[, -1, with = FALSE]
  setres$method = gsub('\\(.*$', '', setres$method)
  setres$fdr = signif(p.adjust(setres$p, 'BH'), 2) ## compute q value
  return(setres)
}


#' @name fftab
#' @title Tabulate data in Rle
#' @description
#'
#' Tabulates data in ffTrack file across a set of interavls (GRanges)
#' by counting the number of positions matching a given "signature" or
#' applying FUN to aggregate data.  Returns the input GRanges populated with one or more meta data columns
#' of counts or averages.
#'
#' Similar to gr.val in gUtils
#'
#' ff can be an ffTrack but also an RleList from same genome as intervals.
#'
#' returns a GRanges with additional columns for metadata counts
#'
#' @param ff  ffTrack or RleList to pull data from
#' @param intervals intervals
#' @param signatures Signatures is a named list that specify what is to be tallied.  Each signature (ie list element)
#' consist of an arbitrary length character vector specifying strings to %in% (grep = FALSE)
#' or length 1 character vector to grepl (if grep = TRUE)
#' or a length 1 or 2 numeric vector specifying exact value or interval to match (for numeric data)
#'
#' Every list element of signature will become a metadata column in the output GRanges
#' specifying how many positions in the given interval match the given query
#' @param FUN function to aggregate with (default is sum)
#' @param grep logical flag (default FALSE), if TRUE will treat the strings in signature as inputs to grep (instead of exact matches if FALSE)
#' @param mc.cores how many cores (default 1)
#' @param chunksize chunk of FF to bring into memory (i.e. the width of interval), decrease if memory becomes an issue
#' @param verbose logical flag
#' @param na.rm logical flag whether to remove na during aggregation.
#' @importFrom data.table rbindlist data.table setkey :=
#' @importFrom gUtils gr.sub seg2gr gr.stripstrand si2gr rle.query gr.fix gr.chr gr.tile grl.unlist gr.findoverlaps gr.dice hg_seqlengths
fftab = function(ff, intervals, signatures = NULL, FUN = sum, grep = FALSE, mc.cores = 1, chunksize = 1e6, verbose = TRUE, na.rm = TRUE)
    {

    id = ix = NULL ## NOTE fix
        if (!is(ff, 'ffTrack') & !is(ff, 'RleList'))
            stop('Input ff should be ffTrack or RleList\n')

        if (length(intervals)==0)
            stop('Must provide non empty interavl input as GRanges')

        if (!is.null(signatures))
            {
                if (!is.list(signatures))
                    stop('Signatures must be a named list of arbitrary length character or length 1 or 2 numeric vectors')

                if (is.null(names(signatures)))
                    names(signatures) = paste('sig', 1:length(signatures), sep = '')

                check = sapply(signatures, function(x)
                    {
                        if (is.numeric(x))
                            return(length(x)>=1 & length(x)<=2)
                        if (is.character(x))
                            if (grep)
                                return(length(x)==1)
                        return(TRUE)
                    })

                if (!all(check))
                    stop('signatures input is malformed, should be either length 1 or 2 numeric, length 1 character (if grep = TRUE), or atbitrary length character otherwise)')

            }
        else
            signatures = list(score = numeric()) ## we are just scoring bases


        ## generate command that will be executed at each access
        cmd = paste('list(', paste(sapply(names(signatures), function(x)
            {
                sig = signatures[[x]]
                if (is.numeric(sig))
                    {
                                if (length(sig)==0)
                                    cmd = sprintf('%s = FUN(dat, na.rm = na.rm)', x)
                                else if (length(sig)==1)
                                    cmd = sprintf('%s = FUN(dat == %s, na.rm = na.rm)', x, sig[1])
                                else
                                    cmd = sprintf('%s = FUN(dat > %s & dat< %s, na.rm = na.rm)', x, sig[1], sig[2])
                            }
                else
                    if (grep)
                        cmd = sprintf('%s = FUN(grepl("%s", dat), na.rm = na.rm) ', x, sig[1])
                    else
                        cmd = paste(x, '= FUN(dat %in%',
                            paste('c(', paste("\"", sig, "\"", sep = '', collapse = ','), '), na.rm = na.rm)', sep = ''))
                    }), collapse = ', ', sep = ''), ')', sep = '')

        val = values(intervals)
        intervals$ix = 1:length(intervals)

        if (verbose)
            cat('Made command\n')
        ## sorting will hopefully make data access more efficient
        gr = sort(intervals[, 'ix'])

        if (verbose)
            cat('Sorted intervals\n')
        ## tailor the chunking to the size of the individual segments
        gr$chunk.id = ceiling(cumsum(as.numeric(width(gr)))/chunksize)
        gr$num = 1:length(gr)

        ## get down to business
        chunks = split(1:length(gr), gr$chunk.id)
        if (verbose)
            cat('Split intervals\n')

        out = rbindlist(parallel::mclapply(chunks, function(ix)
            {
                chunk = gr[ix]
                if (verbose)
                    cat(sprintf('Intervals %s to %s of %s, total width %s, starting\n', chunk$num[1], chunk$num[length(chunk)], length(gr), sum(width(chunk))))

                if (is(ff, 'ffTrack'))
                    tmp = data.table(
                        dat = ff[chunk],
                        id = rep(1:length(chunk), width(chunk))
                    )
                else ## also can handle rle data
                    tmp = data.table(
                        dat = as.numeric(rle.query(ff, chunk)),
                        id = rep(1:length(chunk), width(chunk))
                    )

                setkey(tmp, id)

                if (verbose)
                    cat(sprintf('Intervals %s to %s of %s: read in ff data\n', chunk$num[1], chunk$num[length(chunk)], length(gr)))

                tab = tmp[, eval(parse(text=cmd)), keyby = id]

                if (verbose)
                    cat(sprintf('Intervals %s to %s of %s: tabulated\n', chunk$num[1], chunk$num[length(chunk)], length(gr)))

                ix = chunk$ix
                out = tab[1:length(chunk), ]
                out$id = NULL
                out$ix = ix
                if (verbose)
                    cat(sprintf('Intervals %s to %s of %s: FINISHED\n', chunk$num[1], chunk$num[length(chunk)], length(gr)))

                return(out)
            }, mc.cores = mc.cores))

        setkey(out, ix)
        out = as.data.frame(out[list(1:length(intervals)), ])
#        out = as.data.frame(out)[order(out$ix),]
        out$ix = NULL
        values(intervals) = cbind(val, out)
        return(intervals)
    }

