
library(fishHook)

library(testthat)

Sys.setenv(DEFAULT_BSGENOME = 'BSgenome.Hsapiens.UCSC.hg19::Hsapiens')

# Sample Events
events = readRDS('/Users/ebiederstedt/fishHook/data/events.rds')
## events = readRDS('events.rds')

# Sample Targets
targets = readRDS('/Users/ebiederstedt/fishHook/data/targets.rds')
## targets = readRDS('targets.rds')


# Sample Covariate
replication_timing = readRDS('/Users/ebiederstedt/fishHook/data/covariate.rds')
## replication_timing = readRDS('covariate.rds')


# Same Eligible Subset
eligible = readRDS('/Users/ebiederstedt/fishHook/data/eligible.rds')
## eligible  = readRDS('eligible.rds')


# indexed pathways
indexed_pathways = readRDS('/Users/ebiederstedt/fishHook/data/indexed_pathways.rds')
## indexed_pathways = readRDS('indexed_pathways.rds')


segs = readRDS('/Users/ebiederstedt/fishHook/data/jabba_segs_11517.rds')
## segs = readRDS('jabba_segs_11517.rds')



context('unit testing fishhook operations')






foofoo = function(targets, covered = NULL, events = NULL,  mc.cores = 1, na.rm = TRUE, pad = 0, verbose = TRUE, max.slice = 1e3, 
    ff.chunk = 1e6, max.chunk = 1e11, out.path = NULL, covariates = list(), maxpatientpergene = Inf, ptidcol = NULL, weightEvents = FALSE, ...)
{
    if(weightEvents){
        maxpatientpergene = NULL
    }

    if (is.character(targets)){
        if (grepl('\\.rds$', targets[1])){
            targets = readRDS(targets[1])
        }
        else if (grepl('(\\.bed$)', targets[1])){
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
    }else {
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
                    }
                    else{
                        ev2$ID = mcols(events)[,ptidcol][ev2$query.id]
                    }
                            
                    ev2$target.id = ov$query.id[ev2$subject.id]
                    tab = as.data.table(cbind(ev2$ID, ev2$target.id))
                    counts.unique = tab[, dummy :=1][, .(count = sum(dummy)), keyby =.(V1, V2)][, count := pmin(maxpatientpergene, count)][, .(final_count = sum(count)), keyby = V2]                          
                }
                        
            }
            ## assume it is an Rle of event counts along the genome  
            else{
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
                }
                
                else{
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
            }
            else if (cov$type == 'numeric'){
                if (is.character(cov$track)){
                    if (grepl('.rds$', cov$track)){
                        cov$track = readRDS(cov$track)
                    }
                    ## assume it is a UCSC format
                    else{
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
                }
                ## then must be GRanges
                else{
                    if (is.na(cov$field)){

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
            }
            else if (cov$type == 'interval'){

                if (is.character(cov$track)){

                    if (grepl('.rds$', cov$track)){
                        cov$track = readRDS(cov$track)
                    }
                    ## assume it is a UCSC format
                    else{
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
    }

    ovdt = gr2dt(ov)
        
        
    cmd = 'list(coverage = sum(width), ';

    if (!is.null(events)){
        cmd = paste(cmd, 'count = sum(count)', sep = '')
    }else{
        cmd = paste(cmd, 'count = NA', sep = '')  
    }
    
        
    cov.nm = setdiff(names(values(ov)), c('coverage', 'count', 'query.id', 'subject.id'))

    if (length(ov) > 0){

        if (length(cov.nm) > 0){
            cmd = paste(cmd,  ',', paste(cov.nm, '= mean(', cov.nm, ')', sep = '', collapse = ', '), ')',  sep = '')
        }else{
            cmd = paste(cmd, ')',  sep = '')
        }

        ovdta =  ovdt[, eval(parse(text = cmd)), keyby = query.id]
        values(targets) = as(as.data.frame(ovdta[list(1:length(targets)), ]), 'DataFrame')

        ##if(!is.null())
        if(!is.null(maxpatientpergene)){
            targets$count = 0
            targets$count[as.numeric(counts.unique$V2)] = counts.unique$final_count
        }
                
    }else{
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



foofoo(targets)
## Overlapping with covered intervals
## Finished overlapping with covered intervals
## Error in foofoo(targets) : object 'counts.unique' not found




