
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
