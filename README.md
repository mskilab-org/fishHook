[![Build Status](https://travis-ci.org/mskilab/fish.hook.svg?branch=master)](https://travis-ci.org/mskilab/fish.hook)

fish.hook
======

R package for applying Gamma-Poisson regression to identify statistical enrichment or depletion of somatic mutations in regions after correcting for genomic covariates.

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

Installation
-----------

1. Install devtools from CRAN (if you don't have it already)

  ```
  install.packages('devtools')
  ```

2. Load devtools

  ```
  library(devtools)
  ````

3. Install gUtils (if you don't have it already)

  ```
  install_github('mskilab/gUtils')
  ````


4. Install ffTrack

  ```
  install_github('mskilab/ffTrack')
  ````

5. Install fish.hook

  ```
  install_github('mskilab/fish.hook')
  ````

See Demo

Note: Github does not correctly display the rendered version. In order to properly view the rendered version you must run it through ipynb.

[fish.hook Demo](https://github.com/mskilab/fish.hook/blob/master/Fish.hook.class.demo.ipynb)

See documentation

[R Documentation](https://raw.githubusercontent.com/mskilab/fish.hook/master/fish.hook.pdf)

Description
-----------

Attributions
------------
> Marcin Imielinski - Assistant Professor, Weill-Cornell Medical College. Core Member, New York Genome Center.



Demo
-----------


```R
library(fishhook)
library(skitools)
```


```R
## Set this to the path for the dir "data" that contains the demo data
setwd("~/git/fish.hook/data")
# Sample Events
events = readRDS("events.rds")
# Sample Targets
targets = readRDS("targets.rds")
# Sample Covariate
replication_timing = readRDS("covariate.rds")
# Sample Eligible Subset
eligible = readRDS('eligible.rds')
```


```R
## Check the loaded data to make sure it exists
## events contains the 1,988,162 mutation calls of 8475 individual cancers.
events
```


```R
## Targets contains 19688 annotated human genes
## With meta data columns, gene_name, links, summary
## gene_name is the name of the gene
## links is a link to the gene_cards page for the gene
## This will be used later when creating an interactive plot
## summary contains a breif gene_cards summary of the gene
targets
```


```R
## replcation_timing is a GRanges object that we will use as a Covariate
## replication_timing contains 2,385,966 rows
## The meta data columns are : "score"
## score is a numeric Covariate.
replication_timing
```


```R
## eligible contains the regions of the genome that we will whitelist for fish.hook
## In our case, the events data comes from exon sequencing
## exon sequencing does not afford equal coverage throughout the genome so we
## will want to exclude all regions of poor coverage from the analysis
## we will do this by whitelisting only the regions of the genome that have a score of 0.80 or greater
## eligible is a GRanges object with 2,329,621 rows and one meta data column "score"
## Eligible was previously subset to only include regions with a score > 0.80
eligible
```


```R
## First lets create a covariate from replication timing:

## name can be any string, in our case we will use rept 
## Type must be numeric, sequence or interval, in our case the "score" column will be considered a numeric
## For more information on convariate types check the mskilab/Fishhook.R repo.
## For now, try and keep things to GRanges as these seem to be the most stable
rept = Cov$new(Covariate = replication_timing, type = "numeric",name = "rept")

## As an example we can create a vector of covariates to pass to fish.hook
Covariates = c(rept,rept,rept,rept,rept)

## Note that the underlying data structure of Covariates is a list of Covariates
print("Covariate List with rept repeated 5 times:")
Covariates 
```


```R
## This vector is subsetable
Covariates = Covariates[2:4]

print("Subsetted Covariates")
Covariates
```


```R
print("Single Covariate")
Covariates = Covariates[1]

Covariates
```


```R
## Now lets create a "FishHook" object that will hold and manage our data
## Make sure to keep targets, events and eligible as GRanges
## The only required inputs are targets and events
## can modify elements of FishHook object using replace functions.
fish = FishHook$new(targets = targets, events = events, eligible = eligible, covariates = Covariates)

##Lets take a look at our FishHook object:
fish

```


```R
## annotate targets with the events data only in the eligible regions
## Note that if you are testing and trying to get this to run,
## Covariates add about 2 min to the run time (over a base runtime of ~ 20sec) so get rid of them if debugging
##anno = fish$annotateTargets()

## Note that the default annotate run mode counts each event-target intersection as 1 hit
## If you would like to weigh the hits such that a large hit of which only 10% spans a target contributes a count
## of 0.1, you can set the param weightEvents=TRUE
##anno = fish$annotateTargets(weightEvents = TRUE)

## You can also choose to annotate such that each patient can contribute at most "n" events to each target.
## This is to fix some issuses that can arise when you have a few patients with a large amount of mutations 
## in a given target.
## To do this set maxPtGene = n
## and include a column named ID in your events or specify the column name using the param "PtIDCol".
## That indicates which sample/patient this event comes from
anno = fish$annotateTargets(maxPtGene = 1, verbose = FALSE, PtIDCol = "patient_code")

```


```R
## Lets take a look at our Annotate Object
anno

anno$getTargets()
```


```R
## Score Targets, score the annotated data set
## Now we will score the annotated data set to 
## see which of our targets (in this case genes)
## Are significantly mutated
## This portion is quite quick
score = anno$scoreTargets()

## Lets take a look at our scored targets:

score

## We can look at the data store in score by using
score$getAll()
##query.id is a column used during fish.hook
##coverage is the size of the eligible region of the gene
##count is the number of mutations we see at that gene
##rept is our Covariate that we named rept
##p is the p-value
##q is the q-value
```


```R
## Now that we have all of this data we will want a way to visualize it.
## We can both visualize and assess the merit of our analysis by using a qq_plot
## fish.class comes with a built in plotting function built using ploty
## Note that displaying html widget objects in jupyter requires installing "pandoc"
## Can get from : "sudo apt-get install pandoc"

x = score$qq_plot()


x

```


```R
## Here is a nicer version that lets you annotate the points w/metadata
## First lets get all of out Metadata out from the score object
res = score$getAll()

## Lets add hover annotations
## gene_name -> gene names from targets
## q -> q-value from score
## count -> count from annotate
## summary -> gene summary from targets obtained from genecards 
## p-value is added in automatically so no need to add
x = score$qq_plot(annotations = list(Hypothesis_ID = res$gene_name,Count = res$count, q = res$q,Summary = res$summary))

x
```


```R
## you can also plot so that there are no hover annotations other than p-value with:

x = score$qq_plot(annotations = list())

x
```


```R
## You can also set annotations by specifying the columns of your score object to annotate with

x = score$qq_plot(columns = c("count","q","gene_name"))

x

```


```R
## Now lets look at the aggregation function and how we can use it to study pathways.
## First we will need to get some pathways metadata.
pathways = readRDS("indexed_pathways.rds")

##lets take a look at a few pathways:
pathways[1:3]

#And the number of pathways
length(pathways)
```


```R
##Lets take a look at a series of pathways, for example lets look at the KEGG pathways:

Chosen_Pathways_Index = which(grepl("KEGG",names(pathways)))
Chosen_Pathways = pathways[Chosen_Pathways_Index]

## Lets look at our chosen pathways:
Chosen_Pathways[1:10]
length(Chosen_Pathways)
```


```R
## Now lets aggregate our results by pathways:
anno2 = anno$aggregateTargets(by = Chosen_Pathways)

```


```R
##Lets look at our aggregates:
anno2
```


```R
##Scoring:
score = anno2$scoreTargets()
```


```R
##Plotting
x = score$qq_plot()


x
```


```R
## change annotations

new_meta = as.data.frame(names(Chosen_Pathways))

score$replaceMeta(new_meta)

score$getMeta()

x = score$qq_plot()

x
```


```R
targets[Chosen_Pathways$KEGG_LYSINE_DEGRADATION]$gene_name

```


```R
## Another way to call the plotting fuction is to use the names of the columns:

class(score$getAll())
#score$getAll()

x = score$qq_plot(columns = c("count","q","names(Chosen_Pathways)"))

x
```


```R
##Since the names(Chosen_Pathways) looks a bit funky we can use a mix of the two methods

colnames(new_meta) = "Pathways"
x = score$qq_plot(columns = c("count","q"),annotations = list(new_meta))

x
```
