[![Build Status](https://travis-ci.org/mskilab/fishHook.svg?branch=master)](https://travis-ci.org/mskilab/fishHook)
[![codecov.io](https://img.shields.io/codecov/c/github/mskilab/fishHook.svg)](https://codecov.io/github/mskilab/fishHook?branch=master)


<img src="images/fishhook.png" width = "500">

fishHook
======

R package for applying Gamma-Poisson regression to identify statistical
enrichment or depletion of somatic mutations in arbitrary (sets of) genomic
intervals after correcting for genomic  covariates, e.g. replication timing,
sequence context, chromatin state.

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

# Table of Contents:
-----------
1. [Installation](#installation)
2. [R Documentation](#rdocs)
3. [Attributions](#attributions)
4. [TL;DR](#tldr)
5. [FishHook Demo](#demo)
6. [FishHook Operations](#fishhook_ops)
7. [Active Bindings and FishHook Variables](#active_bindings)
7. [FishHook Functions](#functions)

<div id="installation"/>

Installation 
-----------

1. Install devtools from CRAN (if you don't have it already)

  ```
  install.packages('devtools')
  ```

2. Install gUtils 

  ```
  devtools::install_github('mskilab/gUtils')
  ````


3. Install ffTrack

  ```
  devtools::install_github('mskilab/ffTrack')
  ````

4. Install fishHook

  ```
  devtools::install_github('mskilab/fishHook')
  ````


<p align="center">
<img src="images/qqdreams.jpg" width = "700">
</p>

<div id="rdocs"/>

See Documentation 
------------
[R Documentation](https://raw.githubusercontent.com/mskilab/fishHook/master/fishHook.pdf)



<div id="attributions"/>

Attributions
------------
> Marcin Imielinski - Assistant Professor, Weill-Cornell Medical College. Core Member, New York Genome Center.
> Zoran Gajic - Undergraduate Research Assistant, Rutgers University, Weill Cornell Medicine, New York Genome Center.


<div id="tldr"/>


TL;DR
-----------


### Load the fishHook Package


```R
library(fishHook)
```

###  Load in Your Data
**mutational_events** is a GRanges containing mutations (e.g. snvs/indels or SCNAs)

**targets** is a GRanges containing the start and ends of each gene and some metadata.
This can be other 'targets' such as 1kb tiles of the genome

**eligible** is a GRanges indicating which regions of the genome are captured in our assay (whole exome sequencing)
This can be replace with an array/whole genome sequencing specific track.

**covariate** is a GRanges of replication timing across the genome. This can be replaced with other covariates
such as H3K9me3 profile, chromhmm intervals, etc.


```R
data(events)
data(targets)
data(eligible)
data(replication_timing)
```

### Create a Covariate Object


```R
cov = Cov(data = replication_timing, name = 'rept', type = 'numeric', field = 'score')
```

### Instantiate The FishHook object



```R
fish = Fish(hypotheses = gene_targets, events = mutational_events, eligible = eligible, covariates = cov)
fish
```


### Fit a Gamma Poisson Model of Background Mutational Density and Covariates and Score Your Hypotheses


```R
fish$score()
fish
```

### Statistical Validation of Results with QQ-Plots


```R
plot <- fish$qqp(plotly = F)
```

<div id="demo"/>


[fishHook Tutorial](http://htmlpreview.github.io/?http://github.com/mskilab/fishHook/blob/master/docs/tutorial.html)

[fishHook developer reference](http://htmlpreview.github.io/?http://github.com/mskilab/fishHook/blob/master/docs/reference.md)


