[![Build Status](https://travis-ci.org/mskilab/fishHook.svg?branch=master)](https://travis-ci.org/mskilab/fishHook)
[![codecov.io](https://img.shields.io/codecov/c/github/mskilab/fishHook.svg)](https://codecov.io/github/mskilab/fishHook?branch=master)

fishHook
======

R package for applying Gamma-Poisson regression to identify statistical enrichment or depletion of somatic mutations in regions after correcting for genomic covariates.

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

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

See Demo:

[fishHook Demo](https://github.com/mskilab/fishHook/blob/master/fishHook_Demo.ipynb)

See Documentation

[R Documentation](https://raw.githubusercontent.com/mskilab/fishHook/master/fishHook.pdf)

Description
-----------

Attributions
------------
> Marcin Imielinski - Assistant Professor, Weill-Cornell Medical College. Core Member, New York Genome Center.



Demo
-----------


## Load the Required Packages


```R
library(fishHook)
```

## Now we will need some data
fishHook utilizes Gamma-Poisson regression to identify frequently-mutated or amplified/deleted regions of the genome from sequencing and microarray data. To do this we need to take a set of genomic targets, and test each one against the hypothesis that they are significantly altered in comparison to the other targets. In this first example we will use genes as our targets and use exome data as the mutational events. Since exome sequencing tends to exhibit strong sequencing bias, we want to include this information in our analysis. To do this we constructed a GRanges called eligible that will indicate the regions that have sufficient coverage. 


```R
setwd("~/git/fishHook/data")
```

## Mutational Events


```R
mutational_events = readRDS("events.rds")
mutational_events
```


    GRanges object with 1985704 ranges and 1 metadata column:
                seqnames                 ranges strand   | patient_code
                   <Rle>              <IRanges>  <Rle>   |  <character>
            [1]       10   [52587953, 52587953]      *   | TCGA-D8-A1J8
            [2]       10   [52595854, 52595854]      *   | TCGA-BH-A0HP
            [3]       10   [52595854, 52595854]      *   | TCGA-BH-A0HP
            [4]       10   [52595937, 52595937]      *   | TCGA-BH-A18P
            [5]       10   [52596055, 52596055]      *   | TCGA-AC-A2FB
            ...      ...                    ...    ... ...          ...
      [1985700]       22 [ 39239549,  39239549]      *   | TCGA-ZN-A9VW
      [1985701]       22 [ 50720356,  50720356]      *   | TCGA-ZN-A9VW
      [1985702]        X [ 83362035,  83362035]      *   | TCGA-ZN-A9VW
      [1985703]        X [110005955, 110005955]      *   | TCGA-ZN-A9VW
      [1985704]        X [153588774, 153588774]      *   | TCGA-ZN-A9VW
      -------
      seqinfo: 23 sequences from an unspecified genome; no seqlengths


## Gene Targets


```R
gene_targets = readRDS("targets.rds")
gene_targets
```


    GRanges object with 19688 ranges and 1 metadata column:
              seqnames                 ranges strand   |    gene_name
                 <Rle>              <IRanges>  <Rle>   |  <character>
          [1]       19   [58856544, 58864865]      -   |         A1BG
          [2]       10   [52559169, 52645435]      -   |         A1CF
          [3]       12   [ 9220260,  9268825]      -   |          A2M
          [4]       12   [ 8975068,  9039597]      +   |        A2ML1
          [5]        1   [33772367, 33786699]      -   |      A3GALT2
          ...      ...                    ...    ... ...          ...
      [19684]        7 [143078173, 143088204]      +   |          ZYX
      [19685]       17 [  3907739,   4046314]      -   |        ZZEF1
      [19686]        1 [ 78028101,  78149104]      -   |         ZZZ3
      [19687]       19 [ 14183348,  14185874]      +   | hsa-mir-1199
      [19688]       19 [ 50003781,  50004614]      +   |  hsa-mir-150
      -------
      seqinfo: 25 sequences from an unspecified genome; no seqlengths


## Eligible


```R
eligible = readRDS("eligible.rds")
eligible
```


    GRanges object with 2329621 ranges and 1 metadata column:
                seqnames                 ranges strand   |     score
                   <Rle>              <IRanges>  <Rle>   | <numeric>
            [1]        1       [861300, 861310]      *   |      0.84
            [2]        1       [861311, 861367]      *   |      0.85
            [3]        1       [861368, 861386]      *   |      0.84
            [4]        1       [861387, 861395]      *   |      0.83
            [5]        1       [866417, 866420]      *   |      0.95
            ...      ...                    ...    ... ...       ...
      [2329617]        X [155239804, 155239804]      *   |      0.88
      [2329618]        X [155239805, 155239805]      *   |      0.86
      [2329619]        X [155239806, 155239806]      *   |      0.85
      [2329620]        X [155239807, 155239807]      *   |      0.83
      [2329621]        X [155239808, 155239808]      *   |      0.82
      -------
      seqinfo: 25 sequences from an unspecified genome


## The FishHook Object
All of the data manipulations are handled by the fishHook object given entered data. You can initialize it as follows. 


```R
fish = FishHook$new(targets = gene_targets, events = mutational_events, eligible = eligible)
fish
```


    Contains 19688 hypotheses.
    Contains 1985704 events to map to hypotheses.
    Will map only eliglble regions.
    Covariates:
    No covariates will be used.
    Targets contains 2 metadata columns
    Current State: Initialized
    



## Points about FishHook Object
The FishHook object will take various states during our analysis. You can access this through state or from the output of the fish object. You can also access any of the provided variables.


```R
fish$state
fish$targets
fish$events
fish$eligible
```


    Initialized



    GRanges object with 19688 ranges and 2 metadata columns:
              seqnames                 ranges strand   |    gene_name      name
                 <Rle>              <IRanges>  <Rle>   |  <character> <integer>
          [1]       19   [58856544, 58864865]      -   |         A1BG         1
          [2]       10   [52559169, 52645435]      -   |         A1CF         2
          [3]       12   [ 9220260,  9268825]      -   |          A2M         3
          [4]       12   [ 8975068,  9039597]      +   |        A2ML1         4
          [5]        1   [33772367, 33786699]      -   |      A3GALT2         5
          ...      ...                    ...    ... ...          ...       ...
      [19684]        7 [143078173, 143088204]      +   |          ZYX     19684
      [19685]       17 [  3907739,   4046314]      -   |        ZZEF1     19685
      [19686]        1 [ 78028101,  78149104]      -   |         ZZZ3     19686
      [19687]       19 [ 14183348,  14185874]      +   | hsa-mir-1199     19687
      [19688]       19 [ 50003781,  50004614]      +   |  hsa-mir-150     19688
      -------
      seqinfo: 25 sequences from an unspecified genome; no seqlengths



    GRanges object with 1985704 ranges and 1 metadata column:
                seqnames                 ranges strand   | patient_code
                   <Rle>              <IRanges>  <Rle>   |  <character>
            [1]       10   [52587953, 52587953]      *   | TCGA-D8-A1J8
            [2]       10   [52595854, 52595854]      *   | TCGA-BH-A0HP
            [3]       10   [52595854, 52595854]      *   | TCGA-BH-A0HP
            [4]       10   [52595937, 52595937]      *   | TCGA-BH-A18P
            [5]       10   [52596055, 52596055]      *   | TCGA-AC-A2FB
            ...      ...                    ...    ... ...          ...
      [1985700]       22 [ 39239549,  39239549]      *   | TCGA-ZN-A9VW
      [1985701]       22 [ 50720356,  50720356]      *   | TCGA-ZN-A9VW
      [1985702]        X [ 83362035,  83362035]      *   | TCGA-ZN-A9VW
      [1985703]        X [110005955, 110005955]      *   | TCGA-ZN-A9VW
      [1985704]        X [153588774, 153588774]      *   | TCGA-ZN-A9VW
      -------
      seqinfo: 23 sequences from an unspecified genome; no seqlengths



    GRanges object with 2329621 ranges and 1 metadata column:
                seqnames                 ranges strand   |     score
                   <Rle>              <IRanges>  <Rle>   | <numeric>
            [1]        1       [861300, 861310]      *   |      0.84
            [2]        1       [861311, 861367]      *   |      0.85
            [3]        1       [861368, 861386]      *   |      0.84
            [4]        1       [861387, 861395]      *   |      0.83
            [5]        1       [866417, 866420]      *   |      0.95
            ...      ...                    ...    ... ...       ...
      [2329617]        X [155239804, 155239804]      *   |      0.88
      [2329618]        X [155239805, 155239805]      *   |      0.86
      [2329619]        X [155239806, 155239806]      *   |      0.85
      [2329620]        X [155239807, 155239807]      *   |      0.83
      [2329621]        X [155239808, 155239808]      *   |      0.82
      -------
      seqinfo: 25 sequences from an unspecified genome


## Annotating The FishHook Object
In order to test each hypothesis against the null, which will be the normal mutational load of the genes. To do this we will need to count how many events fall into each target gene. We call this process annotation and can be done as follows. Note that we use verbose=F so as to limit spam. This process should take from a few seconds up to a minute.


```R
fish$annotate(verbose = F)
```

## Note that the State of our FishHook Object is now "Annotated"
You can access the annotation information with anno


```R
fish
fish$anno
```


    Contains 19688 hypotheses.
    Contains 1985704 events to map to hypotheses.
    Will map only eliglble regions.
    Covariates:
    No covariates will be used.
    Targets contains 2 metadata columns
    Current State: Annotated
    




    GRanges object with 19688 ranges and 3 metadata columns:
              seqnames                 ranges strand   |  query.id  coverage
                 <Rle>              <IRanges>  <Rle>   | <integer> <numeric>
          [1]       19   [58856544, 58864865]      -   |         1      1467
          [2]       10   [52559169, 52645435]      -   |         2      2103
          [3]       12   [ 9220260,  9268825]      -   |         3      3884
          [4]       12   [ 8975068,  9039597]      +   |         4      4588
          [5]        1   [33772367, 33786699]      -   |         5         0
          ...      ...                    ...    ... ...       ...       ...
      [19684]        7 [143078173, 143088204]      +   |     19684      1337
      [19685]       17 [  3907739,   4046314]      -   |     19685      9203
      [19686]        1 [ 78028101,  78149104]      -   |     19686      2899
      [19687]       19 [ 14183348,  14185874]      +   |     19687         0
      [19688]       19 [ 50003781,  50004614]      +   |     19688         0
                         count
                     <numeric>
          [1]               50
          [2]              209
          [3]  269.58461538461
          [4] 302.999999999992
          [5]             <NA>
          ...              ...
      [19684] 68.9999999999806
      [19685] 318.636363636364
      [19686]              121
      [19687]             <NA>
      [19688]             <NA>
      -------
      seqinfo: 25 sequences from an unspecified genome; no seqlengths


## Scoring the Targets
Now that we have determined the mutational burden (count) at each target, we can now need to create a null model and test each of our hypothesize against this model. Note that because we are using the targets as their own controls there is an assumption that a majority of the targets will follow the null hypothesis.


```R
fish$score()
```

    GRanges object with 19688 ranges and 3 metadata columns:
              seqnames                 ranges strand   |  query.id  coverage
                 <Rle>              <IRanges>  <Rle>   | <integer> <numeric>
          [1]       19   [58856544, 58864865]      -   |         1      1467
          [2]       10   [52559169, 52645435]      -   |         2      2103
          [3]       12   [ 9220260,  9268825]      -   |         3      3884
          [4]       12   [ 8975068,  9039597]      +   |         4      4588
          [5]        1   [33772367, 33786699]      -   |         5         0
          ...      ...                    ...    ... ...       ...       ...
      [19684]        7 [143078173, 143088204]      +   |     19684      1337
      [19685]       17 [  3907739,   4046314]      -   |     19685      9203
      [19686]        1 [ 78028101,  78149104]      -   |     19686      2899
      [19687]       19 [ 14183348,  14185874]      +   |     19687         0
      [19688]       19 [ 50003781,  50004614]      +   |     19688         0
                         count
                     <numeric>
          [1]               50
          [2]              209
          [3]  269.58461538461
          [4] 302.999999999992
          [5]             <NA>
          ...              ...
      [19684] 68.9999999999806
      [19685] 318.636363636364
      [19686]              121
      [19687]             <NA>
      [19688]             <NA>
      -------
      seqinfo: 25 sequences from an unspecified genome; no seqlengths


    Loading required package: MASS
    
    Attaching package: ‘MASS’
    
    The following object is masked from ‘package:plotly’:
    
        select
    
    The following object is masked from ‘package:VariantAnnotation’:
    
        select
    


    Setting up problem
    Fitting model with 18,418 data points and 0 covariates
    Scoring results


## Note that the State of our FishHook Object is now "Scored"
You can access the scoring information with scores. Or if you want to merge this with the origanl targets data you can use 'all'. This includes the p and q values assigned to each target.


```R
fish
fish$scores[1:10]
fish$all[1:10]
```


    Contains 19688 hypotheses.
    Contains 1985704 events to map to hypotheses.
    Will map only eliglble regions.
    Covariates:
    No covariates will be used.
    Targets contains 2 metadata columns
    Current State: Scored
    




<table>
<thead><tr><th scope=col>seqnames</th><th scope=col>start</th><th scope=col>end</th><th scope=col>width</th><th scope=col>strand</th><th scope=col>query.id</th><th scope=col>coverage</th><th scope=col>count</th><th scope=col>count.pred</th><th scope=col>count.pred.density</th><th scope=col>count.density</th><th scope=col>p</th><th scope=col>q</th><th scope=col>p.neg</th><th scope=col>q.neg</th><th scope=col>effectsize</th></tr></thead>
<tbody>
	<tr><td>19           </td><td> 58856544    </td><td> 58864865    </td><td> 8322        </td><td>-            </td><td> 1           </td><td>1467         </td><td> 50          </td><td> 80.33959    </td><td>0.05476455   </td><td>0.03408316   </td><td>0.800        </td><td>1.0          </td><td>0.21         </td><td>0.85         </td><td>-0.6841830323</td></tr>
	<tr><td>10           </td><td> 52559169    </td><td> 52645435    </td><td>86267        </td><td>-            </td><td> 2           </td><td>2103         </td><td>209          </td><td>115.16984    </td><td>0.05476455   </td><td>0.09938184   </td><td>0.051        </td><td>0.8          </td><td>0.95         </td><td>1.00         </td><td> 0.8597399309</td></tr>
	<tr><td>12           </td><td>  9220260    </td><td>  9268825    </td><td>48566        </td><td>-            </td><td> 3           </td><td>3884         </td><td>270          </td><td>212.70550    </td><td>0.05476455   </td><td>0.06951596   </td><td>0.240        </td><td>1.0          </td><td>0.76         </td><td>0.94         </td><td> 0.3441020452</td></tr>
	<tr><td>12           </td><td>  8975068    </td><td>  9039597    </td><td>64530        </td><td>+            </td><td> 4           </td><td>4588         </td><td>303          </td><td>251.25975    </td><td>0.05476455   </td><td>0.06604185   </td><td>0.280        </td><td>1.0          </td><td>0.72         </td><td>0.93         </td><td> 0.2701382411</td></tr>
	<tr><td>1            </td><td> 33772367    </td><td> 33786699    </td><td>14333        </td><td>-            </td><td> 5           </td><td>   0         </td><td> NA          </td><td>  0.00000    </td><td>       NaN   </td><td>        NA   </td><td>1.000        </td><td>1.0          </td><td>  NA         </td><td>  NA         </td><td>           NA</td></tr>
	<tr><td>22           </td><td> 43088127    </td><td> 43117304    </td><td>29178        </td><td>-            </td><td> 6           </td><td> 737         </td><td> 36          </td><td> 40.36147    </td><td>0.05476455   </td><td>0.04884668   </td><td>0.540        </td><td>1.0          </td><td>0.47         </td><td>0.85         </td><td>-0.1649818741</td></tr>
	<tr><td>3            </td><td>137842560    </td><td>137851229    </td><td> 8670        </td><td>-            </td><td> 7           </td><td>1168         </td><td> 64          </td><td> 63.96499    </td><td>0.05476455   </td><td>0.05479452   </td><td>0.440        </td><td>1.0          </td><td>0.56         </td><td>0.87         </td><td> 0.0007893747</td></tr>
	<tr><td>12           </td><td> 53701240    </td><td> 53718648    </td><td>17409        </td><td>-            </td><td> 8           </td><td>1838         </td><td> 77          </td><td>100.65724    </td><td>0.05476455   </td><td>0.04189336   </td><td>0.670        </td><td>1.0          </td><td>0.34         </td><td>0.85         </td><td>-0.3865205770</td></tr>
	<tr><td>12           </td><td>125549925    </td><td>125627873    </td><td>77949        </td><td>+            </td><td> 9           </td><td>2155         </td><td>102          </td><td>118.01760    </td><td>0.05476455   </td><td>0.04733179   </td><td>0.570        </td><td>1.0          </td><td>0.43         </td><td>0.85         </td><td>-0.2104328784</td></tr>
	<tr><td>3            </td><td>151531825    </td><td>151546276    </td><td>14452        </td><td>+            </td><td>10           </td><td>1328         </td><td>100          </td><td> 72.72732    </td><td>0.05476455   </td><td>0.07530120   </td><td>0.190        </td><td>1.0          </td><td>0.82         </td><td>0.97         </td><td> 0.4594306920</td></tr>
</tbody>
</table>




<table>
<thead><tr><th scope=col>seqnames</th><th scope=col>start</th><th scope=col>end</th><th scope=col>width</th><th scope=col>strand</th><th scope=col>query.id</th><th scope=col>coverage</th><th scope=col>count</th><th scope=col>count.pred</th><th scope=col>count.pred.density</th><th scope=col>count.density</th><th scope=col>p</th><th scope=col>q</th><th scope=col>p.neg</th><th scope=col>q.neg</th><th scope=col>effectsize</th><th scope=col>gene_name</th><th scope=col>name</th></tr></thead>
<tbody>
	<tr><td>19           </td><td> 58856544    </td><td> 58864865    </td><td> 8322        </td><td>-            </td><td> 1           </td><td>1467         </td><td> 50          </td><td> 80.33959    </td><td>0.05476455   </td><td>0.03408316   </td><td>0.800        </td><td>1.0          </td><td>0.21         </td><td>0.85         </td><td>-0.6841830323</td><td>A1BG         </td><td> 1           </td></tr>
	<tr><td>10           </td><td> 52559169    </td><td> 52645435    </td><td>86267        </td><td>-            </td><td> 2           </td><td>2103         </td><td>209          </td><td>115.16984    </td><td>0.05476455   </td><td>0.09938184   </td><td>0.051        </td><td>0.8          </td><td>0.95         </td><td>1.00         </td><td> 0.8597399309</td><td>A1CF         </td><td> 2           </td></tr>
	<tr><td>12           </td><td>  9220260    </td><td>  9268825    </td><td>48566        </td><td>-            </td><td> 3           </td><td>3884         </td><td>270          </td><td>212.70550    </td><td>0.05476455   </td><td>0.06951596   </td><td>0.240        </td><td>1.0          </td><td>0.76         </td><td>0.94         </td><td> 0.3441020452</td><td>A2M          </td><td> 3           </td></tr>
	<tr><td>12           </td><td>  8975068    </td><td>  9039597    </td><td>64530        </td><td>+            </td><td> 4           </td><td>4588         </td><td>303          </td><td>251.25975    </td><td>0.05476455   </td><td>0.06604185   </td><td>0.280        </td><td>1.0          </td><td>0.72         </td><td>0.93         </td><td> 0.2701382411</td><td>A2ML1        </td><td> 4           </td></tr>
	<tr><td>1            </td><td> 33772367    </td><td> 33786699    </td><td>14333        </td><td>-            </td><td> 5           </td><td>   0         </td><td> NA          </td><td>  0.00000    </td><td>       NaN   </td><td>        NA   </td><td>1.000        </td><td>1.0          </td><td>  NA         </td><td>  NA         </td><td>           NA</td><td>A3GALT2      </td><td> 5           </td></tr>
	<tr><td>22           </td><td> 43088127    </td><td> 43117304    </td><td>29178        </td><td>-            </td><td> 6           </td><td> 737         </td><td> 36          </td><td> 40.36147    </td><td>0.05476455   </td><td>0.04884668   </td><td>0.540        </td><td>1.0          </td><td>0.47         </td><td>0.85         </td><td>-0.1649818741</td><td>A4GALT       </td><td> 6           </td></tr>
	<tr><td>3            </td><td>137842560    </td><td>137851229    </td><td> 8670        </td><td>-            </td><td> 7           </td><td>1168         </td><td> 64          </td><td> 63.96499    </td><td>0.05476455   </td><td>0.05479452   </td><td>0.440        </td><td>1.0          </td><td>0.56         </td><td>0.87         </td><td> 0.0007893747</td><td>A4GNT        </td><td> 7           </td></tr>
	<tr><td>12           </td><td> 53701240    </td><td> 53718648    </td><td>17409        </td><td>-            </td><td> 8           </td><td>1838         </td><td> 77          </td><td>100.65724    </td><td>0.05476455   </td><td>0.04189336   </td><td>0.670        </td><td>1.0          </td><td>0.34         </td><td>0.85         </td><td>-0.3865205770</td><td>AAAS         </td><td> 8           </td></tr>
	<tr><td>12           </td><td>125549925    </td><td>125627873    </td><td>77949        </td><td>+            </td><td> 9           </td><td>2155         </td><td>102          </td><td>118.01760    </td><td>0.05476455   </td><td>0.04733179   </td><td>0.570        </td><td>1.0          </td><td>0.43         </td><td>0.85         </td><td>-0.2104328784</td><td>AACS         </td><td> 9           </td></tr>
	<tr><td>3            </td><td>151531825    </td><td>151546276    </td><td>14452        </td><td>+            </td><td>10           </td><td>1328         </td><td>100          </td><td> 72.72732    </td><td>0.05476455   </td><td>0.07530120   </td><td>0.190        </td><td>1.0          </td><td>0.82         </td><td>0.97         </td><td> 0.4594306920</td><td>AADAC        </td><td>10           </td></tr>
</tbody>
</table>



## Visualizing The Data
Grabbing the raw data from the scores field in the fish object is an ok way to manually go through the data but if we are looking to easily identify what is and what is not significant we would have a hard time with the manual inspection. To solve this issue we can utilize a qqplot that will plot the observed distribution of p values versus the expected (uniform) distribution of p values. Significnat hits will be ones that vary greatly from the expected.


```R
plot <- fish$qq_plot(plotly = F)

```



![](images/standard_plot_noplotly.png)

## Visualizing the Data cont.
The above is cool and all but we probably want to annotate the hover text of each point with targets metadata, to do that we can use the columns param in qq_plot(). Note that you can specify any column that is present in the 'all' output. You can also provide your own vectors through annotations. P value will be included in all graphs created but Count, Effectsize, HypothesisID and q will only be added by default if not annotations are specified by the user.


```R
fish$all[1:10]

"Column Annotations"
plot1 <- fish$qq_plot(columns = c("gene_name"))
plot1

"Novel Annotations"
plot2 <- fish$qq_plot(columns = c("gene_name"), annotations = list(test = c("testing", "123")))
plot2
```

<table>
<thead><tr><th scope=col>seqnames</th><th scope=col>start</th><th scope=col>end</th><th scope=col>width</th><th scope=col>strand</th><th scope=col>query.id</th><th scope=col>coverage</th><th scope=col>count</th><th scope=col>count.pred</th><th scope=col>count.pred.density</th><th scope=col>count.density</th><th scope=col>p</th><th scope=col>q</th><th scope=col>p.neg</th><th scope=col>q.neg</th><th scope=col>effectsize</th><th scope=col>gene_name</th><th scope=col>name</th></tr></thead>
<tbody>
	<tr><td>19           </td><td> 58856544    </td><td> 58864865    </td><td> 8322        </td><td>-            </td><td> 1           </td><td>1467         </td><td> 50          </td><td> 80.33959    </td><td>0.05476455   </td><td>0.03408316   </td><td>0.800        </td><td>1.0          </td><td>0.21         </td><td>0.85         </td><td>-0.6841830323</td><td>A1BG         </td><td> 1           </td></tr>
	<tr><td>10           </td><td> 52559169    </td><td> 52645435    </td><td>86267        </td><td>-            </td><td> 2           </td><td>2103         </td><td>209          </td><td>115.16984    </td><td>0.05476455   </td><td>0.09938184   </td><td>0.051        </td><td>0.8          </td><td>0.95         </td><td>1.00         </td><td> 0.8597399309</td><td>A1CF         </td><td> 2           </td></tr>
	<tr><td>12           </td><td>  9220260    </td><td>  9268825    </td><td>48566        </td><td>-            </td><td> 3           </td><td>3884         </td><td>270          </td><td>212.70550    </td><td>0.05476455   </td><td>0.06951596   </td><td>0.240        </td><td>1.0          </td><td>0.76         </td><td>0.94         </td><td> 0.3441020452</td><td>A2M          </td><td> 3           </td></tr>
	<tr><td>12           </td><td>  8975068    </td><td>  9039597    </td><td>64530        </td><td>+            </td><td> 4           </td><td>4588         </td><td>303          </td><td>251.25975    </td><td>0.05476455   </td><td>0.06604185   </td><td>0.280        </td><td>1.0          </td><td>0.72         </td><td>0.93         </td><td> 0.2701382411</td><td>A2ML1        </td><td> 4           </td></tr>
	<tr><td>1            </td><td> 33772367    </td><td> 33786699    </td><td>14333        </td><td>-            </td><td> 5           </td><td>   0         </td><td> NA          </td><td>  0.00000    </td><td>       NaN   </td><td>        NA   </td><td>1.000        </td><td>1.0          </td><td>  NA         </td><td>  NA         </td><td>           NA</td><td>A3GALT2      </td><td> 5           </td></tr>
	<tr><td>22           </td><td> 43088127    </td><td> 43117304    </td><td>29178        </td><td>-            </td><td> 6           </td><td> 737         </td><td> 36          </td><td> 40.36147    </td><td>0.05476455   </td><td>0.04884668   </td><td>0.540        </td><td>1.0          </td><td>0.47         </td><td>0.85         </td><td>-0.1649818741</td><td>A4GALT       </td><td> 6           </td></tr>
	<tr><td>3            </td><td>137842560    </td><td>137851229    </td><td> 8670        </td><td>-            </td><td> 7           </td><td>1168         </td><td> 64          </td><td> 63.96499    </td><td>0.05476455   </td><td>0.05479452   </td><td>0.440        </td><td>1.0          </td><td>0.56         </td><td>0.87         </td><td> 0.0007893747</td><td>A4GNT        </td><td> 7           </td></tr>
	<tr><td>12           </td><td> 53701240    </td><td> 53718648    </td><td>17409        </td><td>-            </td><td> 8           </td><td>1838         </td><td> 77          </td><td>100.65724    </td><td>0.05476455   </td><td>0.04189336   </td><td>0.670        </td><td>1.0          </td><td>0.34         </td><td>0.85         </td><td>-0.3865205770</td><td>AAAS         </td><td> 8           </td></tr>
	<tr><td>12           </td><td>125549925    </td><td>125627873    </td><td>77949        </td><td>+            </td><td> 9           </td><td>2155         </td><td>102          </td><td>118.01760    </td><td>0.05476455   </td><td>0.04733179   </td><td>0.570        </td><td>1.0          </td><td>0.43         </td><td>0.85         </td><td>-0.2104328784</td><td>AACS         </td><td> 9           </td></tr>
	<tr><td>3            </td><td>151531825    </td><td>151546276    </td><td>14452        </td><td>+            </td><td>10           </td><td>1328         </td><td>100          </td><td> 72.72732    </td><td>0.05476455   </td><td>0.07530120   </td><td>0.190        </td><td>1.0          </td><td>0.82         </td><td>0.97         </td><td> 0.4594306920</td><td>AADAC        </td><td>10           </td></tr>
</tbody>
</table>




        Column Annotations

![](images/plotly1.png)

        Novel Annotations

![](images/plotly2.png)




## Covariates
Now we know how to test for which targets are a hotspot for mutations. However, mutational hotspots can be caused by various biological phenomina that are unrelated to cancer. Fore example, replication timing, transcription status, chromatin state and sequence context can all play a role in the formation of mutations. We refer to these biological factors that influence mutation covariates. FishHook has its own object for instantiating covariates, but first lets load up the replication timing covariate as a Genomic Ranges object. It contain a 'score' for each region of the genome.


```R
replication_timing = readRDS("covariate.rds")
replication_timing
```


    GRanges object with 2385966 ranges and 1 metadata column:
                seqnames                ranges strand   |     score
                   <Rle>             <IRanges>  <Rle>   | <numeric>
            [1]        1        [10150, 10275]      *   |  0.117972
            [2]        1        [10275, 10501]      *   |  0.117178
            [3]        1        [10501, 14889]      *   |  0.115744
            [4]        1        [14889, 16295]      *   |  0.088761
            [5]        1        [16295, 17403]      *   |  0.079838
            ...      ...                   ...    ... ...       ...
      [2385962]        Y [59030705,  59031157]      *   |  -1.19679
      [2385963]        Y [59031157,  59031510]      *   | -1.198627
      [2385964]        Y [59031510,  59032054]      *   | -1.200062
      [2385965]        Y [59032054,  59032496]      *   | -1.202272
      [2385966]        Y [59032496, 308283117]      *   | -1.204069
      -------
      seqinfo: 25 sequences from an unspecified genome


## Creating Covariates
**The following information is required when creating covariates:**

Covariate(referenced with cvs): This is meat of the object and in this case will be our replication timing object. It can be of class GRanges, character (file path), RleList or ffTrack object. In this case replication timing is a GRanges Object.

Type: There are three covariate types. Numeric, like replication timing where each region gets a numeric value assigned to it. Interval, where we indicate regions that are "marked" with this covariate. For example, H3K9me3. Sequence, which can be something like GC content.

Name: The name you give to this covariate

**Other Parameters that are not always required:**

Field: This is for numeric covariates and is the column name where the 'score' is held. Note that it is set to 'score' by default.

Signature: This is only required if the Covariate you are using is an ffTrack Object, this is similar to field.

Pad: This indicates how much to the left and to the right of the covariate we should consider its influence. e.g. if a covariate was from position 100-150 with pad = 5 we would consider it for positions 95-155.



```R
rept = Cov_Arr$new(cvs = replication_timing, type = 'numeric', name = 'rept')
rept
```


    Covariate Number: 1
    Name: rept
    type: numeric	signature: NA
    field: NA	pad: 0
    na.rm: NA	grep: NA
    Covariate Class: GRanges



## Covariate Manipulations:
Covariates can be operated on as if they were atomic


```R
rep1 = c(rept,rept,rept)
rep1
```


    Covariate Number: 1
    Name: rept
    type: numeric	signature: NA
    field: NA	pad: 0
    na.rm: NA	grep: NA
    Covariate Class: GRanges
    
    Covariate Number: 2
    Name: rept
    type: numeric	signature: NA
    field: NA	pad: 0
    na.rm: NA	grep: NA
    Covariate Class: GRanges
    
    Covariate Number: 3
    Name: rept
    type: numeric	signature: NA
    field: NA	pad: 0
    na.rm: NA	grep: NA
    Covariate Class: GRanges




```R
rep2 = c(rep1,rept)
rep2
```


    Covariate Number: 1
    Name: rept
    type: numeric	signature: NA
    field: NA	pad: 0
    na.rm: NA	grep: NA
    Covariate Class: GRanges
    
    Covariate Number: 2
    Name: rept
    type: numeric	signature: NA
    field: NA	pad: 0
    na.rm: NA	grep: NA
    Covariate Class: GRanges
    
    Covariate Number: 3
    Name: rept
    type: numeric	signature: NA
    field: NA	pad: 0
    na.rm: NA	grep: NA
    Covariate Class: GRanges
    
    Covariate Number: 4
    Name: rept
    type: numeric	signature: NA
    field: NA	pad: 0
    na.rm: NA	grep: NA
    Covariate Class: GRanges




```R
rep3 = rep2[c(1,3)]
rep3
```


    Covariate Number: 1
    Name: rept
    type: numeric	signature: NA
    field: NA	pad: 0
    na.rm: NA	grep: NA
    Covariate Class: GRanges
    
    Covariate Number: 2
    Name: rept
    type: numeric	signature: NA
    field: NA	pad: 0
    na.rm: NA	grep: NA
    Covariate Class: GRanges




```R
rept = rep3[1]
rept
```


    Covariate Number: 1
    Name: rept
    type: numeric	signature: NA
    field: NA	pad: 0
    na.rm: NA	grep: NA
    Covariate Class: GRanges



## Accessing Covariate Fields
Covariate fields such as type are stored as vectors and when you acess the field you will be returned a vector or list in the case of the Covariates themselves that is the same length as your covariates object.


```R
rep3$cvs
```


    [[1]]
    GRanges object with 2385966 ranges and 1 metadata column:
                seqnames                ranges strand   |     score
                   <Rle>             <IRanges>  <Rle>   | <numeric>
            [1]        1        [10150, 10275]      *   |  0.117972
            [2]        1        [10275, 10501]      *   |  0.117178
            [3]        1        [10501, 14889]      *   |  0.115744
            [4]        1        [14889, 16295]      *   |  0.088761
            [5]        1        [16295, 17403]      *   |  0.079838
            ...      ...                   ...    ... ...       ...
      [2385962]        Y [59030705,  59031157]      *   |  -1.19679
      [2385963]        Y [59031157,  59031510]      *   | -1.198627
      [2385964]        Y [59031510,  59032054]      *   | -1.200062
      [2385965]        Y [59032054,  59032496]      *   | -1.202272
      [2385966]        Y [59032496, 308283117]      *   | -1.204069
      -------
      seqinfo: 25 sequences from an unspecified genome
    
    [[2]]
    GRanges object with 2385966 ranges and 1 metadata column:
                seqnames                ranges strand   |     score
                   <Rle>             <IRanges>  <Rle>   | <numeric>
            [1]        1        [10150, 10275]      *   |  0.117972
            [2]        1        [10275, 10501]      *   |  0.117178
            [3]        1        [10501, 14889]      *   |  0.115744
            [4]        1        [14889, 16295]      *   |  0.088761
            [5]        1        [16295, 17403]      *   |  0.079838
            ...      ...                   ...    ... ...       ...
      [2385962]        Y [59030705,  59031157]      *   |  -1.19679
      [2385963]        Y [59031157,  59031510]      *   | -1.198627
      [2385964]        Y [59031510,  59032054]      *   | -1.200062
      [2385965]        Y [59032054,  59032496]      *   | -1.202272
      [2385966]        Y [59032496, 308283117]      *   | -1.204069
      -------
      seqinfo: 25 sequences from an unspecified genome




```R
rep3$type
rep3$signature
```

	'numeric'       'numeric'


	NA          NA



## Multiple Covariates
In the case that you want to create multiple covariates at a given time, you can pass a list of covariate tracks
to the cvs arguement and a vector of correct type to the other arguements.


```R
multi_cov = Cov_Arr$new(cvs = list(replication_timing, replication_timing),
                       name = c('replication1', 'replication2'),
                       type = c('numeric','numeric'), pad = c(0,20),
                       field = c('score', 'score'))
multi_cov
```


    Covariate Number: 1
    Name: replication1
    type: numeric	signature: NA
    field: score	pad: 0
    na.rm: NA	grep: NA
    Covariate Class: GRanges
    
    Covariate Number: 2
    Name: replication2
    type: numeric	signature: NA
    field: score	pad: 20
    na.rm: NA	grep: NA
    Covariate Class: GRanges



## fishHook Analysis using Covariates
The only difference is that when we initiate the class, we will need to pass in the Covariates. Note that annotating the covariates takes some extra time. You can speed this part up by using mc.cores (set number of cores) or with parameters we will cover in the next section.


```R
fish = FishHook$new(targets = gene_targets, events = mutational_events, eligible = eligible, covariates = rept)
fish
fish$annotate(mc.cores = 3,verbose = F)
fish$score()
plot <- fish$qq_plot(columns = c('gene_name','count','q'))
plot
```



    Contains 19688 hypotheses.
    Contains 1985704 events to map to hypotheses.
    Will map only eliglble regions.
    Covariates:
    rept
    Targets contains 2 metadata columns
    Current State: Initialized
    



    GRanges object with 19688 ranges and 4 metadata columns:
              seqnames                 ranges strand   |  query.id  coverage
                 <Rle>              <IRanges>  <Rle>   | <integer> <numeric>
          [1]       19   [58856544, 58864865]      -   |         1      1467
          [2]       10   [52559169, 52645435]      -   |         2      2103
          [3]       12   [ 9220260,  9268825]      -   |         3      3884
          [4]       12   [ 8975068,  9039597]      +   |         4      4588
          [5]        1   [33772367, 33786699]      -   |         5         0
          ...      ...                    ...    ... ...       ...       ...
      [19684]        7 [143078173, 143088204]      +   |     19684      1337
      [19685]       17 [  3907739,   4046314]      -   |     19685      9203
      [19686]        1 [ 78028101,  78149104]      -   |     19686      2899
      [19687]       19 [ 14183348,  14185874]      +   |     19687         0
      [19688]       19 [ 50003781,  50004614]      +   |     19688         0
                         count               rept
                     <numeric>          <numeric>
          [1]               50 -0.297950732144312
          [2]              209 -0.817199707881931
          [3]  269.58461538461  0.951710314812554
          [4] 302.999999999992    1.0355531646779
          [5]             <NA>               <NA>
          ...              ...                ...
      [19684] 68.9999999999806  -1.08295616530581
      [19685] 318.636363636364   1.26355495111961
      [19686]              121   1.34188474219116
      [19687]             <NA>               <NA>
      [19688]             <NA>               <NA>
      -------
      seqinfo: 25 sequences from an unspecified genome; no seqlengths
    Setting up problem
    Fitting model with 18,418 data points and 1 covariates
    Scoring results



![](images/plotly3.png)



## fishHook Analysis using Covariates (cont.)
Covariates rely on our prior knowledge about mutational processes. However, there are likely factors that influence mutations that are not known as thus it would be impossible for us to define a covariate for them. However, all of the mutational evidence is present in the mutational landscape (events) and as such we can create a covariate from our events that we will call local mutational density that can model the mutational landscape in the area surrounding our targets. We can use the flag 'use_local_mut_density' for this. The bin for this covariate be specified using 'local_mut_density_bin' and is by default set to 1e6.


```R
fish = FishHook$new(targets = gene_targets, events = mutational_events, eligible = eligible, covariates = rept,
                   use_local_mut_density = T, local_mut_density_bin = 1e5, verbose = F)

fish

fish$annotate(mc.cores = 3,verbose = F)
fish$score()
plot <- fish$qq_plot(columns = c('gene_name','count','q'))
plot

```


    GRanges object with 30971 ranges and 3 metadata columns:
              seqnames               ranges strand   |  query.id  coverage
                 <Rle>            <IRanges>  <Rle>   | <integer> <numeric>
          [1]        1     [     1, 100000]      +   |         1         0
          [2]        1     [100001, 200000]      +   |         2         0
          [3]        1     [200001, 300000]      +   |         3         0
          [4]        1     [300001, 400000]      +   |         4         0
          [5]        1     [400001, 500000]      +   |         5         0
          ...      ...                  ...    ... ...       ...       ...
      [30967]        Y [59000001, 59100000]      +   |     30967         0
      [30968]        Y [59100001, 59200000]      +   |     30968         0
      [30969]        Y [59200001, 59300000]      +   |     30969         0
      [30970]        Y [59300001, 59373566]      +   |     30970         0
      [30971]        M [       1,    16571]      +   |     30971         0
                  count
              <numeric>
          [1]      <NA>
          [2]      <NA>
          [3]      <NA>
          [4]      <NA>
          [5]      <NA>
          ...       ...
      [30967]      <NA>
      [30968]      <NA>
      [30969]      <NA>
      [30970]      <NA>
      [30971]      <NA>
      -------
      seqinfo: 25 sequences from an unspecified genome
    Setting up problem
    Fitting model with 15,358 data points and 0 covariates
    Scoring results



    Contains 19688 hypotheses.
    Contains 1985704 events to map to hypotheses.
    Will map only eliglble regions.
    Covariates:
    Local Mutation Density
    rept
    Targets contains 2 metadata columns
    Current State: Initialized
    



    GRanges object with 19688 ranges and 5 metadata columns:
              seqnames                 ranges strand   |  query.id  coverage
                 <Rle>              <IRanges>  <Rle>   | <integer> <numeric>
          [1]       19   [58856544, 58864865]      -   |         1      1467
          [2]       10   [52559169, 52645435]      -   |         2      2103
          [3]       12   [ 9220260,  9268825]      -   |         3      3884
          [4]       12   [ 8975068,  9039597]      +   |         4      4588
          [5]        1   [33772367, 33786699]      -   |         5         0
          ...      ...                    ...    ... ...       ...       ...
      [19684]        7 [143078173, 143088204]      +   |     19684      1337
      [19685]       17 [  3907739,   4046314]      -   |     19685      9203
      [19686]        1 [ 78028101,  78149104]      -   |     19686      2899
      [19687]       19 [ 14183348,  14185874]      +   |     19687         0
      [19688]       19 [ 50003781,  50004614]      +   |     19688         0
                         count Local.Mutation.Density               rept
                     <numeric>              <numeric>          <numeric>
          [1]               50      0.044721257913008 -0.297950732144312
          [2]              209     0.0992943081081466 -0.817199707881931
          [3]  269.58461538461     0.0695159629248198  0.951710314812554
          [4] 302.999999999992     0.0542573619593414    1.0355531646779
          [5]             <NA>                   <NA>               <NA>
          ...              ...                    ...                ...
      [19684] 68.9999999999806     0.0735778091208275  -1.08295616530581
      [19685] 318.636363636364      0.034898612890737   1.26355495111961
      [19686]              121      0.042336217552534   1.34188474219116
      [19687]             <NA>                   <NA>               <NA>
      [19688]             <NA>                   <NA>               <NA>
      -------
      seqinfo: 25 sequences from an unspecified genome; no seqlengths



![](images/plotly4.png)


## FishHook Extras: Subsetting
The fishHook object can be subseted in the following way: 

    fish[i,j,k,l] 

where: 
* i is a vector indicating which targets to keep, 
* j is a vector indicating which events to keep,
* k is a vector indicating which covariates to keep, and
* l is a vector indicating which eligible regions to keep

Here are some examples to play with using the previous fish object


```R
fish
test1 = fish[1:10000,1:100000,c(1,5),1:30] 
test1
```


    Contains 19688 hypotheses.
    Contains 1985704 events to map to hypotheses.
    Will map only eliglble regions.
    Covariates:
    Local Mutation Density
    rept
    Targets contains 2 metadata columns
    Current State: Scored
    




    Contains 10000 hypotheses.
    Contains 100000 events to map to hypotheses.
    Will map only eliglble regions.
    Covariates:
    Local Mutation Density
    NA
    Targets contains 2 metadata columns
    Current State: Initialized
    




```R
fish
test2 = fish[1:10000,1:100000,c(1)]                                                                                                                                                                                                                                                                                                                
test2
```


    Contains 19688 hypotheses.
    Contains 1985704 events to map to hypotheses.
    Will map only eliglble regions.
    Covariates:
    Local Mutation Density
    rept
    Targets contains 2 metadata columns
    Current State: Scored
    




    Contains 10000 hypotheses.
    Contains 100000 events to map to hypotheses.
    Will map only eliglble regions.
    Covariates:
    Local Mutation Density
    Targets contains 2 metadata columns
    Current State: Initialized
    




```R
fish
test3 = fish[,1:100000,,1:30]                                                                                                                                                                                                                                                                                                                        
test3
```


    Contains 19688 hypotheses.
    Contains 1985704 events to map to hypotheses.
    Will map only eliglble regions.
    Covariates:
    Local Mutation Density
    rept
    Targets contains 2 metadata columns
    Current State: Scored
    




    Contains 19688 hypotheses.
    Contains 100000 events to map to hypotheses.
    Will map only eliglble regions.
    Covariates:
    Local Mutation Density
    rept
    Targets contains 2 metadata columns
    Current State: Initialized
    




```R
fish
test4 = fish[1:10000]                                                                                                                                                                                                                                                                                                                                
test4
```


    Contains 19688 hypotheses.
    Contains 1985704 events to map to hypotheses.
    Will map only eliglble regions.
    Covariates:
    Local Mutation Density
    rept
    Targets contains 2 metadata columns
    Current State: Scored
    




    Contains 10000 hypotheses.
    Contains 1985704 events to map to hypotheses.
    Will map only eliglble regions.
    Covariates:
    Local Mutation Density
    rept
    Targets contains 2 metadata columns
    Current State: Initialized
    




```R
fish
test5 = fish[,1:100000]  
test5
```


    Contains 19688 hypotheses.
    Contains 1985704 events to map to hypotheses.
    Will map only eliglble regions.
    Covariates:
    Local Mutation Density
    rept
    Targets contains 2 metadata columns
    Current State: Scored
    




    Contains 19688 hypotheses.
    Contains 100000 events to map to hypotheses.
    Will map only eliglble regions.
    Covariates:
    Local Mutation Density
    rept
    Targets contains 2 metadata columns
    Current State: Initialized
    




```R
fish
test6 = fish[,,1]
test6
```


    Contains 19688 hypotheses.
    Contains 1985704 events to map to hypotheses.
    Will map only eliglble regions.
    Covariates:
    Local Mutation Density
    rept
    Targets contains 2 metadata columns
    Current State: Scored
    




    Contains 19688 hypotheses.
    Contains 1985704 events to map to hypotheses.
    Will map only eliglble regions.
    Covariates:
    Local Mutation Density
    Targets contains 2 metadata columns
    Current State: Initialized
    


