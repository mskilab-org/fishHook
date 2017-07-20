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




```R
library(fishhook)
library(skitools)
```

    Loading required package: GenomicRanges
    Loading required package: BiocGenerics
    Loading required package: parallel
    
    Attaching package: ‘BiocGenerics’
    
    The following objects are masked from ‘package:parallel’:
    
        clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
        clusterExport, clusterMap, parApply, parCapply, parLapply,
        parLapplyLB, parRapply, parSapply, parSapplyLB
    
    The following objects are masked from ‘package:stats’:
    
        IQR, mad, xtabs
    
    The following objects are masked from ‘package:base’:
    
        anyDuplicated, append, as.data.frame, as.vector, cbind, colnames,
        do.call, duplicated, eval, evalq, Filter, Find, get, grep, grepl,
        intersect, is.unsorted, lapply, lengths, Map, mapply, match, mget,
        order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
        rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
        union, unique, unlist, unsplit
    
    Loading required package: S4Vectors
    Loading required package: stats4
    Loading required package: IRanges
    Loading required package: GenomeInfoDb
    Loading required package: gUtils
    Loading required package: data.table
    
    Attaching package: ‘data.table’
    
    The following object is masked from ‘package:GenomicRanges’:
    
        shift
    
    The following object is masked from ‘package:IRanges’:
    
        shift
    
    
    Attaching package: ‘gUtils’
    
    The following object is masked from ‘package:base’:
    
        %o%
    
    Loading required package: VariantAnnotation
    Loading required package: SummarizedExperiment
    Loading required package: Biobase
    Welcome to Bioconductor
    
        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.
    
    Loading required package: Rsamtools
    Loading required package: XVector
    Loading required package: Biostrings
    
    Attaching package: ‘VariantAnnotation’
    
    The following object is masked from ‘package:base’:
    
        tabulate
    
    Loading required package: htmlwidgets
    Loading required package: devtools
    Loading required package: stringr
    
    Attaching package: ‘stringr’
    
    The following object is masked from ‘package:VariantAnnotation’:
    
        fixed
    
    Loading required package: plotly
    Loading required package: ggplot2
    
    Attaching package: ‘ggplot2’
    
    The following object is masked from ‘package:gUtils’:
    
        %+%
    
    
    Attaching package: ‘plotly’
    
    The following object is masked from ‘package:ggplot2’:
    
        last_plot
    
    The following object is masked from ‘package:VariantAnnotation’:
    
        select
    
    The following object is masked from ‘package:XVector’:
    
        slice
    
    The following object is masked from ‘package:IRanges’:
    
        slice
    
    The following object is masked from ‘package:S4Vectors’:
    
        rename
    
    The following object is masked from ‘package:stats’:
    
        filter
    
    The following object is masked from ‘package:graphics’:
    
        layout
    
    Loading required package: reshape2
    
    Attaching package: ‘reshape2’
    
    The following objects are masked from ‘package:data.table’:
    
        dcast, melt
    
    Warning message:
    “replacing previous import by ‘plotly::select’ when loading ‘skitools’”Warning message:
    “replacing previous import by ‘plotly::last_plot’ when loading ‘skitools’”
    Attaching package: ‘skitools’
    
    The following object is masked from ‘package:ggplot2’:
    
        alpha
    
    The following object is masked from ‘package:fishhook’:
    
        qq_pval
    
    The following object is masked from ‘package:gUtils’:
    
        standardize_segs
    
    The following object is masked from ‘package:stats’:
    
        ccf
    
    The following object is masked from ‘package:utils’:
    
        timestamp
    



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



```R
## Targets contains 19688 annotated human genes
## With meta data columns, gene_name, links, summary
## gene_name is the name of the gene
## links is a link to the gene_cards page for the gene
## This will be used later when creating an interactive plot
## summary contains a breif gene_cards summary of the gene
targets
```


    GRanges object with 19688 ranges and 6 metadata columns:
              seqnames                 ranges strand   |      chr      pos1
                 <Rle>              <IRanges>  <Rle>   | <factor> <integer>
          [1]       19   [58856544, 58864865]      -   |       19  58856544
          [2]       10   [52559169, 52645435]      -   |       10  52559169
          [3]       12   [ 9220260,  9268825]      -   |       12   9220260
          [4]       12   [ 8975068,  9039597]      +   |       12   8975068
          [5]        1   [33772367, 33786699]      -   |        1  33772367
          ...      ...                    ...    ... ...      ...       ...
      [19684]        7 [143078173, 143088204]      +   |        7 143078173
      [19685]       17 [  3907739,   4046314]      -   |       17   3907739
      [19686]        1 [ 78028101,  78149104]      -   |        1  78028101
      [19687]       19 [ 14183348,  14185874]      +   |       19  14183348
      [19688]       19 [ 50003781,  50004614]      +   |       19  50003781
                   pos2    gene_name
              <integer>  <character>
          [1]  58864865         A1BG
          [2]  52645435         A1CF
          [3]   9268825          A2M
          [4]   9039597        A2ML1
          [5]  33786699      A3GALT2
          ...       ...          ...
      [19684] 143088204          ZYX
      [19685]   4046314        ZZEF1
      [19686]  78149104         ZZZ3
      [19687]  14185874 hsa-mir-1199
      [19688]  50004614  hsa-mir-150
                                                                       links
                                                                    <factor>
          [1]         http://www.genecards.org/cgi-bin/carddisp.pl?gene=A1BG
          [2]         http://www.genecards.org/cgi-bin/carddisp.pl?gene=A1CF
          [3]          http://www.genecards.org/cgi-bin/carddisp.pl?gene=A2M
          [4]        http://www.genecards.org/cgi-bin/carddisp.pl?gene=A2ML1
          [5]      http://www.genecards.org/cgi-bin/carddisp.pl?gene=A3GALT2
          ...                                                            ...
      [19684]          http://www.genecards.org/cgi-bin/carddisp.pl?gene=ZYX
      [19685]        http://www.genecards.org/cgi-bin/carddisp.pl?gene=ZZEF1
      [19686]         http://www.genecards.org/cgi-bin/carddisp.pl?gene=ZZZ3
      [19687] http://www.genecards.org/cgi-bin/carddisp.pl?gene=hsa-mir-1199
      [19688]  http://www.genecards.org/cgi-bin/carddisp.pl?gene=hsa-mir-150
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     summary
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    <factor>
          [1]                                                                                                                                                                                                                                                                                                           The protein encoded by this gene is a plasma<br>glycoprotein of unknown function. The protein<br>shows sequence similarity to the variable regions<br>of some immunoglobulin supergene family member<br>proteins. [provided by RefSeq, Jul 2008]<br>
          [2]         Mammalian apolipoprotein B mRNA undergoes<br>site-specific C to U deamination, which is<br>mediated by a multi-component enzyme complex<br>containing a minimal core composed of APOBEC-1 and<br>a complementation factor encoded by this gene. The<br>gene product has three non-identical RNA<br>recognition motifs and belongs to the hnRNP R<br>family of RNA-binding proteins. It has been<br>proposed that this complementation factor<br>functions as an RNA-binding subunit and docks<br>APOBEC-1 to deaminate the upstream cytidine. Stud<br>
          [3]                                                                                                                                                                                         Alpha-2-macroglobulin is a protease inhibitor and<br>cytokine transporter. It inhibits many proteases,<br>including trypsin, thrombin and collagenase. A2M<br>is implicated in Alzheimer disease (AD) due to its<br>ability to mediate the clearance and degradation<br>of A-beta, the major component of beta-amyloid<br>deposits. [provided by RefSeq, Jul 2008]<br>
          [4]      This gene encodes a member of the<br>alpha-macroglobulin superfamily. The encoded<br>protein is thought to be an N-glycosylated<br>monomeric protein that acts as an inhibitor of<br>several proteases. It has been shown to form<br>covalent interactions with proteases, and has been<br>reported as the p170 antigen recognized by<br>autoantibodies in the autoimmune disease<br>paraneoplastic pemphigus (PNP; PMID:20805888).<br>Mutations in these gene have also been associated<br>with some cases of Noonan syndrome (NS;<br>PMID:24939586)<br>
          [5] A3GALT2 (Alpha 1,3-Galactosyltransferase 2) is a<br>Protein Coding gene.\r\n                            <br>               \r\n                                 <br>          Among its related pathways are<br>Glycosphingolipid biosynthesis - ganglio series.\r\n<br>                                           GO<br>annotations related to this gene include<br>transferase activity, transferring glycosyl groups<br>and alpha-1,3-galactosyltransferase activity.\r\n   <br>                                        An<br>important paralog of this<br>
          ...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ...
      [19684]         Focal adhesions are actin-rich structures that<br>enable cells to adhere to the extracellular matrix<br>and at which protein complexes involved in signal<br>transduction assemble. Zyxin is a zinc-binding<br>phosphoprotein that concentrates at focal<br>adhesions and along the actin cytoskeleton. Zyxin<br>has an N-terminal proline-rich domain and three<br>LIM domains in its C-terminal half. The<br>proline-rich domain may interact with SH3 domains<br>of proteins involved in signal transduction<br>pathways while the LIM domains <br>
      [19685]                                                                                                                                       ZZEF1 (Zinc Finger ZZ-Type And EF-Hand Domain<br>Containing 1) is a Protein Coding gene.\r\n         <br>                                  \r\n              <br>                             \r\n                   <br>                        GO annotations related to<br>this gene include calcium ion binding.\r\n          <br>                                 An important<br>paralog of this gene is CUL9.<br>
      [19686]                                                                                                                                                                                                    ZZZ3 (Zinc Finger ZZ-Type Containing 3) is a<br>Protein Coding gene.\r\n                            <br>               \r\n                                 <br>          Among its related pathways are Chromatin<br>organization.\r\n                                   <br>        GO annotations related to this gene<br>include chromatin binding.<br>
      [19687]         microRNAs (miRNAs) are short (20-24 nt) non-coding<br>RNAs that are involved in post-transcriptional<br>regulation of gene expression in multicellular<br>organisms by affecting both the stability and<br>translation of mRNAs. miRNAs are transcribed by<br>RNA polymerase II as part of capped and<br>polyadenylated primary transcripts (pri-miRNAs)<br>that can be either protein-coding or non-coding.<br>The primary transcript is cleaved by the Drosha<br>ribonuclease III enzyme to produce an<br>approximately 70-nt stem-loop precurso<br>
      [19688]         microRNAs (miRNAs) are short (20-24 nt) non-coding<br>RNAs that are involved in post-transcriptional<br>regulation of gene expression in multicellular<br>organisms by affecting both the stability and<br>translation of mRNAs. miRNAs are transcribed by<br>RNA polymerase II as part of capped and<br>polyadenylated primary transcripts (pri-miRNAs)<br>that can be either protein-coding or non-coding.<br>The primary transcript is cleaved by the Drosha<br>ribonuclease III enzyme to produce an<br>approximately 70-nt stem-loop precurso<br>
      -------
      seqinfo: 25 sequences from an unspecified genome; no seqlengths



```R
## replcation_timing is a GRanges object that we will use as a Covariate
## replication_timing contains 2,385,966 rows
## The meta data columns are : "score"
## score is a numeric Covariate.
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

    [1] "Covariate List with rept repeated 5 times:"



    $rept
    type:  numeric 	signature:  
    field:  	pad:  
    na.rm:  	grep:  
    track:
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
    
    $rept
    type:  numeric 	signature:  
    field:  	pad:  
    na.rm:  	grep:  
    track:
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
    
    $rept
    type:  numeric 	signature:  
    field:  	pad:  
    na.rm:  	grep:  
    track:
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
    
    $rept
    type:  numeric 	signature:  
    field:  	pad:  
    na.rm:  	grep:  
    track:
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
    
    $rept
    type:  numeric 	signature:  
    field:  	pad:  
    na.rm:  	grep:  
    track:
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



