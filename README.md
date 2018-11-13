# Yet Another Methylation Array Toolkit (YAMAT) - Copy-Number Analysis Package

The package is part of [yamat ecosystem for methylation array data analysis](https://github.com/markgene/yamat). It focuses on copy-number 
analysis.

## Introduction

TBA.

## Appendix: CNV Pipelines in Papers and Bioconductor

I am listing a few copy-number analysis pipelines published in research papers 
or archived in Bioconductor. 
 
* [Conumee](https://www.bioconductor.org/packages/release/bioc/html/conumee.html). 
As described in its document, it "contains a set of processing and plotting 
methods for performing copy-number variation (CNV) analysis using Illumina 
450k or EPIC methylation arrays". The parameters are supposed to be tuned as 
described in the vignette. However, it is unclear how the parameters have been 
tuned and which metrics were used for tuning. Also, it is unclear which 
preprocessing workflow works best with the pipeline.

* [MethylCNV](https://www.tandfonline.com/doi/pdf/10.4161/sysb.25896). 

* [ChAMP CNV pipeline](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r30).

* [CopyNumber450k](https://bioconductor.riken.jp/packages/3.1/bioc/html/CopyNumber450k.html).

* [CopyNumber450kCancer](https://cran.r-project.org/web/packages/CopyNumber450kCancer/index.html).

* [RnBeads CNV pipeline](https://bioconductor.org/packages/release/bioc/vignettes/RnBeads/inst/doc/RnBeads.pdf).

A common step of CNV analyses is the segmentation. I am listing some of the 
segmentation methods:

* [DNAcopy](http://bioconductor.org/packages/release/bioc/html/DNAcopy.html)
* [GLAD](http://bioconductor.org/packages/release/bioc/html/GLAD.html)
* [copynumber](http://bioconductor.org/packages/devel/bioc/html/copynumber.html)


