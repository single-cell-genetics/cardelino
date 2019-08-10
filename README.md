# cardelino: clone identification from single-cell data 

[![Linux Build Status](https://travis-ci.org/PMBio/cardelino.svg?branch=master)](https://travis-ci.org/PMBio/cardelino)
[![codecov.io](https://codecov.io/github/PMBio/cardelino/coverage.svg?branch=master)](https://codecov.io/github/PMBio/cardelino/?branch=master)

<img src=inst/cardelino_sticker.png height="200">

This R package contains a Bayesian method to infer clonal structure for a 
population of cells using single-cell RNA-seq data (and possibly other data 
modalities). 

## Installation

### From R

The **latest** `cardelino` package can be conveniently installed using the 
[`devtools`](https://www.rstudio.com/products/rpackages/devtools/) package thus:

```{R}
devtools::install_github("PMBio/cardelino", build_vignettes = TRUE)
```

### Cardelino in a container

For situations in which installing software is a challenge (for example, on 
institutional HPC clusters or on cloud computing platforms), we provide a 
pre-built [Docker image](https://hub.docker.com/r/davismcc/r-singlecell-img) on
[DockerHub](https://hub.docker.com/). This image contains R version 3.5.0 with 
`cardelino` and other packages (e.g. tidyverse, basic Bioconductor and other
single-cell RNA-seq packages) installed and ready to use with 
[Docker](https://www.docker.com/) or [Singularity](https://www.sylabs.io/).

For example, to build a Singularity image that can be used on an HPC cluster
(with Singularity installed) one simply pulls the image from DockerHub:

```{bash}
singularity build rsc.img docker://davismcc/r-singlecell-img
```

This builds a Singularity image called `rsc.img` in the current working 
directory. We can then run R from the container and use the installed version
of `cardelino`:

```{bash}
singularity exec rsc.img R
```

Equivalent commands enable running R from the container with Docker.

## Getting started

The best place to start are the vignettes. From inside an R session, load 
`cardelino` and then browse the vignettes:

```{r}
library(cardelino)
browseVignettes("cardelino")
```

Vignettes for clone identification use cases are provided. 

Accessing the vignettes from within your R session is recommended, but
you can also [view the clone ID vignette](https://rawgit.com/PMBio/cardelino/master/inst/doc/vignette-cloneid.html).

## Notes for donor deconvolution
The denor demultiplex function, namely Vireo, was supported in this R package 
before, but now has been re-implemented in Python, which is more memory 
efficient and easier to run via a command line. We, therefore, highly recommend 
you switch to the Python version: https://vireoSNP.readthedocs.io

The vireo function is not supported from version >=0.5.0. If you want to use the
R functions, please use the version ==0.4.2 or lower. You can also find it in a
separate branch in this repository: 
[with_vireo branch](https://github.com/PMBio/cardelino/tree/with_vireo)
or use the 
[donor_id.R](https://github.com/PMBio/cardelino/blob/with_vireo/R/donor_id.R) 
file directly.

## About the name

`cardelino` is almost an anagram of "clone ID R" and is almost the same as the 
Italian for "goldfinch", a common and attractive European bird. In the Western 
art canon, the goldfinch is considered a 
["saviour" bird](https://en.wikipedia.org/wiki/European_goldfinch) and appears 
in notable paintings from the 
[Italian renaissance](https://en.wikipedia.org/wiki/Madonna_del_cardellino) and 
the [Dutch Golden Age](https://en.wikipedia.org/wiki/The_Goldfinch_(painting)). 
Perhaps this package may prove a saviour for certain single-cell datasets.

`Vireo`(variational inference for reconstructing ensemebl origin) is a Latin 
word referring to a green migratory bird, perhaps the female golden oriole, 
possibly the [European greenfinch](https://en.wikipedia.org/wiki/European_greenfinch)

<img src=inst/cardelino_med.jpg height="400">

**Acknowledgement:**
The `cardelino` image was produced by [Darren Bellerby](https://www.flickr.com/photos/world-birds/). It was obtained from
[Flickr](https://www.flickr.com/photos/world-birds/18740373165/in/photolist-uy2j3a-uxAdib-aLcHGB-9BjDvc-YkgQg7-QN9Tr1-BVjkHh-8oWiKC-WFkDcS-nhZzXt-Y4zM2h-zULNgX-7uZCFT-f5ghc4-Ugx9pj-UJ5tog-7v4rVy-7wsLpm-bru3Ha-JnmcUQ-frkUqa-bohcgU-KAB14-dieCGY-FJ6n6A-GHJ5UK-X2qjGh-8cAjtw-FshfBi-8cwZst-qEMHSX-dTtAUs-EtqKxo-oZdJB3-8cx1Tn-D1jHjU-PWzWY2-brtKfH-ch2tvW-qEFKTd-wVmxsG-oYZbhP-Aa5cBB-h6aQf6-9Bny23-ayfnFS-dgG2Kn-QUyKgf-bBc31B-cVik3)
and is reproduced here under a CC-BY-2.0
[licence](https://creativecommons.org/licenses/by/2.0/legalcode).


