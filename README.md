# cardelino: clone and donor identification from single-cell data 

[![Linux Build Status](https://travis-ci.org/PMBio/cardelino.svg?branch=master)](https://travis-ci.org/PMBio/cardelino)

<img src=inst/cardelino_sticker.png height="200">

This R package contains methods to assign donor and clone identities to 
individual cells from single-cell RNA-sequencing data.

## Installation

### From R

The `cardelino` package can be conveniently installed using the 
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

Vignettes for the donor identification and clone identification use cases are 
provided. 

Accessing the vignettes from within your R session is recommended, but
you can also [view the clone ID vignette](https://rawgit.com/PMBio/cardelino/master/inst/doc/vignette-cloneid.html) 
and [view the donor ID vignette](https://rawgit.com/PMBio/cardelino/master/inst/doc/vignette-donorid.html) online.



## About the name

`cardelino` is almost an anagram of "clone ID R" and is almost the same as the 
Italian for "goldfinch", a common and attractive European bird. In the Western 
art canon, the goldfinch is considered a 
["saviour" bird](https://en.wikipedia.org/wiki/European_goldfinch) and appears 
in notable paintings from the 
[Italian renaissance](https://en.wikipedia.org/wiki/Madonna_del_cardellino) and 
the [Dutch Golden Age](https://en.wikipedia.org/wiki/The_Goldfinch_(painting)). 
Perhaps this package may prove a saviour for certain single-cell datasets.

<img src=inst/cardelino_med.jpg height="400">

**Acknowledgement:**
The `cardelino` image was produced by [Darren Bellerby](https://www.flickr.com/photos/world-birds/). It was obtained from
[Flickr](https://www.flickr.com/photos/world-birds/18740373165/in/photolist-uy2j3a-uxAdib-aLcHGB-9BjDvc-YkgQg7-QN9Tr1-BVjkHh-8oWiKC-WFkDcS-nhZzXt-Y4zM2h-zULNgX-7uZCFT-f5ghc4-Ugx9pj-UJ5tog-7v4rVy-7wsLpm-bru3Ha-JnmcUQ-frkUqa-bohcgU-KAB14-dieCGY-FJ6n6A-GHJ5UK-X2qjGh-8cAjtw-FshfBi-8cwZst-qEMHSX-dTtAUs-EtqKxo-oZdJB3-8cx1Tn-D1jHjU-PWzWY2-brtKfH-ch2tvW-qEFKTd-wVmxsG-oYZbhP-Aa5cBB-h6aQf6-9Bny23-ayfnFS-dgG2Kn-QUyKgf-bBc31B-cVik3)
and is reproduced here under a CC-BY-2.0
[licence](https://creativecommons.org/licenses/by/2.0/legalcode).


