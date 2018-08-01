## ----knitr-options, echo=FALSE, message=FALSE, warning=FALSE---------------
## To render an HTML version that works nicely with github and web pages, do:
## rmarkdown::render("vignettes/vignette.Rmd", "all")
library(knitr)
opts_chunk$set(fig.align = 'center', fig.width = 6, fig.height = 5, dev = 'png',
    warning = FALSE, error = FALSE, message = FALSE)
library(ggplot2)
library(BiocStyle)
theme_set(theme_bw(12))

## ----load-pkg--------------------------------------------------------------
library(cardelino)
cell_vcf <- read_vcf(system.file("extdata", "cells.donorid.vcf.gz", 
                                 package = "cardelino"))

## ----ppca-raw, fig.height=4, fig.width=6-----------------------------------
raw_geno <-  VariantAnnotation::genotypeToSnpMatrix(
    VariantAnnotation::geno(cell_vcf, "GT"),
    ref = VariantAnnotation::ref(cell_vcf),
    alt = VariantAnnotation::alt(cell_vcf))
pp <- pcaMethods::ppca(as(raw_geno$genotypes, "numeric"))
df <- data.frame(PPCA1 <- pp@scores[, 1], PPCA2 <- pp@scores[, 2])
ggplot(df, aes(PPCA1, PPCA2)) +
        geom_point(alpha = 0.5) +
        theme_bw()

## ----plot-tre--------------------------------------------------------------
ids <- donor_id(cell_vcf,
                system.file("extdata", "donors.donorid.vcf.gz", 
                                 package = "cardelino"),
                 check_doublet = TRUE)
table(ids$assigned$donor_id)

## ----head-assigned---------------------------------------------------------
head(ids$assigned)

## ----plot-doublet, fig.height=4, fig.width=6-------------------------------
ggplot(ids$assigned, aes(n_vars, prob_doublet, colour = donor_id)) +
    geom_point(alpha = 0.5) +
    theme_bw()

## ----plot-postprob, fig.height=4, fig.width=6------------------------------
ggplot(ids$assigned, aes(n_vars, prob_max, colour = donor_id)) +
    geom_point(alpha = 0.5) +
    scale_x_log10() +
    theme_bw()

## ----plot-ppca-donor, fig.height=4, fig.width=6----------------------------
df$donor_id <- ids$assigned$donor_id
ggplot(df, aes(PPCA1, PPCA2, colour = donor_id)) +
        geom_point(alpha = 0.5) +
        theme_bw()

## ----extreme-first---------------------------------------------------------
donor_vcf <- read_vcf(system.file("extdata", "donors.donorid.vcf.gz", 
                                 package = "cardelino"))
ids_sing <- donor_id(cell_vcf, donor_vcf, check_doublet = FALSE)

## ----extreme-second--------------------------------------------------------
ids_doub <- donor_id(cell_vcf, donor_vcf, check_doublet = TRUE,
                        donors = unique(ids_sing$assigned$donor_id))
table(ids_doub$assigned$donor_id)

## --------------------------------------------------------------------------
sessionInfo()

