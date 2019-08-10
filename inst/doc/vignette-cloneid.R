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
data(example_donor)

## ----plot-tre--------------------------------------------------------------
plot_tree(tree, orient = "v")

## ----cell-assign-----------------------------------------------------------
assignments <- clone_id(A_clone, D_clone, Config = tree$Z)
names(assignments)

## ----prob-heatmap----------------------------------------------------------
prob_heatmap(assignments$prob)

## ----assign-cell-clone-easy------------------------------------------------
df <- assign_cells_to_clones(assignments$prob)
head(df)
table(df$clone)

## ----read-vcf-data---------------------------------------------------------
vcf <- read_vcf(system.file("extdata", "cell_example.mpileup.vcf.gz", 
                   package = "cardelino"))
input_data <- get_snp_matrices(vcf)

## ----read-canopy-data------------------------------------------------------
canopy <- readRDS(system.file("extdata", "canopy_results.example.rds", 
                   package = "cardelino"))
Config <- canopy$tree$Z

## ----correct-variant-ids---------------------------------------------------
rownames(Config) <- gsub("chr", "", gsub(":", "_", gsub("_.*", "", 
                                                        rownames(Config))))

## ----run-cell-assign-------------------------------------------------------
assignments <- clone_id(input_data$A, input_data$D, Config)

## ----assign-cell-clone-vcf-------------------------------------------------
df <- assign_cells_to_clones(assignments$prob)
table(df$clone)

## --------------------------------------------------------------------------
sessionInfo()

