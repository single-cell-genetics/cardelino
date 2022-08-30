# News for Package `cardelino`

## Changes in v0.99.1
* Changes to prepare package for submission to Bioconductor. (No functional
  changes)

## Changes in v0.6.4 (21/08/2019)
* Fix bug in devarianceIC() = D_post + 2 * p_D; and add both the orignal p_D, 
  i.g., D_mean - D_post and the Gelman's alternative p_D = 2 * var(D). By 
  default, DIC uses Gelman's method.
* Suggest nerrowing theta1 for cardelino-free, e.g., prior1=c(45, 55)
* Add get_logLik() function to get log likelihood for clone_id_Gibbs(), namely
  P(A, D | C, I, theta0, theta1)

## Changes in v0.6.3 (19/08/2019)
* Add devarianceIC() for Devariance information criterion (DIC) for 
  clone_id_Gibbs()
* Change the default of relax_Config=TRUE for clone_id() and clone_id_Gibbs()

## Changes in version 0.6.2
* Fix the bug in `Geweke_Z` function for convegence diagnostic
* Change the default number of iterations in `clone_id` to 5000
* Add more tests

## Changes in version 0.6.1
* Change the default colMatch to greedy search
* Rename the cell_assign_Gibbs to clone_id_Gibbs
* Rename the cell_assign_EM to clone_id_EM
* Move get_tree from clone_id.R to tree_utils.R
* Remove "Bernoulli" model in both clone_id_EM and clone_id_Gibbs. The 
  user can transfer the binomial to Bernoulli by setting a threhold 
  beforehand
* Remove the A_germ and D_germ parameters, leave it to future 
  development to jointly modelling the germline variants
* Remove fucntion load_vcf_h5 and dependency hdf5r


## Changes in version 0.6.0
* Remove donor id relevant informatio and backup donor_id into `v0.4.2` 
  and a new branch `with_vireo`
* Remove donor_id.R, bin/run_vireo.R
* Remove tests/testthat/test-donorid.R, vignette/vignette-vireo.Rmd
* Remove unnecessary logo photos in inst


## Changes in version 0.4.2
* Add binaryROC function of ROC curve
* Add message for switch vireo to Python implementation


## Changes in version 0.4.1
* Support direct loading sparse matrix from cellSNP file in hdf5 format


## Changes in version 0.4.0
* Support multiple chains when running clone_id (one chain by default)
* Support relabel during Gibbs sampling when running clone_id (not in 
  use by default)
* Change sim_read_count: the default of wise0 is to "element", i.e.,
  sequencing error rate is both variant and cell specific. Also, a new
  parameter sample_cell is added for the option to turn off sampling 
  cells from seed D
* Change heat_matrix for supporting rownames and colnames matching the 
  original order in the matrix instead of alphabetical order


## Changes in version 0.3.9
* Add a runnable R file: run_vireo.R in /bin to run from command line 
* Change load_GT_vcf by supporting GP for genotype probability
* Change load_cellSNP_vcf by supporting GL for genotype likelihood
* Change heat_matrix to for scale_x_discrete by default


## Changes in version 0.3.8
* Change the covergency diagnosis in cell_assign_Gibbs to all cells
* Change the donor_read_simulator to betabinomial for each variant in 
  each cell


## Changes in version 0.3.7
* Fix a bug in colMatch in `force` mode
* Specify the fixed relax rate in cell_assign_Gibbs, which can be used 
  to turn the input Config as a uniform prior when setting 
  relax_rate_fixed=0.5


## Changes in version 0.3.4
* Add averaged `Config_prob` and `relax_rate` in cell_assign_Gibbs 
  outputs
* Add `force` option in colMatch function to force one-to-one match


## Changes in version 0.3.3
* Fix the bug in inferring the relax_rate in cell_assign_Gibbs


## Changes in version 0.3.2
* Supporting learning the relax rate on clone configuration 
  automatically


## Changes in version 0.3.1
* fix the clone id missing for prob_mat.
* change the default the sampling iteration to 3000.


## Changes in version 0.3.0
* cell_assign_Gibbs supports updating clone Configuration now. Set 
  relax_Config between 0 and 1.
* change the default parameters for beta prior to better represent 
  allelic dropout and imblance in scRNA-seq data. It involves fucntions:
  sim_read_count, donor_read_simulator, and cell_assign_Gibbs
* minor change of pub.theme: title will be plain rather than bold.


## Changes in version 0.2.7
* vireo supports match SNP from donor_data to cel_data; change default
  number of processors to n_proc=1
* change load_cellSNP_vcf default paramters to support more general case
* add more dependency to pass tests
* correct test-donor_id.R
* remove vignette-donorid.Rmd vignette and correct vignette-vireo.Rmd


## Changes in version 0.2.5
* change donor_id to vireo
* fix minor bug for n_vars in vireo (i.e., donor_id)
* add vignette for demultiplexing without genotype
* change assessment for doublet detection indicator to prob_doublet


## Changes in version 0.1.0
* add examples and remove unnecessary functions to pass biocCheck
