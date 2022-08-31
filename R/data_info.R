### File: config_all

#' @name Config_all
#' @title A list of tree configuration
#' @description This list of tree configuration between 3 clones to 10 clones,
#' each element is a list with all possible tree matrix
#' @return NULL, but makes available a list
#' @docType data
#' @usage config_all
#' @format a list of list of matrix
#' @source PASTRI Python package
#' @author Yuanhua Huang, Davis McCarthy, 2018-06-25
NULL


### File: simulation_input

#' @name D_input
#' @title A matrix of sequencing depths
#' @description This matrix contains sequencing depths for 439 somatic variants
#' across 151 cells, from one particular scRNA-seq sample, can be used to
#' generate sequencing depths
#' @return NULL, but makes available a matrix
#' @docType data
#' @usage simulation_input
#' @format a matrix of float
#' @source A fibroblast sample from HipSci project
#' @author Yuanhua Huang, Davis McCarthy, 2018-06-25
NULL

#' @name tree_3clone
#' @title A tree object
#' @description This tree object with 3 clones contains clonal tree information,
#' inferred from bulk exome-seq data
#' @return NULL, but makes available a tree object
#' @docType data
#' @usage simulation_input
#' @format a tree object
#' @source A fibroblast sample from HipSci project
#' @author Yuanhua Huang, Davis McCarthy, 2018-06-25
NULL

#' @name tree_4clone
#' @title A tree object
#' @description This tree object with 4 clones contains clonal tree information,
#' inferred from bulk exome-seq data
#' @return NULL, but makes available a tree object
#' @docType data
#' @usage simulation_input
#' @format a tree object
#' @source A fibroblast sample from HipSci project
#' @author Yuanhua Huang, Davis McCarthy, 2018-06-25
NULL

#' @name tree_5clone
#' @title A tree object
#' @description This tree object with 5 clones contains clonal tree information,
#' inferred from bulk exome-seq data
#' @return NULL, but makes available a tree object
#' @docType data
#' @usage simulation_input
#' @format a tree object
#' @source A fibroblast sample from HipSci project
#' @author Yuanhua Huang, Davis McCarthy, 2018-06-25
NULL


### File: example_donor

#' @name tree
#' @title A tree object
#' @description This tree object contains clonal tree information, inferred
#' from bulk exome-seq data
#' @return NULL, but makes available a tree object
#' @docType data
#' @usage example_donor
#' @format a tree object
#' @source A fibroblast sample from HipSci project
#' @author Yuanhua Huang, Davis McCarthy, 2018-06-25
NULL

#' @name D_clone
#' @title A matrix of sequencing depths for clone ID
#' @description This matrix contains sequencing depths for 34 somatic variants
#' across 428 cells, from one example scRNA-seq sample
#' @return NULL, but makes available a matrix
#' @docType data
#' @usage example_donor
#' @format a matrix of float
#' @source A fibroblast sample from HipSci project
#' @author Yuanhua Huang, Davis McCarthy, 2018-06-25
NULL

#' @name A_clone
#' @title A matrix of read numbers of alternative alleles for clone ID
#' @description This matrix contains read numbers of alternative alleles for
#' 34 somatic variants across 428 cells, from one example scRNA-seq sample
#' @return NULL, but makes available a matrix
#' @docType data
#' @usage example_donor
#' @format a matrix of float
#' @source A fibroblast sample from HipSci project
#' @author Yuanhua Huang, Davis McCarthy, 2018-06-25
NULL

#' @name D_germline
#' @title A matrix of sequencing depths
#' @description This matrix contains sequencing depths for 34 germline variants
#' (near the somatic variants) across 428 cells, from one example scRNA-seq
#' sample
#' @return NULL, but makes available a matrix
#' @docType data
#' @usage example_donor
#' @format a matrix of float
#' @source A fibroblast sample from HipSci project
#' @author Yuanhua Huang, Davis McCarthy, 2018-06-25
NULL

#' @name A_germline
#' @title A matrix of read numbers of alternative alleles
#' @description This matrix contains read numbers of alternative alleles for
#' 34 germline variants (near the somatic variants) across 428 cells,
#' from one example scRNA-seq sample
#' @return NULL, but makes available a matrix
#' @docType data
#' @usage example_donor
#' @format a matrix of float
#' @source A fibroblast sample from HipSci project
#' @author Yuanhua Huang, Davis McCarthy, 2018-06-25
NULL
