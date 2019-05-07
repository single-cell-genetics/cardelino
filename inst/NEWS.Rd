\name{NEWS}
\title{News for Package \pkg{cardelino}}

\section{Changes in version 0.3.9}{
  \itemize{
    \item Add a runnable R file: run_vireo.R in /bin to run from command line 
    \item Change load_GT_vcf by supporting GP for genotype probability
    \item Change load_cellSNP_vcf by supporting GL for genotype likelihood
    \item Change heat_matrix to for scale_x_discrete by default
  }
}

\section{Changes in version 0.3.8}{
  \itemize{
    \item Change the covergency diagnosis in cell_assign_Gibbs to all cells
    \item Change the donor_read_simulator to betabinomial for each variant in 
          each cell
  }
}
\section{Changes in version 0.3.7}{
  \itemize{
    \item Fix a bug in colMatch in `force` mode
    \item Specify the fixed relax rate in cell_assign_Gibbs, which can be used 
          to turn the input Config as a uniform prior when setting 
          relax_rate_fixed=0.5
  }
}
\section{Changes in version 0.3.4}{
  \itemize{
    \item Add averaged `Config_prob` and `relax_rate` in cell_assign_Gibbs 
          outputs
    \item Add `force` option in colMatch function to force one-to-one match
  }
}

\section{Changes in version 0.3.3}{
  \itemize{
    \item Fix the bug in inferring the relax_rate in cell_assign_Gibbs
  }
}

\section{Changes in version 0.3.2}{
  \itemize{
    \item Supporting learning the relax rate on clone configuration 
          automatically
  }
}

\section{Changes in version 0.3.1}{
  \itemize{
    \item fix the clone id missing for prob_mat.
    \item change the default the sampling iteration to 3000.
  }
}

\section{Changes in version 0.3.0}{
  \itemize{
    \item cell_assign_Gibbs supports updating clone Configuration now. Set 
          relax_Config between 0 and 1.
    \item change the default parameters for beta prior to better represent 
          allelic dropout and imblance in scRNA-seq data. It involves fucntions:
          sim_read_count, donor_read_simulator, and cell_assign_Gibbs
    \item minor change of pub.theme: title will be plain rather than bold.
  }
}

\section{Changes in version 0.2.7}{
  \itemize{
    \item vireo supports match SNP from donor_data to cel_data; change default
          number of processors to n_proc=1
    \item change load_cellSNP_vcf default paramters to support more general case
    \item add more dependency to pass tests
    \item correct test-donor_id.R
    \item remove vignette-donorid.Rmd vignette and correct vignette-vireo.Rmd
  }
}

\section{Changes in version 0.2.5}{
  \itemize{
    \item change donor_id to vireo
    \item fix minor bug for n_vars in vireo (i.e., donor_id)
    \item add vignette for demultiplexing without genotype
    \item change assessment for doublet detection indicator to prob_doublet
  }
}

\section{Changes in version 0.1.0}{
  \itemize{
    \item add examples and remove unnecessary functions to pass biocCheck
  }
}
