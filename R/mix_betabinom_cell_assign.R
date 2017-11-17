# beta-binomial mixture model to assign cells to clones

#' EM algorithm for assigning cells to cloens and estimating beta-binomial mixture model.
#' The beta-binomial will capture the varibility across cells, which covers false 
#' positive and false negative.
#' As there is no closed-form for MLE estimate on beta-binomial in the M-step, we use a
#' hard assignment of cells in the E-step and optimize the parameters numerically in the
#' M-step by using VGAM package.
#' 
#' @param A A matrix of integers. Number of alteration reads in variant i in cell j
#' @param D A matrix of integers. Number of total reads (depth) in variant i in cell j
#' @param C A matrix of binary values. The clone-variant configuration, whcih encodes
#' the phylogenetic tree structure. This is the output Z from Canopy
#' @param Psi A vector of float. The fractional size of clone, output P from Canopy
#' @param max_iter A integer. The maximum number of iterations in EM algorithm. 
#' The real iteration may finish earlier.
#' 
#' @return a list containing \code{u}, a float denoting the estimated false positive 
#' rate, \code{v}, a float denoting the estimated false negative rate, \code{prob}, the 
#' matrix of fitted probabilities of each cell belonging to each clone, and \code{logLik},
#' the log likelihood of the parameter based on the final cell assignment.
#' 
#' @import stats, VGAM, bbmle
#' 
#' @export
#' 
#' @examples
#' 
mix_betaBinomial_tree <- function(A, D, C, Psi, max_iter=100){
  if (dim(A)[1] != dim(D)[1] || dim(A)[2] != dim(D)[2] || dim(A)[1] != dim(C)[1] || 
      dim(C)[2] != length(Psi)) {
    stop(paste("A and D must have the same size;\n ",
               "A and C must have the same number of variants;\n",
               "C and Psi must have the same number of clones",sep = ""))
  }
  
  ## preprocessing
  N <- dim(A)[1] # number of variants 
  M <- dim(A)[2] # number of cells
  K <- dim(C)[2] # number of clones
  
  ## random initialization for EM
  mu <- c(0.05, 0.45)
  rho <- c(0.1, 0.7)
  mu_new <- stats::runif(2, 0.001, 0.5)
  rho_new <- stats::runif(2, 0.1,  0.7)
    
  lik_mat <- matrix(1.0, M, K)
  prob_mat <- lik_mat / rowSums(lik_mat)
  
  ## EM iterations
  for(t in seq_len(max_iter)){
    #TODO: check convergence by the fold change on logLik
    if (mu[1] == mu_new[1] && mu[2] == mu_new[2] && rho[1] == rho_new[1] &&
        rho[2] == rho_new[2]){
      print(paste("Total iterations: ", t))
      break
    } else{
      mu <- mu_new
      rho  <- rho_new
    }
    
    #E-step
    for(j in seq_len(M)){
      for(k in seq_len(K)){
        lik_mat[j,k] <- Psi[k]
        for(i in seq_len(N)){
          if(is.na(A[i,j])==FALSE && C[i,k]==0){
            lik_mat[j,k] <- lik_mat[j,k] * VGAM::pbetabinom(A[i,j], size=D[i,j], prob=mu[1], rho[1])
          }else if (is.na(A[i,j])==FALSE && C[i,k]==1){
            lik_mat[j,k] <- lik_mat[j,k] * VGAM::pbetabinom(A[i,j], size=D[i,j], prob=mu[2], rho[2])
          }
        }
      }
    }
    
    Assign <- argmax(lik_mat, rows = TRUE)
    Cell_Clone <- C[,Assign]
    
    idx1 <- which((Cell_Clone==0) * (is.na(A)==FALSE) == 1)
    idx2 <- which((Cell_Clone==1) * (is.na(A)==FALSE) == 1)
    
    #M-step with numerical optimization
    x <- A[idx1]
    y <- D[idx1]-A[idx1]
    model1 <- bbmle::mle2(y~dbetabinom.ab(size=x+y,shape1,shape2),
                   data=data.frame(x,y), start=list(shape1=2,shape2=30))
    
    x <- A[idx2]
    y <- D[idx2]-A[idx2]
    model2 <- bbmle::mle2(y~dbetabinom.ab(size=x+y,shape1,shape2),
                   data=data.frame(x,y),
                   start=list(shape1=2,shape2=30))
    
    mu_new[1] <- model1@coef[["shape1"]] / (model1@coef[["shape1"]] + model1@coef[["shape2"]])
    mu_new[2] <- model2@coef[["shape1"]] / (model2@coef[["shape1"]] + model2@coef[["shape2"]])
    
    rho_new[1] <- 1 / (1 + model1@coef[["shape1"]] + model1@coef[["shape2"]])
    rho_new[2] <- 1 / (1 + model2@coef[["shape1"]] + model2@coef[["shape2"]])
  }
  
  prob_mat <- lik_mat / rowSums(lik_mat)
  logLik <- sum(log(rowSums(lik_mat)))
  
  ## return values
  return_list <- list("alpha"=alpha, "beta"=beta, "prob"=prob_mat, "logLik"=logLik)
  return_list
}

