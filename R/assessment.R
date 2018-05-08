# Assessment of cardelino's performance on simulated data

#' Assess performance with simulated data
#' 
#' Return overall accuracy, and fraction of assignable cells accuracy with 
#' assignable cells given a threshold of prob_gap
#' 
#' @param prob_mat A matrix (M cells x K clones), estiamted probability of cell 
#' to clone
#' @param simu_mat A matrix (M cells x K clones), tru cell-clone assignment
#' @param prob_gap_min A float value, the threshold for assignable cells
#' @export
#' 
assign_score <- function(prob_mat, simu_mat, prob_gap_min=0.2){
  assign_0 <- cardelino::get_prob_label(simu_mat)
  assign_1 <- cardelino::get_prob_label(prob_mat)
  prob_gap <- cardelino::get_prob_gap(prob_mat)
  idx <- prob_gap >= prob_gap_min
  
  #return: assignable cells, accuracy with assignable cells, overall accuracy
  c(mean(idx), mean((assign_0 == assign_1)[idx]), mean(assign_0 == assign_1))
}

#' Calculate assignability and accuracy along prob_gap
#' 
#' @param prob_mat A matrix (M cells x K clones), estiamted probability of cell 
#' to clone
#' @param simu_mat A matrix (M cells x K clones), tru cell-clone assignment
#' @export
#' 
assign_curve <- function(prob_mat, simu_mat){
  assign_0 <- cardelino::get_prob_label(simu_mat)
  assign_1 <- cardelino::get_prob_label(prob_mat)
  prob_gap <- cardelino::get_prob_gap(prob_mat)
  
  GAP <- sort(unique(prob_gap))
  ACC <- rep(0, length(GAP))
  ASS <- rep(0, length(GAP))
  for(i in seq_len(length(GAP))){
    idx <- prob_gap >= GAP[i] 
    #if(GAP[i] == 0){ idx <- prob_gap > GAP[i] }
    ASS[i] <- mean(idx)
    ACC[i] <- mean((assign_0 == assign_1)[idx])
  }
  GAP <- c(GAP, 1.0)
  ACC <- c(ACC, 1.0)
  ASS <- c(ASS, 0.0)
  
  AUC <- 0.0
  for(i in seq_len(length(GAP)-1)){
    AUC <- AUC + 0.5 * (ASS[i] - ASS[i+1]) * (ACC[i] + ACC[i+1])
  }
  AUC <- AUC / (ASS[1] - ASS[length(GAP)])
  
  rt_list <- list("GAP"=GAP, "ACC"=ACC, "ASS"=ASS, "AUC"=AUC)
  rt_list
}


#' Calculate macro ROC
#' 
#' @param prob_mat A matrix (M cells x K clones), estiamted probability of cell 
#' to clone
#' @param simu_mat A matrix (M cells x K clones), tru cell-clone assignment
#' @export
#' 
assign_macro_ROC <- function(prob_mat, simu_mat){
  thresholds <- seq(0, 1, 0.001)
  FPR <- rep(0, length(thresholds))
  TPR <- rep(0, length(thresholds))
  for(i in seq_len(length(thresholds))){
    pred_mat <- prob_mat >= thresholds[i]
    FPR[i] <- sum((simu_mat == 0) & (pred_mat == 1)) / sum(simu_mat == 0)
    TPR[i] <- sum((simu_mat == 1) & (pred_mat == 1)) / sum(simu_mat == 1)
  }
  
  AUC <- 0.0
  for(i in seq_len(length(thresholds)-1)){
    AUC <- AUC + 0.5 * (FPR[i] - FPR[i+1]) * (TPR[i] + TPR[i+1])
  }
  
  rt_list <- list("thresholds"=thresholds, "FPR"=FPR, "TPR"=TPR, "AUC"=AUC)
  rt_list
}

