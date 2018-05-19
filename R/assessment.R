# Assessment of cardelino's performance on simulated data

#' Assess performance with simulated data
#' 
#' Return overall accuracy, and fraction of assignable cells accuracy with 
#' assignable cells given a threshold of prob_gap
#' 
#' @param prob_mat A matrix (M cells x K clones), estiamted probability of cell 
#' to clone
#' @param simu_mat A matrix (M cells x K clones), tru cell-clone assignment
#' @param threshold A float value, the threshold for assignable cells
#' @param mode A string, the mothod for defining scores for filtering cells: 
#' best and delta. best: highest probability of a cell to K clones, delta: the 
#' difference between the best and second.
#' @export
#' 
assign_score <- function(prob_mat, simu_mat, threshold=0.2, mode="delta"){
  assign_0 <- cardelino::get_prob_label(simu_mat)
  assign_1 <- cardelino::get_prob_label(prob_mat)
  prob_val <- cardelino::get_prob_value(prob_mat, mode=mode)
  idx <- prob_val >= threshold
  
  rt_list <- list("ass"=mean(idx), 
                  "acc"=mean(assign_0 == assign_1),
                  "acc_ass"=mean((assign_0 == assign_1)[idx]))
  rt_list
}

#' Calculate assignability and accuracy along prob_gap
#' 
#' @param prob_mat A matrix (M cells x K clones), estiamted probability of cell 
#' to clone
#' @param simu_mat A matrix (M cells x K clones), true cell-clone assignment
#' @param mode A string, the mothod for defining scores for filtering cells: 
#' best, second and delta. Best: highest probability of a cell to K clones, 
#' similarly for second. delta is the difference between the best and second.
#' 
#' @export
#' 
assign_curve <- function(prob_mat, simu_mat, mode="delta"){
  assign_0 <- cardelino::get_prob_label(simu_mat)
  assign_1 <- cardelino::get_prob_label(prob_mat)
  prob_val <- cardelino::get_prob_value(prob_mat, mode=mode)
  
  thresholds <- sort(unique(prob_val))
  ACC <- rep(0, length(thresholds))
  ASS <- rep(0, length(thresholds))
  for(i in seq_len(length(thresholds))){
    idx <- prob_val >= thresholds[i] 
    #if(thresholds[i] == 0){ idx <- prob_val > thresholds[i] }
    ASS[i] <- mean(idx)
    ACC[i] <- mean((assign_0 == assign_1)[idx])
  }
  thresholds <- c(thresholds, 1.0)
  ACC <- c(ACC, 1.0)
  ASS <- c(ASS, 0.0)
  
  AUC <- AUC_acc_ass <- 0.0
  for(i in seq_len(length(thresholds)-1)){
    AUC <- AUC + 0.5 * (thresholds[i] - thresholds[i+1]) * (ACC[i] + ACC[i+1])
    AUC_acc_ass <- AUC_acc_ass + 0.5 * (ASS[i] - ASS[i+1]) * (ACC[i] + ACC[i+1])
  }
  AUC <- AUC / (thresholds[1] - thresholds[length(thresholds)])
  AUC_acc_ass <- AUC_acc_ass / (ASS[1] - ASS[length(thresholds)])
  
  
  rt_list <- list("ACC"=ACC, "ASS"=ASS, "AUC"=AUC, "AUC_acc_ass"=AUC_acc_ass,
                  "thresholds"=thresholds)
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
  thresholds <- seq(0, 0.999, 0.001)
  ACC <- rep(0, length(thresholds))
  ASS <- rep(0, length(thresholds))
  FPR <- rep(0, length(thresholds))
  TPR <- rep(0, length(thresholds))
  for(i in seq_len(length(thresholds))){
    idx <- prob_mat >= thresholds[i]
    ASS[i] <- mean(idx) # not very meaningful
    ACC[i] <- mean(simu_mat[idx])
    FPR[i] <- sum(simu_mat[idx] == 0) / sum(simu_mat == 0)
    TPR[i] <- sum(simu_mat[idx] == 1) / sum(simu_mat == 1)
  }
  
  AUC <- 0.0
  for(i in seq_len(length(thresholds)-1)){
    AUC <- AUC + 0.5 * (FPR[i] - FPR[i+1]) * (TPR[i] + TPR[i+1])
  }
  
  rt_list <- list("FPR"=FPR, "TPR"=TPR, "AUC"=AUC,
                  "thresholds"=thresholds, "ACC"=ACC, "ASS"=ASS)
  rt_list
}

