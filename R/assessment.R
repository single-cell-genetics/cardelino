# A few basic functions for assessment

#' Convert a matrix to data frame
#' @param X A matrix of values
#' @export
#' @example 
#' mtx_to_df(matrix(seq(12), nrow=3))
mtx_to_df <- function(X) {
  if (is.null(row.names(X))) {
    row.names(X) <- seq(nrow(X))
  }
  if (is.null(colnames(X))) {
    colnames(X) <- seq(ncol(X))
  }
  data.frame(Var1 = rep(row.names(X), ncol(X)),
             Var2 = rep(colnames(X), each = nrow(X)),
             value = c(X))
}

#' Maximum value for each row in a matrix
#'
#' @param X A matrix of floats. 
#' @param mode A string, the mothod for defining scores for filtering cells:
#' best, second and delta. \code{best}: highest value for each row, similarly 
#' for the \code{second}. \code{delta} is the difference between the best and 
#' the second.
#' @return a vector of the collapsed value for each row, depending
#' on the mode used.
#'
#' @export
#'
rowMax <- get_prob_value <- function(X, mode = "best") {
    max_val <- rep(0, nrow(X))
    for (i in seq_len(nrow(X))) {
        sorted_val <- sort(X[i,], decreasing = TRUE)
        if (mode == "delta") {
            max_val[i] <- sorted_val[1] - sorted_val[2]
        } else if (mode == "second") {
            max_val[i] <- sorted_val[2]
        } else {#default mode: best
            max_val[i] <- sorted_val[1]
        }
    }
    max_val
}


#' Column index of the maximum value for each row in a matrix
#'
#' @param X A matrix of floats. 
#' @return a vector of the index of column for each row. Note, when multiple
#' columns have the same value, only the earliest column will be
#' returned.
#'
#' @export
#'
rowArgmax <- get_prob_label <- function(X){
    max_idx <- rep(0, nrow(X))
    for (j in seq_len(nrow(X))) {
        max_idx[j] <- which.max(X[j, ])
    }
    max_idx
}

#' Column match between two matrices by the minimum mean abosolute difference
#' @param A The first matrix which will be matched
#' @param B The second matrix, the return index will be used on
#' @param force bool(1), If TRUE, force traversing all permutations of B to 
#' find the optimised match to A with computing cost of O(n!). Otherwise, use
#' greedy search with computing cost of O(n^2).
#' @return \code{idx}, the column index of B to be matched to A
#' @export
#' 
colMatch <- function(A, B, force = FALSE) {
    if (nrow(A) != nrow(B)) {
        stop("Error: A and B have different rows.")
    }
    if (ncol(A) > ncol(B)) {
        stop("Error: A must have equal or smaller columns than B.")
    }
    if (force == TRUE) {
        idx_list <- combinat::permn(ncol(B))
        MAE_list <- rep(NA, length(idx_list))
        for (ii in seq_len(length(idx_list))) {
            MAE_list[ii] <- mean(abs(A - B[, idx_list[[ii]][1:ncol(A)]]))
        }
        idx_out <- idx_list[[which.min(MAE_list)]][1:ncol(A)]
    } else {
        MAE_mat = matrix(0, ncol(A), ncol(B))
        for (i in seq_len(ncol(A))){
            for (j in seq_len(ncol(B))){
                MAE_mat[i, j] = mean(abs(A[, i] - B[, j]))
            }
        }
        MAE_tmp = MAE_mat
        idx_out = rep(-1, ncol(A))
        while (-1 %in% idx_out) {    
            idx_tmp = which(MAE_tmp == min(MAE_tmp), arr.ind = TRUE)
            idx_row = idx_tmp[1]
            idx_col = idx_tmp[2]
            
            idx_out[idx_row] = idx_col
            MAE_tmp[idx_row, ] = max(MAE_mat) + 1
            MAE_tmp[, idx_col] = max(MAE_mat) + 1
        }
    }
    idx_out
}


#' Precision-recall curve for binary label prediction
#' @param scores Prediction score for each sample
#' @param labels True labels for each sample, e.g., from simulation
#' @param cutoff A list of cutoff; if NULL use all unique scores
#' @param cut_direction A string to compare with cutoff: >=, >, <=, <
#' @param add_cut1 Logical value; if True, manually add a cutoff of 1
#' @param empty_precision Float value for default precision if no any recall
#' @export
#' 
binaryPRC <- function(scores, labels, cutoff=NULL, cut_direction=">=", 
                      add_cut1=FALSE, empty_precision=1) {
    if (is.null(cutoff)) {
        cutoff <- sort(unique(scores))
    }
    n_positive <- sum(labels == 1 | labels)
    
    Precision <- rep(0, length(cutoff))
    Recall <- rep(0, length(cutoff))
    for (i in seq_len(length(cutoff))) {
        if (cut_direction == "<=") {
            idx <- scores <= cutoff[i]
        } else if (cut_direction == "<") {
            idx <- scores < cutoff[i]
        } else if (cut_direction == ">") {
            idx <- scores > cutoff[i]
        } else {
            idx <- scores >= cutoff[i]
        }
        
        is_idx_true <- (labels[idx] == 1 | labels[idx])
        Recall[i] <- sum(is_idx_true) / n_positive
        Precision[i] <- mean(is_idx_true)
    }
    if (!is.null(empty_precision)) { 
            Precision[Recall == 0] <- empty_precision
    }
    if (add_cut1) {
        cutoff <- c(cutoff, 1.0)
        Recall <- c(Recall, 0.0)
        Precision <- c(Precision, 1.0)
    }
    
    AUC <- 0.0
    for (i in seq_len(length(cutoff) - 1)) {
        AUC <- ((Recall[i] - Recall[i + 1]) * 
                (Precision[i] + Precision[i + 1]) * 0.5 + AUC)
    }
    AUC <- AUC / (Recall[1] - Recall[length(Recall)])
    
    df <- data.frame("cutoff" = cutoff, "Recall" = Recall, 
                     "Precision" = Precision)
    list("df" = df, "AUC" = AUC)
}


binaryROC <- function(scores, labels, cutoff=NULL, cut_direction=">=",
                      add_cut1=TRUE, cutoff_point=0.9) {
  if (is.null(cutoff)) {
    cutoff <- sort(unique(scores))
  }
  cutoff <- sort(unique(c(cutoff, cutoff_point, 0, 1)))
  
  TPR <- rep(0, length(cutoff))
  FPR <- rep(0, length(cutoff))
  Precision <- rep(0, length(cutoff))
  for (i in seq_len(length(cutoff))) {
    if (cut_direction == "<=") {
      idx <- scores <= cutoff[i]
    } else if (cut_direction == "<") {
      idx <- scores < cutoff[i]
    } else if (cut_direction == ">") {
      idx <- scores > cutoff[i]
    } else {
      idx <- scores >= cutoff[i]
    }
    
    FPR[i] <- sum(labels[idx] == 0) / sum(labels == 0) ## FPR
    TPR[i] <- sum(labels[idx] == 1) / sum(labels == 1) ## TPR
    Precision[i] <- mean(labels[idx] == 1)
  }    
  Precision[(TPR == 0) & is.na(Precision)] <- 1.0
  
  if (add_cut1) {
    cutoff <- c(cutoff, 1.0)
    FPR <- c(FPR, 0.0)
    TPR <- c(TPR, 0.0)
    Precision <- c(Precision, 1.0)
  }
  
  AUC <- AUPRC <- 0.0
  for (i in seq_len(length(cutoff) - 1)) {
    AUC <- ((FPR[i] - FPR[i + 1]) * 
              (TPR[i] + TPR[i + 1]) * 0.5 + AUC)
    AUPRC <- ((TPR[i] - TPR[i + 1]) * 
                (Precision[i] + Precision[i + 1]) * 0.5 + AUPRC)
  }
  AUC <- AUC / (FPR[1] - FPR[length(FPR)])
  AUPRC <- AUPRC / (TPR[1] - TPR[length(TPR)])
  
  df <- data.frame("cutoff" = cutoff, "TPR" = TPR, 
                   "FPR" = FPR, Precision = Precision)
  list("df" = df, "AUC" = AUC, "AUPRC" = AUPRC)
}


#' Precision-recall curve for multi-class prediction
#' @param prob_mat Probability matrix for each cell to each component
#' @param simu_mat The true identity of assignment from simulation
#' @param marginal_mode A string for the mode to marginalize the column: best, 
#' second, or delta
#' @param cutoff A list of cutoff; if NULL use all unique scores
#' @param add_cut1 Logical value; if True, manually add a cutoff of 1
#' @param multiLabel.rm Logical value; if True, remove the samples with 
#' multiple labels
#' @export
#' 
multiPRC <- function(prob_mat, simu_mat, marginal_mode="best", 
                     cutoff=NULL, multiLabel.rm = TRUE, add_cut1=FALSE){
    if (nrow(prob_mat) != nrow(simu_mat) || 
        ncol(prob_mat) != ncol(simu_mat)) {
        stop("The shapes are not matched: prob_mat, simu_mat.")
    }
    
    idx_multiLabel <- rowSums(simu_mat > 0) > 1
    if (sum(idx_multiLabel)) {
        print(paste(sum(idx_multiLabel), "samples have multiple labels."))
    }
    
    assign_0 <- rowArgmax(simu_mat)
    assign_1 <- rowArgmax(prob_mat)
    assign_0[idx_multiLabel] <- -1
    
    if (marginal_mode == "sum") {
        prob_val <- rowSums(prob_mat)
    } else {
        prob_val <- rowMax(prob_mat, mode = marginal_mode)
    }
    
    if (multiLabel.rm) {
        prob_val <- prob_val[assign_0 > 0]
        assign_1 <- assign_1[assign_0 > 0]
        assign_0 <- assign_0[assign_0 > 0]
    }
    
    if (is.null(cutoff)) {
        cutoff <- sort(unique(prob_val))
    }
    
    Precision <- rep(0, length(cutoff))
    Recall <- rep(0, length(cutoff))
    for (i in seq_len(length(cutoff))) {
        idx <- prob_val >= cutoff[i]
        Recall[i] <- mean(idx)
        Precision[i] <- mean((assign_0 == assign_1)[idx])
    }
    if (add_cut1) {
        cutoff <- c(cutoff, 1.0)
        Recall <- c(Recall, 0.0)
        Precision <- c(Precision, 1.0)
    }
    
    AUC <- 0.0
    for (i in seq_len(length(cutoff) - 1)) {
        AUC <- ((Recall[i] - Recall[i + 1]) * 
                (Precision[i] + Precision[i + 1]) * 0.5 + AUC)
    }
    AUC <- AUC / (Recall[1] - Recall[length(Recall)])
    
    df <- data.frame("cutoff" = cutoff, "Recall" = Recall, 
                     "Precision" = Precision)
    list("df" = df, "AUC" = AUC)
}

#' Scoring the simulation in assignment of singlets and doublets
#' @param prob Probability matrix for each cell to each component
#' @param I_sim The true identity of assignment from simulation
#' @param cutoff A list of cutoff from 0 to 1
#' @export
#' 
assign_scores <- function(prob, I_sim, cutoff=seq(0, 1, 0.001)) {
    col_idx <- colMatch(I_sim, prob)
    print(t(col_idx))

    ass_sg <- multiPRC(prob[, col_idx], I_sim, multiLabel.rm = TRUE, 
                       cutoff = cutoff)
    ass_sg$AUC
    
    ass_db <- binaryPRC(rowSums(prob), rowSums(I_sim > 0) > 1,
                        cut_direction = "<=", cutoff = cutoff)
    ass_db$AUC
    
    list("df_sg" = ass_sg$df, "df_db" = ass_db$df, 
         "AUC_sg" = ass_sg$AUC, "AUC_db" = ass_db$AUC)
}

