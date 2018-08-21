# A few basic utility functions

#' Get the collapsed probability value
#'
#' @param prob_assign A matrix of floats. Clonal assignment probability of M
#' cells to K clones.
#' @param mode A string, the mothod for defining scores for filtering cells:
#' best, second and delta. \code{best}: highest probability of a cell to K
#' clones, similarly for the \code{second}. \code{delta} is the difference
#' between the best and the second.
#' @return a vector of the collapsed probability value for each cell, depending
#' on the mode used.
#'
#' @export
#' @examples
#' data(example_donor)
#' assignments <- clone_id(A_clone, D_clone, Config = tree$Z, inference = "EM")
#' prob_val <- get_prob_value(assignments$prob, mode = "best")
#'
get_prob_value <- function(prob_assign, mode = "best") {
    prob_val <- rep(0, nrow(prob_assign))
    for (i in seq_len(length(prob_val))) {
        prob_sorted <- sort(prob_assign[i,], decreasing = TRUE)
        if (mode == "delta") {
            prob_val[i] <- prob_sorted[1] - prob_sorted[2]
        } else if (mode == "second") {
            prob_val[i] <- prob_sorted[2]
        } else {#default mode: best
            prob_val[i] <- prob_sorted[1]
        }
    }
    prob_val
}


#' Get the clone label from the assignmnet probabilities.
#'
#' @param prob_assign A matrix of floats. Clonal assignment probability of N
#' cells to K clones.
#' @return a vector of the index of clone for each cell. Note, when multiple
#' clones have the same assignment probability, only the earliest clone will be
#' return. Normally, these cells should be filtered for analysis.
#'
#' @export
#' @examples
#' data(example_donor)
#' assignments <- clone_id(A_clone, D_clone, Config = tree$Z, inference = "EM")
#' clone_idx <- get_prob_label(assignments$prob)
#'
get_prob_label <- function(prob_assign){
    assign_clone <- rep(0, nrow(prob_assign))
    for (j in seq_len(length(assign_clone))) {
        assign_clone[j] <- which.max(prob_assign[j, ])
    }
    assign_clone
}
