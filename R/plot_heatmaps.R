# Functions for heatmap plots

#' Plot a heatmap for probability of clone assignment
#' 
#' @param prob_mat A matrix (M x K), the probability of cell j to clone k
#' @param prob_gap_min A float, the threshold for assignable cells, i.e., the 
#' minimum probability gap between the highest and the second clones. It is used
#' for calculating the fraction of clone
#' @param cell_idx A vector the indices of the input cells. If NULL, order by 
#' the probabilty of each clone
#' 
#' @export
prob_heatmap <- function(prob_mat, prob_gap_min=0.2, cell_idx=NULL){
  cell_label <- cardelino::get_prob_label(prob_mat)
  prob_gap <- cardelino::get_prob_gap(prob_mat)
  # add clone id
  for (i in seq_len(ncol(prob_mat))){
    conf_frac <- mean(cell_label[prob_gap>=prob_gap_min] == i)        
    colnames(prob_mat)[i] <- paste0("C", i, ": ", 
                                    round(conf_frac*100,digits=1), "%")
  }
  
  if(is.null(cell_idx)){
    cell_idx <- order(cell_label - diag(prob_mat[, cell_label]))
  }
  nba.m <- reshape2::melt(prob_mat[cell_idx,])
  colnames(nba.m) <- c("Cell", "Clone", "Prob")

  fig_assign <- ggplot(nba.m, aes(Clone, Cell, fill = Prob)) +
    geom_tile(show.legend=T) +
    scale_fill_gradient(low = "white", high = "firebrick4") +
    ylab(paste("Clonal assignment:", nrow(prob_mat), "cells")) +
    cardelino::heatmap.theme() # + cardelino::pub.theme()
  
  fig_assign
}


#' Plot a heatmap for number of mutation sites in each cell
#' 
#' @param Config A matrix (N x K), clonal genotype configuration
#' @param prob_mat A matrix (M x K), the probability of cell j to clone k
#' @param A A matrix (N x M), the present of alternative reads. NA means missing
#' @param mode A string: present or absent
#' 
#' @export
sites_heatmap <- function(Config, A, prob_mat, mode="present"){
  mut_label <- Config %*% (2**seq(ncol(Config),1))
  mut_label <- seq(length(unique(mut_label)),1)[as.factor(mut_label)]
  mut_uniq <- sort(unique(mut_label))

  A_cnt <- A
  if(mode=="absent"){
    A_cnt[is.na(A_cnt)] <- 1
    A_cnt[which(A_cnt>0)] <- 1
    A_cnt <- 1 - A_cnt
  }else{
    A_cnt[is.na(A_cnt)] <- 0
    A_cnt[which(A_cnt>0)] <- 1
  }
  
  mut_mat <- matrix(0, nrow=ncol(A_cnt), ncol=length(mut_uniq))
  for(i in seq_len(length(mut_uniq))){
    idx.tmp <- mut_label == mut_uniq[i]
    mut_mat[,i] <- colSums(A_cnt[idx.tmp, ])
  }
  colnames(mut_mat) <- paste0("Mut", mut_uniq, ": ", table(mut_label))
  print(table(mut_label))
  
  #order by assignment probability
  cell_label <- cardelino::get_prob_label(prob_mat)
  idx <- order(cell_label - diag(prob_mat[, cell_label]))
  mut_mat <- mut_mat[idx,]
  
  nba.m <- reshape2::melt(mut_mat)
  colnames(nba.m) <- c("Cell", "Mut", "Sites")

  fig_sites <- ggplot(nba.m, aes(Mut, Cell, fill = Sites)) +
    geom_tile(show.legend=T) +
    scale_fill_gradient(low = "white", high = "darkblue") +
    ylab(paste("Total sites: ", nrow(prob_mat), "cells")) +
    cardelino::heatmap.theme()# + cardelino::pub.theme()
  
  fig_sites
}


#' Plot confusion heatmap for assessing simulation accuracy
#' 
#' @param prob_mat A matrix (M x K), the estimated probability of cell j to 
#' clone k
#' @param prob_mat A matrix (M x K), the true assignment of cell to clone
#' @param prob_gap_min A float value, only show confusion matrix on cells with 
#' prob_gap >= prob_gap_min
#' 
#' @export
confusion_heatmap <- function(prob_mat, sim_mat, prob_gap_min=0.0,
                              pre_title=""){
  assign_0 <- cardelino::get_prob_label(sim_mat)
  assign_1 <- cardelino::get_prob_label(prob_mat)
  prob_gap <- cardelino::get_prob_gap(prob_mat)
  idx <- prob_gap >= prob_gap_min
  
  print(paste("assignable:", mean(idx)))
  print(paste("accuracy:", mean((assign_0 == assign_1)[idx])))
  print(paste("overall acc:", mean(assign_0 == assign_1)))
  
  acc = mean((assign_0 == assign_1)[idx])
  confusion_matrix <- as.data.frame(table(assign_0[idx], assign_1[idx]))
  colnames(confusion_matrix) <- c("Var1", "Var2", "Freq")
  
  confusion.plot <- ggplot(data=confusion_matrix, mapping=aes(x=Var1, y=Var2)) +
    geom_tile(aes(fill = Freq), colour = "grey") + 
    xlab("True clone") + ylab("Estimated clone") + 
    geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 0.5) + 
    ggtitle(paste0(pre_title, sprintf("Acc=%.1f%%", acc*100))) + 
    scale_fill_gradient(low = "white", high = "steelblue") +  
    theme_grey(base_size = 12) + pub.theme() + 
    theme(legend.position="none", 
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank())
  confusion.plot
}


#' The theme of heatmaps for prob_heatmap and sites_heatmap
#' 
#' @export
heatmap.theme <- function(legend.position="bottom") {
  ggplot2::theme_gray() + ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    legend.position=legend.position)
}


#' Plot a variant-cell heatmap for cell clonal assignment
#'
#' @param mat A matrix for heatmap: N variants x M cells. row and column will be 
#' sorted automatically.
#' @param prob A matrix of probability of clonal assignment: M cells x K clones
#' @param Config A binary matrix of clonal Configuration: N variants x K clones
#' @param show_legend A bool value: if TRUE, show the legend
#' 
#' @return a pheatmap object
#' 
#' @import pheatmap
#' @import ggplot2
#' 
#' @export
#' 
#' @references 
#' This function makes use of the \code{\link{pheatmap}} packages
#' 
#' @examples
vc_heatmap <- function(mat, prob, Config, show_legend=F){
  # sort variants
  mut_label <- Config %*% (2**seq(ncol(Config),1))
  mut_label <- seq(length(unique(mut_label)),1)[as.factor(mut_label)]
  idx_row <- order(mut_label - rowMeans(mat, na.rm = T)*0.9 + 0.05)
  anno_row <- data.frame(Mut=as.factor(mut_label))
  row.names(anno_row) <- row.names(Config)
  
  # sort cells
  cell_label <- cardelino::get_prob_label(prob)
  idx_col <- order(cell_label - diag(prob[, cell_label]))
  anno_col <- data.frame(Clone=as.factor(cell_label),
                         Prob=diag(prob[, cell_label]))
  row.names(anno_col) <- row.names(prob)
  
  # order the variant-cell matrix
  mat <- mat[idx_row, ][, idx_col]
  anno_row <- anno_row[idx_row, , F]
  anno_col <- anno_col[idx_col, , F]
  gaps_col <- cumsum(table(anno_col[[1]]))
  gaps_row <- cumsum(table(anno_row[[1]]))
  
  # pheatmap
  fig <- pheatmap::pheatmap(mat, legend = show_legend,
                            cluster_rows = F, cluster_cols = F, 
                            gaps_row = gaps_row, gaps_col = gaps_col,
                            annotation_row=anno_row, annotation_col = anno_col,
                            show_rownames = F, show_colnames = F)
  fig
}
