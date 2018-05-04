## Functions for plotting phylogenetic trees

#' Plot a phylogenetic tree
#'
#' @param tree A phylgenetic tee object of class "phylo"
#' @param orient A string for the orientation of the tree: "v" (vertical) or "h"
#' (horizontal)
#' 
#' @details This function plots a phylogenetic tree from an object of class "phylo", as
#' produced, for example, by the Canopy package. 
#' 
#' @return a ggtree object
#' 
#' @import ggtree
#' @import ggplot2
#' 
#' @author Davis McCarthy and Yuanhua Huang
#' 
#' @export
#' 
#' @references 
#' This function makes use of the \code{\link{ggtree}} packages: 
#' 
#' Guangchuang Yu, David Smith, Huachen Zhu, Yi Guan, Tommy Tsan-Yuk Lam. 
#' ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. Methods in Ecology and Evolution 2017, 8(1):28-36, doi:10.1111/2041-210X.12628
#' 
#' @examples 
plot_tree <- function(tree, orient="h") {
  node_total <- max(tree$edge)
  node_shown <- length(tree$P[, 1])
  node_hidden <- node_total - node_shown
  
  prevalence <- c(tree$P[, 1]*100, rep(0, node_hidden))
  # node_size <- c(rep(20, node_shown), rep(0, node_hidden))
  
  mut_ids <- 0
  mut_id_all <- tree$Z %*% (2**seq(ncol(tree$Z),1))
  mut_id_all <- seq(length(unique(mut_id_all)),1)[as.factor(mut_id_all)]
  
  branch_ids <- NULL
  for (i in seq_len(node_total)) {
    if (i <= node_shown) {
      tree$tip.label[i] = paste0("C", i, ": ", round(prevalence[i], digits = 0),
                                 "%")
    }
    mut_num = sum(tree$sna[,3] == i)
    if (mut_num == 0) {
      if(i == node_shown+1){branch_ids = c(branch_ids, "Root")}
      else{branch_ids = c(branch_ids, "")} #NA
    }
    else {
      vaf <- mean(tree$VAF[tree$sna[,3] == i])
      mut_ids <- mut_ids + 1
      mut_ids <- mean(mut_id_all[tree$sna[,3] == i])
      branch_ids = c(branch_ids, paste0("Mut", mut_ids, ": ", mut_num, 
                                        sprintf("; %.2f", vaf)))
    }
  }
  pt <- ggtree::ggtree(tree)
  pt <- pt + ggplot2::geom_label(ggplot2::aes_string(x = "branch"), label = branch_ids, 
                                 color = "firebrick")
  pt <- pt + ggplot2::xlim(-0, node_hidden + 0.5) + ggplot2::ylim(0.8, node_shown + 0.5) #the degree may not be 3
  if (orient == "v") {
    pt <- pt + ggtree::geom_tiplab(hjust = 0.39, vjust = 1.0) + 
        ggplot2::scale_x_reverse() + ggplot2::coord_flip() 
  } else {
    pt <- pt + ggtree::geom_tiplab(hjust = 0.0, vjust = 0.5)
  }
  pt
}

mut.label <- function(tree){
  SNA.label <- tree$sna[, 3]
  mut_ids <- 0
  branch_ids <- NULL
  for (i in seq_len(max(tree$edge))) {
    mut_num = sum(tree$sna[,3] == i)
    if (mut_num == 0) {
      branch_ids = c(branch_ids, "") #NA
    }
    else{
      mut_ids <- mut_ids + 1
      branch_ids = c(branch_ids, paste0("Mut", mut_ids, ": ", mut_num))
      SNA.label[tree$sna[,3] == i] = paste0("Mut", mut_ids, ": ", mut_num)
    }
  }
  SNA.label
}


#' @export
pub.theme <- function(size = 12){
    ggplot2::theme(axis.text = ggplot2::element_text(size = size),
                   axis.title = ggplot2::element_text(
                       face = "bold", size = size),
                   plot.title = ggplot2::element_text(
                       face = "bold", size = size * 1.3, hjust = 0.5))
}
