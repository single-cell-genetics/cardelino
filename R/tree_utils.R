#' Get a clonal tree from a configuration matrix
#'
#' @param Config variant x clone matrix of binary values. The clone-variant
#' configuration, which encodes the phylogenetic tree structure. This is the
#' output Z of Canopy
#' @param P a one-column numeric matrix encoding the (observed or estimated)
#' prevalence (or frequency) of each clone
#' @param strictness character(1), a character string defining the strictness of
#' the function if there are all-zero rows in the Config matrix. If \code{"lax"}
#' then the function silently drops all-zero rows and proceeds. If \code{"warn"}
#' then the function warns of dropping all-zero rows and proceeds. If
#' \code{"error"} then the function throws an error is all-zero rows are
#' detected.
#'
#' @return
#' An object of class "phylo" describing the tree structure. The output object
#' also contains an element "sna" defining the clustering of variants onto the
#' branches of the tree, and if \code{P} is non-null it also contains VAF
#' (variant allele frequency), CCF (cell clone fraction) and clone prevalence
#' values (computed from the supplied \code{P} argument).
#'
#' @details
#' Output tree may be nonsensical if the input \code{Config} matrix does not
#' define a coherent tree structure.
#'
#' @author Davis McCarthy
#'
#' @import utils
#' @export
#'
#' @examples
#' Configk3 <- matrix(c(rep(0, 15), rep(1, 8), rep(0, 7), rep(1, 5), rep(0, 3), rep(1, 7)), ncol = 3)
#' tree_k3 <- get_tree(Config = Configk3, P = matrix(rep(1/3, 3), ncol = 1))
#' plot_tree(tree_k3)
get_tree <- function(Config, P = NULL, strictness = "lax") {
  if (!is.null(P)) {
    if (ncol(P) != 1)
      stop("P must be a matrix with one column encoding clone prevalence values")
  }
  all_zero_rows <- rowSums(Config) == 0
  strictness <- match.arg(strictness, c("lax", "warn", "error"))
  if (any(all_zero_rows)) {
    if (strictness == "error")
      stop("Config matrix contains all-zero rows.")
    else {
      if (strictness == "warn")
        warning("Dropped ", sum(all_zero_rows), " all-zero rows from Config matrix.")
      else
        message("Dropped ", sum(all_zero_rows), " all-zero rows from Config matrix.")
      Config <- Config[!all_zero_rows,]
    }
  }
  k <- ncol(Config) # number of clones
  varnames <- rownames(Config)
  sna <- matrix(nrow = nrow(Config), ncol = 3)
  sna[, 1] <- seq_len(nrow(sna))
  rownames(sna) <- varnames
  colnames(sna) <- c("sna", "sna.st.node", "sna.ed.node")
  tip_label <- seq_len(k)
  tip_vals <- 2^seq_len(k)
  ## Need to determine number of internal nodes in the tree
  Config2 <- t(t(Config) * tip_vals)
  var_bin_vals <- rowSums(Config2)
  node_vals <- names(table(var_bin_vals))
  node_vals <- node_vals[!(node_vals %in% tip_vals)]
  node_num <- sum(!(node_vals %in% tip_vals))
  ## define a list with subsets of edge matrices
  ## start with root node (k+1), which always connects to tip 1 (base clone)
  if (node_num > 0.5) {
    node_def_list <- list()
    edge_list <- list()
    for (i in seq_len(k - 1)) {
      clone_combos <- utils::combn(2:k, (k - i), simplify = FALSE)
      for (j in seq_len(length(clone_combos))) {
        test_sum <- sum(tip_vals[clone_combos[[j]]])
        if (test_sum %in% node_vals)
          node_def_list[[
            paste0("node", paste0(clone_combos[[j]]), collapse = "_")]] <-
            clone_combos[[j]]
      }
    }
    if (node_num != length(node_def_list))
      stop("Conflict in computed number of internal nodes.")
    ## Sort out edges for the root node
    tip_nodes <- seq_len(k)
    root_to_tip <- tip_nodes[!(tip_nodes %in% unique(unlist(node_def_list)))]
    edge_list[["root_node_tips"]] <- matrix(
      c(rep(k + 1, length(root_to_tip)), root_to_tip),
      nrow = length(root_to_tip))
    el_counter <- 1
    for (i in seq_len(length(node_def_list))) {
      ## add edge from root to internal node if not already done
      if (i < 1.5) {
        el_counter <- el_counter + 1
        edge_list[[el_counter]] <- matrix(c(k + 1, k + 1 + i), nrow = 1)
        sna[var_bin_vals == sum(2^node_def_list[[i]]), 2] <- k + 1
        sna[var_bin_vals == sum(2^node_def_list[[i]]), 3] <- k + 1 + i
      } else {
        clones_in_this_node <- node_def_list[[i]]
        clones_in_prev_nodes <- unique(unlist(node_def_list[seq_len(i - 1)]))
        if (!any(clones_in_this_node %in% clones_in_prev_nodes)) {
          el_counter <- el_counter + 1
          edge_list[[el_counter]] <- matrix(c(k + 1, k + 1 + i), nrow = 1)
          sna[var_bin_vals == sum(2^node_def_list[[i]]), 2] <- k + 1
          sna[var_bin_vals == sum(2^node_def_list[[i]]), 3] <- k + 1 + i
        }
      }
      ## add edge from internal node to internal node
      ## if all of the clones for the node are present in the previous node in
      ## the tree (immediately above in the hierarchy), then add the edge
      ## check the size of previous nodes, and select the node that has minimum
      ## number of clones that is more than the number in this node
      if (i > 1.5) {
        prev_nodes <- seq_len(i - 1)
        prev_node_sizes <- vapply(node_def_list[prev_nodes], length, numeric(1))
        prev_nodes <- prev_nodes[prev_node_sizes > length(node_def_list[[i]])]
        min_prev_node_size <- min(prev_node_sizes[prev_nodes])
        prev_nodes <- prev_nodes[prev_node_sizes[prev_nodes] ==
                                   min_prev_node_size]
        for (j in prev_nodes) {
          if (all(node_def_list[[i]] %in% node_def_list[[j]])) {
            el_counter <- el_counter + 1
            edge_list[[el_counter]] <- matrix(
              c(k + 1 + j, k + 1 + i), nrow = 1)
            sna[var_bin_vals == sum(2^node_def_list[[i]]), 2] <- k + 1 + j
            sna[var_bin_vals == sum(2^node_def_list[[i]]), 3] <- k + 1 + i
          }
        }
      }
      ## add edge from internal node to tip
      ## (if clone not present in any subsequent nodes)
      if (node_num < 1.5) {
        ## if only one internal node, there are edges from this node to all tips
        node_to_tip <- tip_nodes[-1]
        el_counter <- el_counter + 1
        edge_list[[el_counter]] <- matrix(
          c(rep(k + 1 + i, length(node_to_tip)), node_to_tip),
          nrow = length(node_to_tip))
        for (m in node_to_tip) {
          sna[var_bin_vals == sum(2^m), 2] <- k + 1 + i
          sna[var_bin_vals == sum(2^m), 3] <- m
        }
      } else {
        ## if more than one internal node, need to check if tips mentioned in
        ## this node appear in any subsequent nodes
        node_to_tip <- node_def_list[[i]]
        if (i < node_num) {
          node_to_tip <- node_to_tip[
            !(node_to_tip %in% unique(unlist(node_def_list[(i + 1):node_num])))]
        } ## else this is the last node; just connect edges from node to tips
        if (length(node_to_tip) > 0.5) {
          el_counter <- el_counter + 1
          edge_list[[el_counter]] <- matrix(
            c(rep(k + 1 + i, length(node_to_tip)), node_to_tip),
            nrow = length(node_to_tip))
          for (m in node_to_tip) {
            sna[var_bin_vals == sum(2^m), 2] <- k + 1 + i
            sna[var_bin_vals == sum(2^m), 3] <- m
          }
        }
      }
    }
  } else {
    edge_list <- list("root_node" = matrix(c(rep(k + 1, k), seq_len(k)),
                                           ncol = 2))
    for (j in 2:k) {
      sna[var_bin_vals == 2^j, 2] <- k + 1
      sna[var_bin_vals == 2^j, 3] <- j
    }
  }
  # node_def_list
  edge_mat <- do.call(rbind, edge_list)
  tree_out <- list(edge = edge_mat, Nnode = node_num + 1, tip.label = tip_label)
  class(tree_out) <- "phylo"
  tree_out$Z <- Config
  if (!is.null(P)) {
    tree_out$P <- P
    tree_out$VAF <- tree_out$Z %*% tree_out$P / 2
    tree_out$CCF <- tree_out$Z %*% tree_out$P
  }
  tree_out$sna <- sna
  tree_out
}


### Create a tree from a config matrix
# Configk3 <- matrix(c(rep(0, 15), rep(1, 8), rep(0, 7), rep(1, 5), rep(0, 3),
#                      rep(1, 7)), ncol = 3)
# Configk4_1 <- matrix(c(rep(0, 15), rep(1, 8), rep(0, 7), rep(1, 5), rep(0, 3),
#                        rep(1, 4), rep(0, 3), rep(1, 5), rep(0, 7), rep(1, 3)),
#                      ncol = 4)
# Configk4_2 <- matrix(c(rep(0, 15), rep(1, 8), rep(0, 7), rep(1, 5), rep(0, 3),
#                        rep(1, 4), rep(0, 3), rep(1, 5), rep(0, 3), rep(1, 2),
#                        rep(0, 2), rep(1, 3)),
#                      ncol = 4)
# Configk4_2_bad <- matrix(c(rep(0, 15), rep(1, 8), rep(0, 2), rep(1, 2),
#                            rep(0, 3), rep(1, 5), rep(0, 3),
#                            rep(1, 4), rep(0, 3), rep(1, 5), rep(0, 3), rep(1, 2),
#                            rep(0, 2), rep(1, 3)),
#                          ncol = 4)
#
# Configk5 <- matrix(c(rep(0, 15), rep(1, 8), rep(0, 7), rep(1, 5), rep(0, 3),
#                      rep(1, 4), rep(0, 3), rep(1, 5), rep(0, 3), rep(1, 2),
#                      rep(0, 2), rep(0, 3), rep(1, 8), rep(0, 4), rep(1, 3)),
#                    ncol = 5)
#
# tree_k3 <- get_tree(Config = Configk3, P = matrix(rep(1/3, 3), ncol = 1))
# tree_k3$sna
# tree_k3$Z
# tree_k3$edge
# cardelino::plot_tree(tree_k3)
#
# tree_k4_1 <- get_tree(Config = Configk4_1, P = matrix(rep(1/4, 4), ncol = 1))
# tree_k4_1$sna
# tree_k4_1$Z
# tree_k4_1$edge
# plot_tree(tree_k4_1)
#
# tree_k4_2 <- get_tree(Config = Configk4_2, P = matrix(rep(1/4, 4), ncol = 1))
# tree_k4_2$sna
# tree_k4_2$Z
# tree_k4_2$edge
# plot_tree(tree_k4_2)
#
# tree_k4_2_bad <- get_tree(Config = Configk4_2_bad, P = matrix(rep(1/4, 4), ncol = 1))
# tree_k4_2_bad$sna
# tree_k4_2_bad$Z
# tree_k4_2_bad$edge
# plot_tree(tree_k4_2_bad)
#
# tree_k5 <- get_tree(Config = Configk5, P = matrix(rep(1/5, 5), ncol = 1))
# tree_k5$sna
# tree_k5$Z
# tree_k5$edge
# plot_tree(tree_k5)
#
# ## joxm - works
# card_joxm <- readRDS("../fibroblast-clonality/data/cell_assignment/cardelino_results_carderelax.joxm.filt_lenient.cell_coverage_sites.rds")
# names(card_joxm)
# card_joxm_Config_prob <- card_joxm$tree$Z
# card_joxm_Config_prob[,] <- colMeans(card_joxm$Config_all)
# card_joxm_Config_best <- round(card_joxm_Config_prob)
#
# tree_joxm <- get_tree(card_joxm_Config_best,
#                       P = matrix(colMeans(card_joxm$prob_mat > 0.5), ncol = 1))
# p1 <- cardelino::plot_tree(tree_joxm, orient = "v") +
#   ggtitle("joxm: Cardelino output tree") +
#   theme(plot.title = element_text(hjust = 0.5))
# p2 <- cardelino::plot_tree(card_joxm$tree, orient = "v") +
#   ggtitle("joxm: Canopy tree") +
#   theme(plot.title = element_text(hjust = 0.5))
# cowplot::plot_grid(p2, p1, ncol = 1)
#
# ## zoxy - fails
# card_zoxy <- readRDS("../fibroblast-clonality/data/cell_assignment/cardelino_results_carderelax.zoxy.filt_lenient.cell_coverage_sites.rds")
# names(card_zoxy)
# card_zoxy_Config_prob <- card_zoxy$tree$Z
# card_zoxy_Config_prob[,] <- colMeans(card_zoxy$Config_all)
# card_zoxy_Config_best <- round(card_zoxy_Config_prob)
#
# tree_zoxy <- get_tree(card_zoxy_Config_best,
#                       P = matrix(colMeans(card_zoxy$prob_mat > 0.5), ncol = 1))
# tree_zoxy$edge
# tree_zoxy$sna
# tree_zoxy <- get_tree(card_zoxy_Config_best,
#                       P = matrix(colMeans(card_zoxy$prob_mat > 0.5), ncol = 1),
#                       strictness = "warn")
# tree_zoxy <- get_tree(card_zoxy_Config_best,
#                       P = matrix(colMeans(card_zoxy$prob_mat > 0.5), ncol = 1),
#                       strictness = "e")
# ## handles all-zero rows
# ## looks like a problem caused by no shared variants between cl2 and cl3
# p1 <- cardelino::plot_tree(tree_zoxy, orient = "v") +
#   ggtitle("zoxy: Cardelino output tree") +
#   theme(plot.title = element_text(hjust = 0.5))
# p2 <- cardelino::plot_tree(card_zoxy$tree, orient = "v") +
#   ggtitle("zoxy: Canopy tree") +
#   theme(plot.title = element_text(hjust = 0.5))
# cowplot::plot_grid(p2, p1, ncol = 1)
#
# ## lexy - some all-zero rows
# card_lexy <- readRDS("../fibroblast-clonality/data/cell_assignment/cardelino_results_carderelax.lexy.filt_lenient.cell_coverage_sites.rds")
# names(card_lexy)
# card_lexy_Config_prob <- card_lexy$tree$Z
# card_lexy_Config_prob[,] <- colMeans(card_lexy$Config_all)
# card_lexy_Config_best <- round(card_lexy_Config_prob)
#
# tree_lexy <- get_tree(card_lexy_Config_best,
#                       P = matrix(colMeans(card_lexy$prob_mat > 0.5), ncol = 1))
# tree_lexy$edge
# ##
# p1 <- cardelino::plot_tree(tree_lexy, orient = "v") +
#   ggtitle("lexy: Cardelino output tree") +
#   theme(plot.title = element_text(hjust = 0.5))
# p2 <- cardelino::plot_tree(card_lexy$tree, orient = "v") +
#   ggtitle("lexy: Canopy tree") +
#   theme(plot.title = element_text(hjust = 0.5))
# cowplot::plot_grid(p2, p1, ncol = 1)
#
# ## vils - works
# card_vils <- readRDS("../fibroblast-clonality/data/cell_assignment/cardelino_results_carderelax.vils.filt_lenient.cell_coverage_sites.rds")
# names(card_vils)
# card_vils_Config_prob <- card_vils$tree$Z
# card_vils_Config_prob[,] <- colMeans(card_vils$Config_all)
# card_vils_Config_best <- round(card_vils_Config_prob)
#
# tree_vils <- get_tree(
#   card_vils_Config_best,
#   P = matrix(colSums(card_vils$prob_mat > 0.5) /
#                sum(rowMax(card_vils$prob_mat) > 0.5), ncol = 1))
# tree_vils$edge
# tree_vils$P
# ##
# p1 <- cardelino::plot_tree(tree_vils, orient = "v") +
#   ggtitle("vils: Cardelino output tree") +
#   theme(plot.title = element_text(hjust = 0.5))
# p2 <- cardelino::plot_tree(card_vils$tree, orient = "v") +
#   ggtitle("vils: Canopy tree") +
#   theme(plot.title = element_text(hjust = 0.5))
# cowplot::plot_grid(p2, p1, ncol = 1)
#
# ## sehl - fails
# card_sehl <- readRDS("../fibroblast-clonality/data/cell_assignment/cardelino_results_carderelax.sehl.filt_lenient.cell_coverage_sites.rds")
# names(card_sehl)
# card_sehl_Config_prob <- card_sehl$tree$Z
# card_sehl_Config_prob[,] <- colMeans(card_sehl$Config_all)
# card_sehl_Config_best <- (card_sehl_Config_prob > 0.49) * 1L
#
# tree_sehl <- get_tree(
#   card_sehl_Config_best,
#   P = matrix(colSums(card_sehl$prob_mat > 0.5) /
#                sum(rowMax(card_sehl$prob_mat) > 0.5), ncol = 1))
# tree_sehl$edge
# tree_sehl$P
# ## Too many edges - looks like it finds an internal node that does not make sense
# p1 <- cardelino::plot_tree(tree_sehl, orient = "v") +
#   ggtitle("sehl: Cardelino output tree") +
#   theme(plot.title = element_text(hjust = 0.5))
# p2 <- cardelino::plot_tree(card_sehl$tree, orient = "v") +
#   ggtitle("sehl: Canopy tree") +
#   theme(plot.title = element_text(hjust = 0.5))
# cowplot::plot_grid(p2, p1, ncol = 1)


#
#
# ## example data - fails; finds too many internal nodes
# data(example_donor)
# plot_tree(get_tree(tree$Z))
# assignments_binom <- cell_assign_Gibbs(A_clone, D_clone, Config = tree$Z,
#                                        relax_Config = TRUE, model = "binomial")
# eg_Config_prob <- assignments_binom$Config_prob
# eg_Config_best <- round(eg_Config_prob)
#
# tree_eg <- get_tree(eg_Config_best,
#                       P = matrix(colMeans(assignments_binom$prob > 0.5), ncol = 1))
# tree_eg$edge
# ## too many internal nodes
# ggtree::ggtree(tree_eg)
# plot_tree(tree_eg)
# ## gives a new and different error *shrug*
