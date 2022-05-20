# Plotting functions

# Functions for heatmap plots

#' Plot a heatmap for probability of clone assignment
#'
#' @param prob_mat A matrix (M x K), the probability of cell j to clone k
#' @param threshold A float value, the threshold for assignable cells
#' @param mode A string, the mothod for defining scores for filtering cells:
#' best and delta. best: highest probability of a cell to K clones, delta: the
#' difference between the best and second.
#'
#' @param cell_idx A vector the indices of the input cells. If NULL, order by
#' the probabilty of each clone
#' @return a ggplot object
#'
#' @export
#'
#' @examples
#' data(example_donor)
#' assignments <- clone_id(A_clone, D_clone, Config = tree$Z, inference = "EM")
#' fig <- prob_heatmap(assignments$prob)
prob_heatmap <- function(prob_mat, threshold = 0.5, mode = "best",
                         cell_idx = NULL) {
    cell_label <- cardelino::rowArgmax(prob_mat)
    prob_value <- cardelino::rowMax(prob_mat, mode = mode)
    # add clone id
    colnames(prob_mat) <- paste0("C", seq_len(ncol(prob_mat)))
    for (i in seq_len(ncol(prob_mat))) {
        conf_frac <- mean(cell_label[prob_value >= threshold] == i)
        colnames(prob_mat)[i] <- paste0(
            "C", i, ": ",
            round(conf_frac * 100, digits = 1), "%"
        )
    }

    if (is.null(cell_idx)) {
        cell_idx <- order(cell_label - diag(prob_mat[, cell_label]))
    }
    nba.m <- data.frame(
        Cell = rep(seq_len(nrow(prob_mat)), ncol(prob_mat)),
        Clone = rep(colnames(prob_mat), each = nrow(prob_mat)),
        Prob = as.vector(prob_mat[cell_idx, ]),
        stringsAsFactors = FALSE
    )

    fig_assign <- ggplot(nba.m, aes_string("Clone", "Cell", fill = "Prob")) +
        geom_tile(show.legend = TRUE) +
        scale_fill_gradient(low = "white", high = "firebrick4") +
        ylab(paste("Clonal assignment:", nrow(prob_mat), "cells")) +
        heatmap.theme()

    fig_assign
}

#' Plot heatmap from a matrix
#'
#' @param mat A matrix to show, column by x-axis and row by y-axis
#' @param base_size Numeric value for the base size in theme_bw
#' @param digits Integer value for the number of digits to show
#' @param show_value Logical value for showing the value for each element or not
#' @export
heat_matrix <- function(mat, base_size = 12, digits = 2, show_value = FALSE) {
    df <- mtx_to_df(mat)
    if (!is.null(rownames(mat))) {
        df$Var1 <- factor(df$Var1, rownames(mat))
    }
    if (!is.null(colnames(mat))) {
        df$Var2 <- factor(df$Var2, colnames(mat))
    }

    df$value <- round(df$value, digits = digits)
    heat.plot <- ggplot(df, aes_string(x = "Var1", y = "Var2")) +
        geom_tile(aes_string(fill = "value"), colour = "grey") +
        scale_fill_gradient(low = "white", high = "steelblue") +
        # theme_grey(base_size = base_size) +
        theme_bw(base_size = base_size) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_x_discrete(expand = c(0, 0)) +
        theme(
            legend.position = "none",
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()
        )
    if (show_value) {
        heat.plot <- heat.plot +
            geom_text(aes_string(label = "value"),
                vjust = 0.5, size = base_size * 0.25
            )
    }
    heat.plot
}


#' The theme of heatmaps for prob_heatmap and sites_heatmap
#'
#' @param legend.position character, describes where to place legend on plot
#' (passed to \code{\link[ggplot2]{theme_gray}})
#' @param size numeric, base font size for plot (passed to
#' \code{\link[ggplot2]{theme_gray}})
#' @return a ggplot theme based on theme_gray
#'
heatmap.theme <- function(legend.position = "bottom", size = 12) {
    ggplot2::theme_gray(base_size = size) + ggplot2::theme(
        axis.text = ggplot2::element_text(size = size),
        axis.title = ggplot2::element_text(face = "bold", size = size),
        axis.title.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(
            face = "bold", size = size * 1.3,
            hjust = 0.5
        ),
        panel.grid.major = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        legend.position = legend.position,
        legend.title = ggplot2::element_text(size = size * 1.1)
    )
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
#' @return a ggplot object
#'
#' @export
#'
#' @examples
#' data(example_donor)
#' assignments <- clone_id(A_clone, D_clone, Config = tree$Z)
#' fig <- vc_heatmap(assignments$prob_variant, assignments$prob, tree$Z)
#' @references
#' This function makes use of the \code{\link{pheatmap}} packages
#'
vc_heatmap <- function(mat, prob, Config, show_legend = FALSE) {
    # sort variants
    mut_label <- Config %*% (2**seq(ncol(Config), 1))
    mut_label <- seq(length(unique(mut_label)), 1)[as.factor(mut_label)]
    idx_row <- order(mut_label - rowMeans(mat, na.rm = TRUE) * 0.9 + 0.05)
    anno_row <- data.frame(Mut = as.factor(mut_label))
    row.names(anno_row) <- row.names(Config)

    # sort cells
    cell_label <- cardelino::rowArgmax(prob)
    idx_col <- order(cell_label - diag(prob[, cell_label]))
    anno_col <- data.frame(
        Clone = as.factor(cell_label),
        Prob = diag(prob[, cell_label])
    )
    row.names(anno_col) <- row.names(prob)

    # order the variant-cell matrix
    mat <- mat[idx_row, ][, idx_col]
    anno_row <- anno_row[idx_row, , drop = FALSE]
    anno_col <- anno_col[idx_col, , drop = FALSE]
    gaps_col <- cumsum(table(anno_col[[1]]))
    gaps_row <- cumsum(table(anno_row[[1]]))

    # pheatmap
    fig <- pheatmap::pheatmap(mat,
        legend = show_legend,
        cluster_rows = FALSE, cluster_cols = FALSE,
        gaps_row = gaps_row, gaps_col = gaps_col,
        annotation_row = anno_row,
        annotation_col = anno_col,
        show_rownames = FALSE, show_colnames = FALSE
    )
    fig
}


## Functions for plotting phylogenetic trees

#' Plot a phylogenetic tree
#'
#' @param tree A phylgenetic tee object of class "phylo"
#' @param orient A string for the orientation of the tree: "v" (vertical) or "h"
#' (horizontal)
#'
#' @details This function plots a phylogenetic tree from an object of class
#' "phylo", as produced, for example, by the Canopy package.
#'
#' @return a ggtree object
#'
#' @importFrom ggtree ggtree geom_tiplab
#' @import ggplot2
#'
#' @author Davis McCarthy and Yuanhua Huang
#'
#' @export
#'
#' @references
#' This function makes use of the \code{\link{ggtree}} package:
#'
#' Guangchuang Yu, David Smith, Huachen Zhu, Yi Guan, Tommy Tsan-Yuk Lam.
#' ggtree: an R package for visualization and annotation of phylogenetic trees
#' with their covariates and other associated data. Methods in Ecology and
#' Evolution 2017, 8(1):28-36, doi:10.1111/2041-210X.12628
#'
#' @examples
#' data(example_donor)
#' plot_tree(tree, orient = "v")
plot_tree <- function(tree, orient = "h") {
    node_total <- max(tree$edge)
    node_shown <- ncol(tree$Z)
    node_hidden <- node_total - node_shown
    if (!is.null(tree$P)) {
        tree$tip.label[seq_len(node_shown)] <- paste0(
            "C", seq_len(node_shown),
            ": ", round(tree$P[, 1] * 100, digits = 0), "%"
        )
    }

    mut_ids <- 0
    mut_id_all <- tree$Z %*% (2**seq(ncol(tree$Z), 1))
    mut_id_all <- seq(length(unique(mut_id_all)), 1)[as.factor(mut_id_all)]

    branch_ids <- NULL
    for (i in seq_len(node_total)) {
        mut_num <- sum(tree$sna[, 3] == i)
        if (mut_num == 0) {
            if (i == node_shown + 1) {
                branch_ids <- c(branch_ids, "Root")
            } else {
                branch_ids <- c(branch_ids, "")
            }
        } else {
            mut_ids <- mut_ids + 1
            branch_ids <- c(
                branch_ids,
                paste0("M", mut_ids, ": ", mut_num, " SNVs")
            )
        }
    }
    pt <- ggtree::ggtree(tree)
    pt <- pt + ggplot2::geom_label(ggplot2::aes_string(x = "branch"),
        label = branch_ids, color = "firebrick"
    )
    pt <- pt + ggplot2::xlim(-0, node_hidden + 0.5) +
        ggplot2::ylim(0.8, node_shown + 0.5) # the degree may not be 3
    if (orient == "v") {
        pt <- pt + ggtree::geom_tiplab(hjust = 0.39, vjust = 1.0) +
            ggplot2::scale_x_reverse() + ggplot2::coord_flip()
    } else {
        pt <- pt + ggtree::geom_tiplab(hjust = 0.0, vjust = 0.5)
    }
    pt
}


#' Define a publication-style plot theme
#'
#' @param size numeric, base font size for adapted ggplot2 theme
#' @return a ggplot theme based on theme_classic
#'
#' @details This theme modifies the \code{\link[ggplot2]{theme_classic}} theme
#' in ggplot2.
#'
#' @import ggplot2
#' @export
#' @examples
#' library(ggplot2)
#' x <- sample(10)
#' y <- x + runif(10) - 0.5
#' df <- data.frame(x = x, y = y)
#' fig <- ggplot(df, aes(x = x, y = y)) +
#'     geom_point() +
#'     pub.theme()
pub.theme <- function(size = 12) {
    theme_classic(base_size = size) +
        ggplot2::theme(
            axis.text = ggplot2::element_text(size = size),
            axis.title = ggplot2::element_text(
                face = "bold", size = size
            ),
            plot.title = ggplot2::element_text(
                size = size * 1.2, hjust = 0.5
            ),
            legend.title = ggplot2::element_text(size = size * 1.1),
            legend.text = ggplot2::element_text(size = size),
            panel.grid.major = ggplot2::element_line(
                size = 0.1, colour = "#d3d3d3"
            ),
            panel.grid.minor = ggplot2::element_line(
                size = 0.05, colour = "#d3d3d3"
            )
        )
}


#' Define a publication-style plot theme
#'
#' @param Config1 variant by clone matrix defining the first clonal structure
#' @param Config2 variant by clone matrix defining the second clonal structure
#' @param show_variant_names logical(1), should the variant names (rownames of
#' Config matrices) be shown on the plot? Default is \code{FALSE}.
#'
#' @return a ggplot heatmap style plot showing the differences between the two
#' Config matrices, specifically the differences Config1 - Config2.
#'
#' @import ggplot2
#' @export
#' @examples
#' Config1 <- matrix(c(
#'     rep(0, 15), rep(1, 8), rep(0, 7),
#'     rep(1, 5), rep(0, 3), rep(1, 7)
#' ), ncol = 3)
#' Config2 <- matrix(c(
#'     rep(0, 15), rep(1, 8), rep(1, 7),
#'     rep(0, 5), rep(1, 3), rep(1, 7)
#' ), ncol = 3)
#' rownames(Config1) <- rownames(Config2) <- paste0("var", 1:nrow(Config1))
#' colnames(Config1) <- colnames(Config2) <- paste0("clone", 1:ncol(Config1))
#' plot_config_diffs(Config1, Config2)
plot_config_diffs <- function(Config1, Config2, show_variant_names = FALSE) {
    if (!identical(rownames(Config1), rownames(Config2))) {
          stop("Config matrices must have identical rownames.")
      }
    if (!identical(colnames(Config1), colnames(Config2))) {
          stop("Config matrices must have identical colnames.")
      }
    diffs <- Config1 - Config2
    df <- data.frame(
        variant = factor(rep(rownames(Config1), ncol(Config1)),
            levels = rev(rownames(Config1))
        ),
        clone = rep(colnames(Config1), each = nrow(Config1)),
        diffs = as.vector(diffs)
    )
    p_out <- ggplot(df,
                    aes_string(x = "clone", y = "variant", fill = "diffs")) +
        geom_raster() +
        scale_x_discrete(position = "top") +
        scale_fill_gradient2(
            low = "dodgerblue4", mid = "white",
            high = "firebrick4", midpoint = 0, space = "Lab",
            na.value = "grey50", guide = "colourbar",
            aesthetics = "fill", name = "Config\ndifferences",
            limits = c(-1, 1)
        ) +
        heatmap.theme() +
        theme(
            legend.position = "bottom",
            legend.key.width = unit(0.5, "in")
        )
    if (!show_variant_names) {
          p_out <- p_out + theme(
              axis.text.y = element_blank(),
              axis.title.y = element_blank()
          )
      }
    p_out
}

# ## joxm - works
# card_joxm <- readRDS(paste0("../fibroblast-clonality/data/cell_assignment/",
#   "cardelino_results_carderelax.joxm.filt_lenient.cell_coverage_sites.rds"))
# names(card_joxm)
# card_joxm_Config_prob <- card_joxm$tree$Z
# card_joxm_Config_prob[,] <- colMeans(card_joxm$Config_all)
# card_joxm_Config_best <- round(card_joxm_Config_prob)
#
# plot_config_diffs(card_joxm_Config_prob, card_joxm$tree$Z) + ggtitle("joxm")
# plot_config_diffs(card_lexy_Config_prob, card_lexy$tree$Z)
# plot_config_diffs(card_zoxy_Config_prob, card_zoxy$tree$Z)
