# get_clone_score <- function(x, var_meta, clone = "clone1") {
#     locs <- dplyr::filter(var_meta, clone == clone)[["locName"]]
#     score <- rowMeans(as.matrix(x[, locNames(x) %in% locs]) > 0, na.rm = TRUE)
#     score[is.nan(score)] <- 0
#     score[is.na(score)] <- 0
#     score
# }
# 
# compute_relatedness_vec <- function(geno_vec, geno_mat, vaf = NULL) {
#     ## geno_mat should be samples x variants
#     if ( !is.null(vaf) ) {
#         ## standardise genotypes
#         message("Standardise genotypes\n")
#         stand_geno_vec <- ( (geno_vec - 2 * vaf)
#             / sqrt(2 * vaf * (1 - vaf)) )
#         stand_geno_mat<- t( (t(geno_mat) - 2 * vaf)
#                                  / sqrt(2 * vaf * (1 - vaf)))
#         message("Compute genotypic correlations\n")
#         corr_mat <- t(t(stand_geno_mat) * stand_geno_vec)
#         #bad_snp <- apply(corr_mat, 2, function(x) any(is.na(x)))
#         #corr_mat <- corr_mat[, !bad_snp, drop = FALSE]
#         ## sum correlations across donors
#         ave_allelic_corr <- rowMeans(corr_mat, na.rm = TRUE)
#     } else {
#         ave_allelic_corr <- rep(NA, length(geno_vec))
#     }
#     ave_allelic_corr
# }
# 
# simil_realrel <- function(geno_mat, vaf) {
#     ## geno_mat should be samples x variants
#     simil_out <- apply(geno_mat, 1, compute_relatedness_vec, geno_mat, vaf)
#     simil_out
# }
# 
# 
# 
# dist_shared_alleles_vec <- function(geno_vec, geno_mat) {
#     ## compute distance between a sample and a set of samples
#     ## based on alternative allele count
#     raw_mat <-  abs(t(t(geno_mat) - geno_vec)) / 2
#     ave_allelic_dist <- rowMeans(raw_mat, na.rm = TRUE)
# }
# 
# dist_shared_alleles <- function(geno_mat) {
#     dist_out <- apply(geno_mat, 1, dist_shared_alleles_vec, geno_mat)
#     as.dist(dist_out)
# }
# 
# 
# plot_sil_width <- function(sill) {
#     df <- data_frame(cluster = as.factor(sill[, 1]),
#                      neighbor = as.factor(sill[, 2]),
#                      sil_width = sill[, 3])
#     df <- dplyr::arrange(df, cluster, sil_width) %>%
#         dplyr::mutate(y = 1:nrow(.))
#     ave_wid <- mean(df[["sil_width"]])
#     df_ave <- group_by(df, cluster) %>%
#         summarise(ave_sil_wid = mean(sil_width),
#                   y_pos = mean(y)) %>%
#         dplyr::mutate(
#             lab = paste0("Ave. ", format(ave_sil_wid, digits = 2)))
#     ggplot(df, aes(y = y, x = sil_width, colour = cluster)) +
#         geom_segment(aes(xend = 0, yend = y), size = 0.5) +
#         geom_text(aes(x = -0.7, y = y_pos, label = lab), data = df_ave,
#                   show.legend = FALSE) +
#         geom_vline(xintercept = 0, linetype = 2) +
#         scale_color_tableau(palette = "tableau20") +
#         xlab("Silhouette width") +
#         theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
#               axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
#         ggtitle(paste0("Ave. silh. width: ", format(ave_wid, digits = 3)))
# }
# 
# ggmds <- function(dist_mat, df, colour = "cluster") {
#     mds_out <- cmdscale(dist_mat, k = 6)
#     df <- bind_cols(
#         df, data_frame(MDS_1 = mds_out[, 1], MDS_2 = mds_out[, 2],
#                        MDS_3 = mds_out[, 3], MDS_4 = mds_out[, 4],
#                        MDS_5 = mds_out[, 5], MDS_6 = mds_out[, 6]))
#     p <- ggplot(df, 
#            aes_string(x = "MDS_1", y = "MDS_2", colour = colour)) +
#         geom_point(alpha = 0.5) + theme(legend.position = "bottom") +
#         guides(size = FALSE)
#     print(p)
#     df
# }
# 
# 
# gensimil <- function(dat, lambda = mean(dat > 0, na.rm = TRUE), delta = 1,
#                      nan_replace = 0.001) {
#     bindat <- (dat > 0)
#     m11 <- (replace(bindat, is.na(bindat), 0) %*%
#             t(replace(bindat, is.na(bindat), 0)))
#     m00 <- (replace(!bindat, is.na(bindat), 0) %*%
#             t(replace(!bindat, is.na(bindat), 0)))
#     m01 <- (replace(!bindat, is.na(bindat), 0) %*%
#             t(replace(bindat, is.na(bindat), 0)))
#     m10 <- (replace(bindat, is.na(bindat), 0) %*%
#             t(replace(!bindat, is.na(bindat), 0)))
#     sim <- ( (m11 + lambda * m00) /
#              (m11 + lambda * m00 + delta * (m10 + m01)) )
#     sim[is.nan(sim)] <- nan_replace
#     diag(sim) <- 1
#     sim
# }
# 
# get_ap_clusters <- function(ap) {
#     out <- as.factor(ap@idx)
#     for (i in seq_along(levels(out)))
#          levels(out) <- gsub(levels(out)[i], i, levels(out))
#     out
# }
# 
# 
# canopy_post <- function(sampchain, projectname, K, numchain, burnin, thin, 
#                         optK, C = NULL, post.config.cutoff = NULL) {
#     if (is.null(C)) {
#         C = diag(nrow(sampchain[[1]][[1]][[1]]$cna))
#         colnames(C) = rownames(C) = rownames(sampchain[[1]][[1]][[1]]$cna)
#     }
#     if (is.null(post.config.cutoff)) {
#         post.config.cutoff = 0.05
#     }
#     if (post.config.cutoff > 1 | post.config.cutoff <= 0) {
#         stop("post.config.cutoff has to be between 0 and 1!")
#     }
#     sampchaink = sampchain[[which(K == optK)]]
#     numchain = length(sampchaink)
#     samptreenew = sampchaink[[1]][(burnin + 1):length(sampchaink[[1]])]
#     numpostburn = length(samptreenew)
#     temp <- thin * c(1:(numpostburn/thin))
#     samptreethin = samptreenew[temp]
#     length(samptreethin)
#     for (numi in 2:numchain) {
#         samptreenew = sampchaink[[numi]][(burnin + 1):length(sampchaink[[numi]])]
#         numpostburn = length(samptreenew)
#         temp <- thin * c(1:(numpostburn/thin))
#         samptreethin = c(samptreethin, samptreenew[temp])
#     }
#     samptreethin.lik = rep(NA, length(samptreethin))
#     for (treei in 1:length(samptreethin)) {
#         samptreethin.lik[treei] = samptreethin[[treei]]$likelihood
#     }
#     samptreethin = samptreethin[which((rank(-1 * samptreethin.lik, 
#                                             ties.method = "first")) <= 5 * (length(samptreethin)/numchain))]
#     samptreethin.lik = rep(NA, length(samptreethin))
#     for (treei in 1:length(samptreethin)) {
#         samptreethin.lik[treei] = samptreethin[[treei]]$likelihood
#     }
#     if (!is.null(sampchain[[1]][[1]][[1]]$cna)) {
#         for (i in 1:length(samptreethin)) {
#             samptreethin[[i]] = sortcna(samptreethin[[i]], C)
#         }
#     }
#     for (i in 1:length(samptreethin)) {
#         samptreethin[[i]]$clonalmut = getclonalcomposition(samptreethin[[i]])
#     }
#     config = rep(NA, length(samptreethin))
#     config[1] = 1
#     categ = 1
#     for (i in 2:length(samptreethin)) {
#         for (categi in 1:categ) {
#             list.a = samptreethin[[i]]$clonalmut
#             list.b = samptreethin[[which(config == categi)[1]]]$clonalmut
#             if ((sum(is.element(list.a, list.b)) == optK) &
#                 (sum(is.element(list.b, list.a)) == optK)) {
#                 config[i] = categi
#             }
#         }
#         if (is.na(config[i])) {
#             config[i] = categ + 1
#             categ = categ + 1
#         }
#     }
#     z.temp = (samptreethin.lik - mean(samptreethin.lik))/sd(samptreethin.lik)
#     samptreethin = samptreethin[z.temp <= 1.5 & z.temp >= -1.5]
#     samptreethin.lik = samptreethin.lik[z.temp <= 1.5 & z.temp >= 
#                                         -1.5]
#     config = config[z.temp <= 1.5 & z.temp >= -1.5]
#     config.summary = matrix(nrow = length(unique(config)), ncol = 3)
#     colnames(config.summary) = c("Configuration", "Post_prob", 
#                                  "Mean_post_lik")
#     config.summary[, 1] = unique(config)
#     for (i in 1:nrow(config.summary)) {
#         configi = config.summary[i, 1]
#         configi.temp = which(config == configi)
#         config.summary[i, 2] = round(length(configi.temp)/length(config), 
#                                      3)
#         config.summary[i, 3] = round(max(samptreethin.lik[which(config == 
#                                                                 configi)]), 2)
#     }
#     minor.config = which(config.summary[, 2] < post.config.cutoff)
#     if (length(minor.config) > 0) {
#         config.sel = rep(TRUE, length(config))
#         for (i in minor.config) {
#             config.sel[which(config == config.summary[i, 1])] = FALSE
#         }
#         samptreethin = samptreethin[config.sel]
#         samptreethin.lik = samptreethin.lik[config.sel]
#         config = config[config.sel]
#         config.summary = config.summary[-minor.config, , drop = FALSE]
#         for (i in 1:nrow(config.summary)) {
#             config[which(config == config.summary[i, 1])] = i
#         }
#         config.summary[, 1] = 1:nrow(config.summary)
#         config.summary[, 2] = round(config.summary[, 2]/sum(config.summary[, 
#                                                                            2]), 3)
#     }
#     for (treei in 1:length(samptreethin)) {
#         output.tree = samptreethin[[treei]]
#         output.tree.Z = output.tree$Z[, 2:ncol(output.tree$Z)]
#         output.tree.P = apply(
#             output.tree$P[2:nrow(output.tree$P),, drop = FALSE], 2,
#             function(x) { x/sum(x) })
#         output.tree$CCF = output.tree.Z %*% output.tree.P
#         samptreethin[[treei]] = output.tree
#     }
#     return(list(samptreethin, samptreethin.lik, config, config.summary))
# }
# 
# 
# canopy_output_to_df <- function(tree) {
#     p <- data_frame(clone = rownames(tree$P), prevalence = tree$P[, 1])
#     outdf <- as_data_frame(tree$Z) %>%
#         dplyr::mutate(var_id = rownames(tree$Z), vaf = tree$VAF[, 1],
#                       ccf = tree$CCF[, 1])
#     outdf[["clone"]] <- colnames(tree$Z)[1]
#     for (i in 2:ncol(tree$Z))
#         outdf[["clone"]][as.logical(tree$Z[, i])] <- colnames(tree$Z)[i]
#     outdf
# }
# 
# plot_tree <- function(tree, title = NULL, size = 10) {
#     ptree <- ggtree(tree, ladderize = FALSE, layout = "slanted") +
#         geom_tiplab(size = 10, color = "firebrick")
#     df_prev <- as.data.frame(tree$P)
#     df_prev[["clone"]] <- rownames(tree$P)
#     df_prev[["branch.y"]] <- ptree$data$branch.y[1:sum(!is.na(ptree$data$label))]
#     df_prev[["branch.y.f"]] <- as.factor(as.numeric(as.factor(df_prev[["branch.y"]])))
#     df_prev <- gather(df_prev, key = "sample", value = "prevalence", -clone,
#                       -branch.y, -branch.y.f)
#     df_prev <- dplyr::mutate(df_prev, prevalence = signif(prevalence, digits = 3))
#     ptile <- ggplot(df_prev, aes(x = sample, y = branch.y.f)) +
#         geom_tile(aes(fill = prevalence)) +
#         geom_label(aes(label = prevalence), colour = "firebrick", size = 5) +
#         scale_fill_viridis(option = "B") +
#         scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
#         theme_minimal() + theme(axis.title = element_blank(),
#                                 axis.text.y = element_blank())
#     pt <- plot_grid(ptree, ptile, nrow = 1, rel_widths = c(1, 0.5))
#     if (!is.null(title)) {
#         title <- ggdraw() + draw_label(title, fontface='bold')
#         pt <- plot_grid(title, pt, ncol = 1, rel_heights = c(0.1, 1))
#     }
#     pt
# }
