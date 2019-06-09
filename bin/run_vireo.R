### Rscript for donor deconvolution from command line
# This file allows to directly run from command line via:
# Rscript run_vireo.R $CELL_FILE $OUT_DIR $N_DONOR $GT_VCF_FILE

# Note, $GT_VCF_FILE is optional. Both $CELL_FILE and $GT_VCF_FILE can be either
# vcf.gz or rds.
# It will save results in to the $OUT_DIR folder with three files: 
# 1) summary.tsv: number of cells for each donor
# 2) donor_ids.tsv: the donor id for each cell
# 3) res_ids.rds: the full data set in as an R data set. This file can be big, 
#    and can be removed if not needed.


args <- commandArgs(TRUE)
if (length(args) < 3) {
    stop("require at least 3 parameters")
}
cell_file <- args[1]
out_dir <- args[2]
n_donor <- as.numeric(args[3])
if (length(args) >= 4) {
    GT_file <- args[4]
} else {
    GT_file <- NULL
}

library(cardelino)

### check out_dir
if (!file.exists(out_dir)) {
    dir.create(out_dir)
}
setwd(out_dir)
set.seed(1)

### load data
if (tools::file_ext(cell_file) == "rds") {
    cell_data <- readRDS(cell_file)
} else if (tools::file_ext(cell_file) == "h5" || 
           tools::file_ext(cell_file) == "hdf5") {
    cell_data <- load_vcf_h5(cell_file)
} else {
    cell_data <- load_cellSNP_vcf(cell_file)
}

### load GT data and match variants
if (is.null(GT_file)) {
    GT_data <- NULL
} else {
    if (tools::file_ext(GT_file) == "rds") {
        GT_data <- readRDS(GT_file)
    } else {
        GT_data <- load_GT_vcf(GT_file)
    }

    mm <- match(rownames(GT_data$GT), rownames(cell_data$A))
    print(paste(sum(!is.na(mm)), "out of", length(mm), 
                "variants mathed between GT and cells."))

    idx <- which(!is.na(mm))
    GT_use <- GT_data$GT[idx, ]
    
    ## check if GP exists
    if (!is.null(GT_data$GP)) {
        idx_GP <- c()
        for (ii in seq_len(ncol(GT_data$GT))) {
            idx_GP <- c(idx_GP, idx + nrow(GT_data$GT) * (ii - 1))
        }
        GP_use <- GT_data$GP[idx_GP, ]
    }
    
    cell_data$A <- cell_data$A[mm[!is.na(mm)], ]
    cell_data$D <- cell_data$D[mm[!is.na(mm)], ]
}

### run Vireo
if (is.null(GT_data)) {
    ids_res <- vireo(cell_data = cell_data, n_donor = n_donor, 
                     n_proc = 4, n_init = 10)
} else {
    if (!is.null(GT_data$GP)) {
        ids_res <- vireo(cell_data = cell_data, n_donor = n_donor, K_amplify = 1, 
                         GT_prior = GP_use, n_proc = 2, n_init = 2)
        for (ii in seq_len(ncol(GT_use))) {
            idx_don <- ids_res$assigned$donor_id == paste0("donor", ii)
            ids_res$assigned$donor_id[idx_don] <- colnames(GT_use)[ii]
        }
        
    } else {
        ids_res <- vireo(cell_data = cell_data, donor_data = GT_data,
                         n_proc = 2, n_init = 2)
    }
}
table(ids_res$assigned$donor_id)

### save results
write.table(table(ids_res$assigned$donor_id), 
            file = paste0(out_dir, "/summary.tsv"), 
            quote = FALSE, sep = "\t", row.names = FALSE)

write.table(format(ids_res$assigned, digits = 3), 
            file = paste0(out_dir, "/donor_ids.tsv"), 
            quote = FALSE, sep = "\t", row.names = FALSE)

### this file may be very big
ids_res$GT_input <- GT_data
saveRDS(ids_res, file = paste0(out_dir, "/res_ids.rds"))
