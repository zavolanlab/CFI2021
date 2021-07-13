#!/usr/bin/env Rscript

# Creation date: 2021-06-28
# Last modified: 2021-06-28
# Version: 1.0.0
# Author: Alexander Kanitz (alexander.kanitz@alumni.ethz.ch)

# Imports
print("Starting DGEA analysis.")
suppressMessages(library(edgeR))

# Get directory of this file
args <- commandArgs(trailingOnly = FALSE)
this_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))

# Create output directory
outpath <- file.path(this_dir, "..", "results", "dgea")
dir.create(outpath, showWarnings=FALSE)

# Load and transform data
# Needs to be adapted for each sample table
infile <- file.path(this_dir, "..", "results", "counts", "merged.all_experiments.counts")
dat <- read.delim(infile)
mt <- as.matrix(dat[, seq(4, ncol(dat), 2)])  # keep only exhjjjkjjionic
colnames(mt) <- gsub(".Aligned.out.sorted_ReadCount_Exonic", "", colnames(mt))  # simplify column names
rownames(mt) <- dat[[1]]  # set gene IDs as row names
mt[is.na(mt)] <- 0  # NAs result from zero width after collapsing intronic/exonic
mt_filt <- mt[, c(25, 27, 29, 1, 2, 6, 7, 8, 3, 4, 5, 9, 10, 11)]  # subset samples of interest
colnames(mt_filt)[1:3] <- c("WT_1", "WT_2", "WT_3")  # consistent sample naming scheme
outfile <- file.path(outpath, paste("dgea", "count_table", "tab", sep="."))
write.table(mt_filt, outfile, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

# Set sample groups
# Needs to be adapted for each sample table
grps <- factor(c(rep("1_wildtype", 3), rep("2_CF1m25_kd", 2), rep("2_CF1m68_kd", 3), rep("2_CF1m25_oe", 3), rep("2_CF1m68_oe", 3)))
cmps <- list(
  c("1_wildtype", "2_CF1m25_kd"),
  c("1_wildtype", "2_CF1m68_kd"),
  c("1_wildtype", "2_CF1m25_oe"),
  c("1_wildtype", "2_CF1m68_oe")
)

# edgeR analysis
dge_ls_all <- DGEList(mt_filt, group=grps)
dge_ls_filt <- dge_ls_all[filterByExpr(dge_ls_all), , keep.lib.sizes=FALSE]  # filter miRs by expression, across all groups
results <- lapply(cmps, function(cmp) {
    print(paste0("Comparing ", cmp[[2]], " to ", cmp[[1]], "..."))
    res <- list()
    res$dge_ls <- dge_ls_filt[, which(dge_ls_filt$samples$group %in% cmp)]
    res$dge_ls <- calcNormFactors(res$dge_ls)
    suppressMessages(res$dge_ls <- estimateDisp(res$dge_ls))
    res$et <- exactTest(res$dge_ls)
    res$diff_tab <- topTags(res$et, n=nrow(res$et$table))
    outfile <- file.path(outpath, paste("dgea", cmp[[1]], "vs", cmp[[2]], "tab", sep="."))
    write.table(res$diff_tab$table, outfile, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
    return(res)
})

# Save image
print("Storing session info and analysis image...")
outfile <- file.path(outpath, paste("dgea", "session_info", "txt", sep="."))
writeLines(capture.output(sessionInfo()), outfile)
outfile <- file.path(outpath, paste("dgea", "Rdata", sep="."))
save.image(file=outfile)
print("Done.")
