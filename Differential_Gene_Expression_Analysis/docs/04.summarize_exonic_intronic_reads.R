#!/usr/bin/env Rscript
# (c) Alexander Kanitz, Biozentrum, University of Basel
# alexander.kanitz@alumni.ethz.ch
# 13-DEC-2018

# Merge count tables according for downstream application

# Global paramters (DO NOT CHANGE!)
args <- commandArgs(trailingOnly = FALSE)
this_dir <- dirname(sub("--file=", "", args[grep("--file=", args)]))
res_dir <- file.path(this_dir, "..", "results", "counts")
input_dir <- file.path(res_dir, "per_sample")
output_file <- file.path(res_dir, "merged.all_experiments.counts")
pattern <- ".counts"
gene_column <- 1
width_columns <- c(2, 3)
count_columns <- c(4, 5)
replacement_suffix <- "_ReadCount"
suffixes <- c("Exonic", "Intronic")
sep="_"

# Merging function
merge_df_from_ls <- function(ls, by=0, all=TRUE) {
    mdf <- ls[[1]]
    for (df in 2:length(ls)) {
        mdf <- merge(mdf, ls[[df]], by=by, all=all)
        if (by == 0) {
            rownames(mdf) <- mdf$Row.names
            mdf <- mdf[,-1]
        }
    }
    return(mdf)
}

# Find, load & rename files
files <- dir(input_dir, pattern=pattern, full.names=TRUE)
dfl <- lapply(files, read.delim, header=TRUE, stringsAsFactors=FALSE)
names(dfl) <- sub(pattern, replacement_suffix , basename(files))

# Extract features widths & counts
widths <- dfl[[1]][, c(gene_column, width_columns)]
counts <- lapply(names(dfl), function(nm) {
    tmp <- dfl[[nm]][, c(gene_column, count_columns)]
    colnames(tmp)[2:3] <- paste(nm, suffixes, sep=sep)
    return(tmp)
})
ls <- c(list(widths), counts)

# Merge data
df.merged <- merge_df_from_ls(ls, by=1, all=TRUE)

# Write output
write.table(
    df.merged,
    output_file,
    col.names=TRUE,
    row.names=FALSE,
    quote=FALSE,
    sep="\t"
)
