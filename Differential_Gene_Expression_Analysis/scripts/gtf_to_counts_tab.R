#!/usr/bin/env Rscript
# (c) Alexander Kanitz, Biozentrum, University of Basel
# alexander.kanitz@alumni.ethz.ch
# 12-DEC-2018
# Dependencies:
# - R v3.5.0
# - rtracklayer_1.40.0


# Converts a GTF file of "intronic" and "exonic" features ("type" column) with column "count" to a 
# TAB count table of format "gene"/"width_exonic"/"width_intronic"/"count_exonic"/"count_intronic"


#==========#
#   MAIN   #
#==========#

# Import packages
if ( suppressWarnings(suppressPackageStartupMessages(require("rtracklayer"))) == FALSE ) { stop("Package 'rtracklayer' required!\nExecution aborted.") }

# Process CLI arguments
args = commandArgs(trailingOnly=TRUE)
gtf <- args[[1]]

# Import GTF file
print("blabla")
gr <- import(gtf, format="gtf")
print("blabla2")
print(head(gr))
# Cast counts to integer
#values(gr)[["count"]] <- as.integer(values(gr)[["count"]])
values(gr)[["score"]] <- as.integer(values(gr)[["score"]])

# Convert GRanges to data frame
df <- as.data.frame(gr)

# Separate exonic and intronic features
df.ex <- df[df[["type"]] == "exonic", ]
df.int <- df[df[["type"]] == "intronic", ]

# Split by gene ID
ls.ex = split(df.ex, f=df.ex$gene_id)
ls.int = split(df.int, f=df.int$gene_id)

# Get data frames of widths and counts
widths_counts.ex <- t(sapply(ls.ex, function(df) setNames(c(sum(df[["width"]]), sum(df[["count"]])), c("width_exonic", "counts_exonic"))))
widths_counts.int <- t(sapply(ls.int, function(df) setNames(c(sum(df[["width"]]), sum(df[["count"]])), c("width_intronic", "counts_intronic"))))

# Merge exonic and intronic features
df.merge <- merge(widths_counts.ex, widths_counts.int, by=0, all=TRUE)

# Rename and reorder columns
colnames(df.merge)[[1]] <- "gene"
df.out <- df.merge[, c(1, 2, 4, 3, 5)]

# Write out results
write.table(df.out, "", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
