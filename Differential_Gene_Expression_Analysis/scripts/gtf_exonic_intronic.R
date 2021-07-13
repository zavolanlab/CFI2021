#!/usr/bin/env Rscript
# (c) Alexander Kanitz, Biozentrum, University of Basel
# alexander.kanitz@alumni.ethz.ch
# 04-DEC-2018
# Dependencies:
# - R v3.5.0
# - rtracklayer_1.40.0
# - GenomicFeatures_1.32.0


#=============#
#   OPTIONS   #
#=============#

# Load option parser
if ( suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE ) { stop("Package 'optparse' required!\nExecution aborted.") }

# Get script name
script <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

# Compile help message
description <- 
"From a GTF gene annotation file, generates a GTF file containing exclusively
'exonic' and 'intronic' regions generated from intersecting the exons/introns
of all transcripts of a given gene.

Genic regions that are annotated as exonic in some and as intronic in other
transcripts are discarded.

First (most 5') and/or last (most 3') exons of each transcript can optionally
be discarded before compiling exonic and intronic regions.

The input file can be either uncompressed or in GZIP format.
"
author <- "Author: Alexander Kanitz, Biozentrum, University of Basel"
created <- "Created: 04-DEC-2018"
version <- "Version: 1.1.1 (20-SEP-2019)"
requirements <- "Requires: optparse, rtracklayer"
msg <- paste(description, author, created, version, requirements, sep="\n")

# List arguments
option_list <- list(
        make_option(c("-i", "--gtf"), action="store", type="character", default="", help="Path to GTF input file (required).", metavar="file"),
        make_option("--remove-first-exon", action="store_true", default=FALSE, help="Remove first (most 5') exon of each transcript."),
        make_option("--remove-last-exon", action="store_true", default=FALSE, help="Remove last (most 3') exon of each transcript."),
        make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die!"),
        make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die!"),
        make_option(c("-v", "--verbose"), action="store_true", default=FALSE, help="Print log messages to STDOUT.")
)

# Parse arguments
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --gtf PATH\n", option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

# Die if required arguments are missing...
if      ( opt$`gtf`  == "" ) { 
    write("[ERROR] Required argument(s) missing!\n\n", stderr())
    stop(print_help(opt_parser))
}


#==========#
#   MAIN   #
#==========#

if ( opt$verbose ) cat("Starting '", script, "'...\n", sep="", file = stderr())

### DEBUGGING ###
#opt <- list()
#opt$gtf <- "annotations.test.gtf"
#opt$verbose <- TRUE
#opt$`remove-last-exon` <- FALSE
#opt$`remove-first-exon` <- TRUE
#################

# Load packages
if ( opt$verbose ) cat("Loading required packages...\n", file = stderr())
if ( suppressWarnings(suppressPackageStartupMessages(require("rtracklayer"))) == FALSE ) { stop("Package 'rtracklayer' required!\nExecution aborted.") }

# Import annotations
if ( opt$verbose ) cat("Importing annotation from file '", opt$`gtf`, "'...\n", sep="", file = stderr())
gr.raw <- gr <- import(opt$`gtf`, format="gtf")

# Keep only exon, transcript and gene entries
gr.filt <- gr <- gr[values(gr)[["type"]] %in% c("exon", "transcript")]

# Extract transcripts
trx <- gr[values(gr)[["type"]] == "transcript"]

# Remove first and/or last exon of each transcript
if ( opt[["remove-first-exon"]] || opt[["remove-last-exon"]] ) {
    if ( opt$verbose ) cat("Removing first and/or last exons of each transcript...\n", file = stderr())
    ex.raw <- ex <- gr[values(gr)[["type"]] == "exon"]
    ex_by_trx <- split(ex, f=values(ex)[["transcript_id"]])
    ex_by_trx_filt <- endoapply(ex_by_trx, function(trx) {
        strand <- unique(as.character(strand(trx)))
        if (length(strand) != 1) {
            cat("Transcripts have no or multiple strands annotated. Check your input file. Execution aborted!\n", file = stderr())
        }
        if ( ( strand == "+" && opt[["remove-first-exon"]] ) || ( strand == "-" && opt[["remove-last-exon"]] ) ) {
            trx <- trx[-1]
        }
        if ( ( strand == "+" && opt[["remove-last-exon"]] ) || ( strand == "-" && opt[["remove-first-exon"]] ) ) {
            trx <- trx[-length(trx)]
        }
        return(trx)
    })
    ex_filt <- unlist(ex_by_trx_filt)
    trx_filt <- trx[values(trx)[["transcript_id"]] %in% unique(values(ex_filt)[["transcript_id"]])]
    gr_filt <- sort(c(trx_filt, ex_filt))
    names(gr_filt) <- NULL
    gr <- gr_filt
}

# Extract exonic regions (all regions/fragments of a gene that are annotated as an exon in at least one transcript)
if ( opt$verbose ) cat("Extracting exonic regions...\n", file = stderr())
ex <- gr[values(gr)[["type"]] == "exon"]
ex_by_gene <- split(ex, f=values(ex)[["gene_id"]])
exonic <- unlist(reduce(ex_by_gene)) # GRanges

# Extract gene boundaries from remaining exons
# Boundaries are not taken from the gene entries because previous filtering may not have considered resetting boundaries
if ( opt$verbose ) cat("Finding gene boundaries...\n", file = stderr())
gene_boundaries <- range(ex_by_gene) # GRangesList

# Extract intronic regions (all regions/fragments within a gene boundary that are annotated as an intron in at least one transcript)
# Obtained by determining set difference between gene boundaries and exons of each transcript
if ( opt$verbose ) cat("Extracting intronic regions...\n", file = stderr())
ex_by_gene_and_trx <- lapply(ex_by_gene, function(gr) split(gr, f=values(gr)[["transcript_id"]])) # List of GRangesList (gene -> trx -> ex); takes long
int <- mapply(function(boundaries, exons_by_trx) endoapply(exons_by_trx, function(ex) setdiff(boundaries, ex)), gene_boundaries, ex_by_gene_and_trx) # takes long
intronic <- unlist(reduce(GRangesList(lapply(int, unlist)))) # takes long

# Merging exclusive exonic/intronic regions
if ( opt$verbose ) cat("Merging exclusive exonic/intronic regions...\n", file = stderr())
exonic.only <- setdiff(exonic, intronic)
intronic.only <- setdiff(intronic, exonic)
mcols(exonic.only)[["type"]] <- "exonic"
mcols(intronic.only)[["type"]] <- "intronic"
features <- sort(c(exonic.only, intronic.only))
if ( length(features) == 0 ) stop("No exonic/intronic regions found! Check file.\nExecution halted.\n")

# Annotate exonic/intronic regions with gene IDs
if ( opt$verbose ) cat("Annotate exonic/intronic regions with gene identifiers...\n", file = stderr())
overlaps <- findOverlaps(features, trx, type="within")
features.out.not_unique <- features.out <- features[from(overlaps)]
mcols(features.out)[["gene_id"]] <- mcols(trx)[["gene_id"]][to(overlaps)]
features.out <- unique(features.out)

# Write out regions
if ( opt$verbose ) cat("Writing regions to GTF file...\n", sep="", file = stderr())
export(features.out, "", format="gtf")

### DEBUGGING ###
# Save session
session_path <- paste(script, "Rimage", sep=".")
if ( opt$verbose ) cat("Writing session to '", session_path, "'...\n", sep="", file = stderr())
save.image(session_path)
#################

if ( opt$verbose ) cat("Done.\n", file = stderr())

