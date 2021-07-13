#!/usr/bin/env Rscript
# (c) 2016, Alexander Kanitz, Biozentrum, Universiry of Basel
# email: alexander.kanitz@alumni.ethz.ch

#######################
###  PARSE OPTIONS  ###
#######################

#---> LOAD OPTION PARSER <---#
if ( suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE ) { stop("[ERROR] Package 'optparse' required! Aborted.") }

#---> GET SCRIPT NAME <---#
script <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

#---> DESCRIPTION <---#
description <- "Merges a set of data tables based on a common ID column.\n"
author <- "Author: Alexander Kanitz, Biozentrum, University of Basel"
version <- "Version: 1.0.0 (07-DEC-2016)"
requirements <- "Requires: optparse"
msg <- paste(description, author, version, requirements, sep="\n")

#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
                make_option(
                    c("-i", "--input-directory"),
                    action="store",
                    type="character",
                    default=NULL,
                    help="Directory containing data tables. Required.",
                    metavar="directory"
                ),
                make_option(
                    "--output-table",
                    action="store",
                    type="character",
                    default=".",
                    help="Filename of the merged table. Required.",
                    metavar="file"
                ),
                make_option(
                    "--glob",
                    action="store",
                    type="character",
                    default="*",
                    help="Glob for selecting files in input directory. Default: \"*\".",
                    metavar="glob"
                ),
                make_option(
                    "--id-column",
                    action="store",
                    type="integer",
                    default=1,
                    help="The field/column containing the row/feature names for the merged table(s). Default: 1.",
                    metavar="int"
                ),
                make_option(
                    "--no-header",
                    action="store_true",
                    default=FALSE,
                    help="Specify if the data tables do *not* have a header line."
                ),
                make_option(
                    "--all-rows",
                    action="store_true",
                    default=FALSE,
                    help="Indicate whether all rows of all tables shall be included in the output. In that case, missing data points are filled up with NAs. By default, only rows common to all data tables are included."
                ),
                make_option(
                    c("-h", "--help"),
                    action="store_true",
                    default=FALSE,
                    help="Show this information and die."
                ),
                make_option(
                    c("-u", "--usage"),
                    action="store_true",
                    default=FALSE,
                    dest="help",
                    help="Show this information and die."
                ),
                make_option(
                    c("-v", "--verbose"),
                    action="store_true",
                    default=FALSE,
                    help="Print log messages to STDOUT."
                )
)

## Parse command-line arguments
opt_parser <- OptionParser(usage=paste("Usage:", script, "[OPTIONS] --input-directory <DIRECTORY>\n", sep=" "), option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

#---> ARGUMENT VALIDATION <---#

## Die if any required arguments are missing...
if ( is.null(opt$`input-directory`) || is.null(opt$`output-table`) ) {
        print_help(opt_parser)
        stop("[ERROR] Required argument missing! Aborted.")
}


###################
###  FUNCTIONS  ###
###################

# Merge multiple dataframes
merge_multi <- function(ls, by=0, all=FALSE) {

    # If the list is of zero length, return NULL
    if (! length(ls) ) return(NULL)

    # Initialize merged dataframe with first dataframe
    df.merged <- ls[[1]]
    for ( idx in seq_len(length(ls))[-1] ) {
        tmp <- merge(df.merged, ls[[idx]], by=by, all=all)
        rownames(tmp) <- tmp[, 1]
        df.merged <- tmp[, -1]
    }

    # Return merged dataframe
    return(df.merged)

}


##############
###  MAIN  ###
##############

#---> START MESSAGE <---#
if ( opt$`verbose` ) cat("Starting ", script, "...\n", sep="'")

#---> IMPORT DATA TABLES <---#

    #---> Non-recursively find files of specified name <---#
    if ( opt$`verbose` ) cat("Finding data tables...\n")
    file.paths <- sort(dir(opt$`input-directory`, pattern=glob2rx(opt$`glob`), recursive=FALSE, full.names=TRUE))

    #---> Read data tables <---#
    if ( opt$`verbose` ) cat("Loading data tables (may take long)...\n")
    df.ls <- lapply(file.paths, read.delim, header=!opt$`no-header`, row.names=opt$`id-column`, stringsAsFactors=FALSE)


#---> MERGING DATA TABLES <---#

    #---> Log message <---#
    if ( opt$`verbose` ) cat("Merging data tables (may take long)...\n", sep="'")

    #---> Merge data tables <---#
    merged <- merge_multi(df.ls, by=0, all=opt$`all-rows`)

#---> WRITING OUTPUT <---#

    #---> Log message <---#
    if ( opt$`verbose` ) cat("Writing merged data table to file ", opt$`output-table`, "...\n", sep="'")

    #---> Write tables <---#
    write.table(merged, opt$`output-table`, col.names=!opt$`no-header`, row.names=TRUE, sep="\t", quote=FALSE)

#---> END MESSAGE <---#
if ( opt$`verbose` ) cat("Done.\n")
