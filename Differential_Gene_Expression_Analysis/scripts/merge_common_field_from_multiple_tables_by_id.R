#!/usr/bin/env Rscript
# (c) 2016-2017, Alexander Kanitz, Biozentrum, University of Basel
# email: alexander.kanitz@alumni.ethz.ch

#######################
###  PARSE OPTIONS  ###
#######################

#---> LOAD OPTION PARSER <---#
if ( suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE ) { stop("[ERROR] Package 'optparse' required! Aborted.") }

#---> GET SCRIPT NAME <---#
script <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

#---> DESCRIPTION <---#
description <- "Merges a field of a set of data tables based on a common ID column.\n"
author <- "Author: Alexander Kanitz, Biozentrum, University of Basel"
version <- "Version: 1.1.0 (10-JAN-2017)"
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
                    "--recursive",
                    action="store_true",
                    default=FALSE,
                    help="Specify if subdirectories shall be searched for data tables."
                ),
                make_option(
                    "--output-directory",
                    action="store",
                    type="character",
                    default=".",
                    help="Directory where merged tables should be written. Default: \".\".",
                    metavar="directory"
                ),
                make_option(
                    "--out-file-suffix",
                    action="store",
                    type="character",
                    default="",
                    help="Output filenames are constructed from category names and this suffix. Default: \"\".",
                    metavar="string"
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
                    "--id-prefix",
                    action="store",
                    type="character",
                    default="",
                    help="Common prefix of data table filenames that is not part of the data table IDs. Default: \"\".",
                    metavar="string"
                ),
                make_option(
                    "--id-suffix",
                    action="store",
                    type="character",
                    default="",
                    help="Common suffix of data table filenames that is not part of the data table IDs. Default: \"\".",
                    metavar="string"
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
                    "--data-column",
                    action="store",
                    type="integer",
                    default=2,
                    help="The field/column containing the data to be merged. Default: 2",
                    metavar="int"
                ),
                make_option(
                    "--has-header",
                    action="store_true",
                    default=FALSE,
                    help="Indicates whether the data tables have a header line. Default: FALSE."
                ),
                make_option(
                    "--annotation-table",
                    action="store",
                    type="character",
                    default=NULL,
                    help="A table containing information for creating data table IDs and categories. Merged output tables will be generated for each unique category. Default: NULL.",
                    metavar="file"
                ),
                make_option(
                    "--anno-id-columns",
                    action="store",
                    type="character",
                    default="1",
                    help="One or more annotation table fields of which IDs matching the data tables can be constructed. Separate multiple fields/columns by comma, with no whitespace. Default: \"1\".",
                    metavar="string"
                ),
                make_option(
                    "--category-columns",
                    action="store",
                    type="character",
                    default="2",
                    help="One or more annotation table fields indicating data categories. Separate multiple fields/columns by comma, with no whitespace. Default: \"2\".",
                    metavar="string"
                ),
                make_option(
                    "--anno-id-separator",
                    action="store",
                    type="character",
                    default="_",
                    help="Character separating multiple annotation table fields when constructing data table IDs. Default: \"_\".",
                    metavar="char"
                ),
                make_option(
                    "--category-separator",
                    action="store",
                    type="character",
                    default=".",
                    help="Character for separating joined categories. Default: \".\"",
                    metavar="char"
                ),
                make_option(
                    "--anno-has-header",
                    action="store_true",
                    default=FALSE,
                    help="Indicates whether the annotation table has a header line. Default: FALSE."
                ),
                make_option(
                    "--whitespace-replacement",
                    action="store",
                    type="character",
                    default="_",
                    help="Replacement character for whitespaces in category names. Default: \"_\".",
                    metavar="char"
                ),
                make_option(
                    "--include-ungrouped",
                    action="store_true",
                    default=FALSE,
                    help="Generate merged table including all data (i.e. not filtered by categories). Will always be generated if annotation table is not provided."
                ),
                make_option(
                    "--prefix-ungrouped",
                    action="store",
                    type="character",
                    default="all",
                    help="Prefix for output table including all data. Default: \"all\"",
                    metavar="char"
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
if ( is.null(opt$`input-directory`) ) {
        print_help(opt_parser)
        stop("[ERROR] Required argument missing! Aborted.")
}

## Set option dependencies
if ( is.null(opt$`annotation-table`) ) {
    opt$`include-ungrouped` <- TRUE
}

##############
###  MAIN  ###
##############

#---> START MESSAGE <---#
if ( opt$`verbose` ) cat("Starting ", script, "...\n", sep="'")

#---> FIND DATA TABLES, EXTRACT TABLE IDS AND IMPORT DATA <---#

    #---> Non-recursively find files of specified name <---#
    if ( opt$`verbose` ) cat("Finding data tables...\n")
    file.paths <- sort(dir(opt$`input-directory`, pattern=glob2rx(opt$`glob`), recursive=opt$`recursive`, full.names=TRUE))

    #---> Extract sample IDs by removing specified pre- and suffix from file basenames <---#
    if ( opt$`verbose` ) cat("Extracting table IDs...\n")
    sample.ids <- sort(gsub(opt$`id-suffix`, "", gsub(opt$`id-prefix`, "", basename(file.paths))))

    #---> Read data tables <---#
    if ( opt$`verbose` ) cat("Loading data tables (may take long)...\n")
    df.ls <- lapply(file.paths, read.table, header=opt$`has-header`, sep="\t", row.names=opt$`id-column`, stringsAsFactors=FALSE)
    vec.ls <- lapply(df.ls, function(df) {
        setNames(df[, opt$`data-column`-1], rownames(df))
    })
    names(vec.ls) <- sample.ids

#---> COMPILE LIST OF CATEGORIES AND CORRESPONDING TABLE IDS <---#

    #---> Log message <---#
    if ( opt$`verbose` ) cat("Categorizing data tables...\n")

    #---> Initialize empty category list... <---#
    if ( is.null(opt$`annotation-table`) ) {

        # Initialize container list
        cat.ls <- list()

    #---> ...or get categories from provided table <---#
    } else {
        #---> Parse ID and category columns <---#
        cols.id <- as.integer(unlist(strsplit(opt$`anno-id-columns`, ",")))
        cols.cats <- as.integer(unlist(strsplit(opt$`category-columns`, ",")))

        #---> Import annotation table <---#
        anno.table <- read.delim(opt$`annotation-table`, header=opt$`anno-has-header`)

        #---> Get and verify IDs <---#
        ids <- as.character(apply(anno.table[, cols.id, drop=FALSE], 1, paste, collapse=opt$`anno-id-separator`))
        if ( ! all(sample.ids %in% ids) ) {
            stop("[ERROR] Annotations not available for all data tables. Aborted.")
        } else {
            ids <- ids[ids %in% sample.ids]
        }

        #---> Order annotation table by IDs and sort IDs <---#
        anno.table <- anno.table[order(ids), ]
        ids <- sort(ids)

        #---> Get categories <---#
        cats <- as.factor(gsub("\\s+", opt$`whitespace-replacement`, apply(anno.table[, cols.cats, drop=FALSE], 1, paste, collapse=opt$`category-separator`)))

        #---> Compile category list <---#
        cat.ls <- split(ids, cats)

    }

    #---> If specifically requested or no annotations provided...
    if ( opt$`include-ungrouped` ) {

        #---> ...add all sample IDs <---#
        cat.ls$all <- sample.ids

    }

#---> MERGING AND WRITING OUT DATA BY CATEGORY <---#

    #---> Iterate over categories in category list... <---#
    for ( cat in names(cat.ls) ) {

        #---> Log message <---#
        if ( opt$`verbose` ) cat("Merging data for category ", cat, " (may take long)...\n", sep="'")

        #---> Filter data <---#
        vec.ls.filt <- vec.ls[cat.ls[[cat]]]

        #---> Create dummy dataframe with all IDs <---#
        ids <- sort(unique(unlist(lapply(vec.ls.filt, names))))
        merged <- data.frame(row.names=ids)

        #---> Iterate over data vectors  <---#
        for (vec in vec.ls.filt) {

            # Sort data by ID
            vec <- vec[order(names(vec))]

            # Get data for new column
            add <- ifelse(rownames(merged) %in% names(vec), vec, NA)

            # Add new column
            merged <- cbind(merged, add)

        }

        #---> Set column/sample names <---#
        colnames(merged) <- cat.ls[[cat]]

        #---> Build output filename <---#
        if ( cat == "all" ) cat <- opt$`prefix-ungrouped`
        out.filename <- file.path(opt$`output-directory`, paste(cat, opt$`out-file-suffix`, sep=""))

        #---> Log message <---#
        if ( opt$`verbose` ) cat("Writing merged data to file ", out.filename, "...\n", sep="'")

        #---> Write tables <---#
        write.table(merged, out.filename, col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)

    }

#---> END MESSAGE <---#
if ( opt$`verbose` ) cat("Done.\n")
