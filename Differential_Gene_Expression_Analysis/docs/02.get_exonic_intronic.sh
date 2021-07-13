#!/bin/bash
# (c) Alexander Kanitz, Biozentrum, University of Basel
# alexander.kanitz@alumni.ethz.ch
# 04-DEC-2018

#####################
###  DESCRIPTION  ###
#####################

# Prepare GTF file of pure exon and intron parts of genes.

####################
###  PARAMETERS  ###
####################

# Set parameters (DO NOT CHANGE!)
root="$(dirname $(cd "$(dirname "$0" )" && pwd))"
resDir="${root}/annotations"
logDir="${root}/docs"
script="${root}/scripts/gtf_exonic_intronic.R"
fileNamePrefix="hsa.GRCh38_90"
geneAnno="${root}/annotations/${fileNamePrefix}.gene_annotations.filtered.gtf.gz"
outFile="${resDir}/${fileNamePrefix}.gene_annotations.exonic_intronic.gtf.gz"

########################
###  PRE-REQUISITES  ###
########################

# Shell options
set -e
set -u
set -o pipefail

# Load modules
# NOTE: This is specific to the Biozentrum/sciCORE infrastructure. You need
# to make sure that R (ideally of the same version, v3.5.0) is in your path.
# module purge
# module load "R/3.5.0-goolf-1.7.20"

# Create directories
mkdir --parents "$resDir"
mkdir --parents "$logDir"

# Create log file
logFile="${logDir}/$(basename $0 ".sh").log"
rm -f "$logFile"; touch "$logFile"
>&2 echo "Log written to '$logFile'..."

##############
###  MAIN  ###
##############

# Extract exonic and intronic regions for each gene
echo "Extracting exonic and intronic regions for each gene..." >> "$logFile"
"$script" --gtf "$geneAnno" --remove-last-exon --verbose 2>> "$logFile" | gzip > "$outFile"

#############
###  END  ###
#############

echo "Original gene annotations in: $geneAnno" >> "$logFile"
echo "Filtered gene annotations in: $outFile" >> "$logFile"
echo "Done. No errors." >> "$logFile"
>&2 echo "Done. No errors."
