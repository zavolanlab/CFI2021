#!/bin/bash
# (c) Alexander Kanitz, Biozentrum, University of Basel
# alexander.kanitz@alumni.ethz.ch
# 12-DEC-2018

#####################
###  DESCRIPTION  ###
#####################

# Count overlaps of alignment (BAM) start positions with exonic/intronic
# features (GTF) per gene. Generate count table in TAB format with columns:
# - gene
# - width_exonic
# - width_intronic
# - count_exonic
# - count_intronic

####################
###  PARAMETERS  ###
####################

# Set parameters (DO NOT CHANGE!)
root="$(dirname $(cd "$(dirname "$0" )" && pwd))"
alignmentsDir="${root}/alignments"
resDir="${root}/results/counts/per_sample"
tmpDir="${root}/.tmp/counting"
logDir="${root}/docs"
script="${root}/scripts/gtf_to_counts_tab.R"
exonic_intronic_anno="${root}/annotations/hsa.GRCh38_90.gene_annotations.exonic_intronic..gtf.gz"

########################
###  PRE-REQUISITES  ###
########################

# Shell options
set -e
set -u
set -o pipefail

# Load modules
# NOTE: This is specific to the Biozentrum/sciCORE infrastructure. You need
# to make sure that R (ideally v3.5.0) and BEDTools (ideally v2.26.0) are in
# your path.
# module purge
# module load "R/3.5.0-goolf-1.7.20"
# module load "BEDTools/2.26.0-goolf-1.7.20"

# Create directories
mkdir --parents "$resDir"
mkdir --parents "$tmpDir"
mkdir --parents "$logDir"

# Create log file
logFile="${logDir}/$(basename $0 ".sh").log"
rm -f "$logFile"; touch "$logFile"
>&2 echo "Log written to '$logFile'..."

##############
###  MAIN  ###
##############

# Iterate over alignment files
for bam in "${alignmentsDir}/"*".bam"; do

    # Status message
    echo "Processing file '$bam'..." >> "$logFile"

    # Get filename prefix
    prefix=$(basename $bam ".bam")

    # Convert BAM to BED
    echo "Converting BAM file to BED12 format..." >> "$logFile"
    bed12="${tmpDir}/${prefix}.alignments.bed12.gz"
    bedtools bamtobed -bed12 -cigar -tag NH -i "$bam" | gzip > "$bed12" 2>> "$logFile"

    # Filter unique
    echo "Filtering unique mappers..." >> "$logFile"
    unique_mappers="${tmpDir}/${prefix}.unique_mappers.bed12.gz"
    awk '$5 == 1' <(zcat "$bed12") | gzip > "$unique_mappers" 2>> "$logFile"

    # Extract alignment start position
    echo "Extracting alignment start positions..." >> "$logFile"
    start_positions="${tmpDir}/${prefix}.start_positions.bed.gz"
    awk -v OFS='\t' '{$3 = $2 + 1; print $1,$2,$3,$4,$5,$6}' <(zcat "$unique_mappers") | gzip > "$start_positions" 2>> "$logFile"

    # Intersect with GTF file
    echo "Counting overlaps between alignment start positions and exonic/intronic regions..." >> "$logFile"
    overlaps="${tmpDir}/${prefix}.overlaps.gtf.gz"
    bedtools intersect -c -s -a "$exonic_intronic_anno" -b <(zcat "$start_positions") | awk -v FS='\t' -v OFS='\t' '{$9 = $9"; count \""$10"\";"; print}' | cut -f10 --complement | gzip > "$overlaps" 2>> "$logFile"

    # Convert to count table
    echo "Preparing count table..." >> "$logFile"
    counts="${resDir}/${prefix}.counts"
    "$script" "$overlaps" > "$counts" 2>> "$logFile"

done

#############
###  END  ###
#############

echo "Raw data in: $alignmentsDir" >> "$logFile"
echo "Temporary data in: $tmpDir" >> "$logFile"
echo "Result files in: $resDir" >> "$logFile"
echo "Done. No errors." >> "$logFile"
>&2 echo "Done. No errors."
