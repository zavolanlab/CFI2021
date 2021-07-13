#!/bin/bash
# (c) Alexander Kanitz, Biozentrum, University of Basel
# alexander.kanitz@alumni.ethz.ch
# 04-DEC-2018

#####################
###  DESCRIPTION  ###
#####################

# Filter gene annotations.

####################
###  PARAMETERS  ###
####################

# Set parameters (DO NOT CHANGE!)
root="$(dirname $(cd "$(dirname "$0" )" && pwd))"
geneAnno="../annotations/Homo_sapiens.GRCh38.90.chr.gtf"
fileNamePrefix="hsa.GRCh38_90"
resDir="${root}/annotations"
tmpDir="${root}/.tmp/annotations"
logDir="${root}/docs"

# Filters
# - All filters are positive filters, i.e. entries meeting the specified
# filters are kept
# - Separate multiple entries by a single space and quote the whole string
# - Set to empty string '""' if no filtering is desired
# - Transcriptome sequences are filtered according to the transcript
# annotations remaining after applying gene annotation filters
# - Warnings are issued if sequences for annotated transcripts are absent in
# the transcriptome
# Gene annotation / transcriptome filters
geneAnnoFilterChromosomes=""
geneAnnoFilterGeneBiotypes="lincRNA protein_coding"
geneAnnoFilterTranscriptBiotypes=""
geneAnnoFilterTranscriptSupportLevels="1"

########################
###  PRE-REQUISITES  ###
########################

# Shell options
set -e
set -u
set -o pipefail

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

# Get & filter gene annotations

# Compress gene annotations
geneAnnoGzip="${tmpDir}/${fileNamePrefix}.gene_annotations.gtf.gz"
echo "Compressing gene annotations..." >> "$logFile"
gzip --stdout "$geneAnno" > "$geneAnnoGzip"

# Filter gene annotations
geneAnnoFilt="${resDir}/${fileNamePrefix}.gene_annotations.filtered.gtf.gz"
geneAnnoFiltTmp="${tmpDir}/${fileNamePrefix}.gene_annotations.filtered.tmp.gtf.gz"
cp "$geneAnnoGzip" "$geneAnnoFiltTmp"

    # Filter requested chromosomes
    # ----------------------------
    # - If filter provided, filters comments and matching chromosomes
    if [ "$geneAnnoFilterChromosomes" != ""  ]; then
        echo "Filtering gene annotations by chromosomes..." >> "$logFile"
        perl -ane 'if(!@ARGV){if(/^#\!/){print}else{$keep=$chr{$F[0]}}}$keep?print:chomp;$chr{$_}=1 if @ARGV' <(echo "$geneAnnoFilterChromosomes" | sed 's/ /\n/g') <(zcat $geneAnnoFiltTmp) | gzip > "$geneAnnoFilt"
        cp "$geneAnnoFilt" "$geneAnnoFiltTmp"
    fi

    # Filter requested gene biotypes
    # ------------------------------
    # - If filter provided, filters comments and matching gene biotypes
    if [ "$geneAnnoFilterGeneBiotypes" != ""  ]; then
        echo "Filtering gene annotations by gene biotypes..." >> "$logFile"
        perl -ne 'if(/^#\!/){print;$keep=0}elsif(/gene_biotype\s\"(\S+)\"/){$keep=$type{$1}}else{$keep=0}$keep?print:chomp;$type{$_}=1 if @ARGV' <(echo "$geneAnnoFilterGeneBiotypes" | sed 's/ /\n/g') <(zcat $geneAnnoFiltTmp) | gzip > "$geneAnnoFilt"
        cp "$geneAnnoFilt" "$geneAnnoFiltTmp"
    fi

    # Filter requested transcript biotypes
    # ------------------------------------
    # - If filter provided, filters 'gene' entries, commentss and matching transcript biotypes
    if [ "$geneAnnoFilterTranscriptBiotypes" != ""  ]; then
        echo "Filtering annotations by transcript biotypes..." >> "$logFile"
        perl -ane 'if(/^#\!/||$F[2] eq "gene"){print;$keep=0}elsif(/transcript_biotype\s\"(\S+)\"/){$keep=$type{$1}}else{$keep=0}$keep?print:chomp;$type{$_}=1 if @ARGV' <(echo "$geneAnnoFilterTranscriptBiotypes" | sed 's/ /\n/g') <(zcat $geneAnnoFiltTmp) | gzip > "$geneAnnoFilt"
        cp "$geneAnnoFilt" "$geneAnnoFiltTmp"
    fi

    # Filter requested transcript support levels
    # ------------------------------------------
    # - If filter provided, filters 'gene' entries, comments and matching transcript support levels
    if [ "$geneAnnoFilterTranscriptSupportLevels" != ""  ]; then
        echo "Filtering annotations by transcript support levels..." >> "$logFile"
        perl -ane 'if(/^#\!/||$F[2] eq "gene"){print;$keep=0}elsif(/transcript_support_level\s\"(\S+?)\"?/){$keep=$level{$1}}else{$keep=0}$keep?print:chomp;$level{$_}=1 if @ARGV' <(echo "$geneAnnoFilterTranscriptSupportLevels" | sed 's/ /\n/g') <(zcat $geneAnnoFiltTmp) | gzip > "$geneAnnoFilt"
        cp "$geneAnnoFilt" "$geneAnnoFiltTmp"
    fi

    # Remove orphan 'genes' (i.e. 'genes' with all child entries removed) & temporary file
    echo "Removing 'orphan' genes..." >> "$logFile"
    perl -ane 'if ($F[2] eq "gene"){$prev=$_}else{print $prev,$_; $prev=""}' <(zcat $geneAnnoFiltTmp) | gzip > "$geneAnnoFilt"
    rm "$geneAnnoFiltTmp"

#############
###  END  ###
#############

echo "Original gene annotations in: $geneAnno" >> "$logFile"
echo "Filtered gene annotations in: $geneAnnoFilt" >> "$logFile"
echo "Done. No errors." >> "$logFile"
>&2 echo "Done. No errors."
