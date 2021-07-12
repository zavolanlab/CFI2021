#
# A short piece of code to create the final plots for the paper
# with the coverage profiles.
#

#---> LOAD LIBRARIES <---#
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("Gviz"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("GenomicFeatures"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("GenomicAlignments"))

# mimic the CLargs:
opt_parser <- OptionParser()
opt <- parse_args(opt_parser)
opt$verbose = TRUE
opt$gtf = "Homo_sapiens.GRCh38.90.chr.gtf" # assume that genomic annotation is copied alongside the script into the same dir!
opt$bed = "regions.bed"
opt$design_table = "design_table.tsv"
opt$output_dir = "plots"

#---> Configuration for non ucsc chromosomes (chomosome that start with chr) <---
options(ucscChromosomeNames=FALSE)

#---> IMPORT GTF <---#
# Use rtracklayer::import method to import GTF file to GRanges object
txdb <- makeTxDbFromGFF(opt$gtf)
geneTrack <- GeneRegionTrack(txdb,fontcolor="black",collapseTranscripts="meta",transcriptAnnotation="gene",col="#a9a9a9",fill="#a9a9a9")
transcriptTrack <- GeneRegionTrack(txdb,fontcolor="black",collapseTranscripts=FALSE,transcriptAnnotation="transcript",col="#a9a9a9",fill="#a9a9a9")
names(geneTrack) = ""
names(transcriptTrack) = ""

#---> IMPORT bed regions <---#
bed <- read.table(opt$bed,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

#---> Create output directory <---#
dir.create(opt$output_dir, showWarnings = FALSE)

#---> Get the gene information <---#
# Use rtracklayer::import method to import GTF file to GRanges object 
gr <- import(con=opt$gtf, format="gtf")
gr <- gr[values(gr)[["type"]] == "gene"]
gr <- gr[gr$gene_id %in% bed[["V9"]] ]

#---> Read the design table <---#
design_table <- read.table(opt$design_table,header = TRUE, sep="\t",stringsAsFactors=FALSE, quote="",comment.char="")

# Get the coverage counts
coverages = apply(design_table, 1, function(x){
	return(coverage(readGAlignments(x["bam"]))[gr])
	})

### TIMP2
svg(file.path(opt$output_dir,"TIMP.svg"), width=20, height=12)
alignment_tracks = apply(design_table, 1, function(x){
  track = AlignmentsTrack(x["bam"], isPaired=FALSE, fill.coverage=x["color"], lwd.coverage="2", col.coverage=x["color"], alpha=1) 
  names(track) = x["sample"]
  return(track)
})
plotTracks(c(alignment_tracks,transcriptTrack),
           chromosome=bed[1,"V1"],
           from=bed[1,"V2"],
           to=bed[1,"V3"],
           type="coverage",
           main=bed[1,"V4"],
           ylim=c(0,250),
           cex.main=1.0,
           sizes=c(rep(2,length(alignment_tracks)),10)
)
dev.off()

### AGO2
svg(file.path(opt$output_dir,"AGO2.svg"), width=20, height=12)
alignment_tracks = apply(design_table, 1, function(x){
  track = AlignmentsTrack(x["bam"], isPaired=FALSE, fill.coverage=x["color"], lwd.coverage="2", col.coverage=x["color"], alpha=1) 
  names(track) = x["sample"]
  return(track)
})
plotTracks(c(alignment_tracks,transcriptTrack),
           chromosome=bed[2,"V1"],
           from=bed[2,"V2"],
           to=bed[2,"V3"],
           type="coverage",
           main=bed[2,"V4"],
           ylim=c(0,320),
           cex.main=1.0,
           sizes=c(rep(2,length(alignment_tracks)),10)
)
dev.off()

### DICER1
svg(file.path(opt$output_dir,"DICER1.svg"), width=20, height=12)
alignment_tracks = apply(design_table, 1, function(x){
  track = AlignmentsTrack(x["bam"], isPaired=FALSE, fill.coverage=x["color"], lwd.coverage="2", col.coverage=x["color"], alpha=1) 
  names(track) = x["sample"]
  return(track)
})
plotTracks(c(alignment_tracks,transcriptTrack),
           chromosome=bed[3,"V1"],
           from=bed[3,"V2"],
           to=bed[3,"V3"],
           type="coverage",
           main=bed[3,"V4"],
           ylim=c(0,500),
           cex.main=1.0,
           sizes=c(rep(2,length(alignment_tracks)),10)
)
dev.off()
