# Differential gene expression analysis

Follow these instructions to recreate the differential gene expression analyses
described in the study.

## Requirements

The following software needs to be available on your system:

- GNU bash 5.0.17(1)
- GNU coreutils 8.3.0
- `awk` (GNU Awk 4.0.2)
- `bedtools` (v2.26.0)
- `gzip` (v1.10)
- `perl` (v5.16.3)
- `R` (v3.5.0)
  - package `edgeR` (v3.34.0)
- `sed` (GNU sed 4.2.2)

> Note: Indicated versions are the ones that were used for the analyses as
> described in the study. Other versions, particularly those with different
> major version numbers, may lead to errors or different results.

## Preparing the annotations

Firstly, download human gene annotations from Ensembl (release 90):

```bash
wget \
  -O annotations/Homo_sapiens.GRCh38.90.chr.gtf \
  http://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.chr.gtf.gz
```

> If gene annotations are not available at
> `annotations/Homo_sapiens.GRCh38.90.chr.gtf`, relative to _this_ directory,
> you will need to edit the script called in the next step.

Then filter gene annotations by executing:

```bash
bash ./docs/01.filter_annotations.sh
```

Finally, discard 3'-terminal exons and compute exclusively "exonic" regions:

```bash
bash ./docs/02.get_exonic_intronic.sh
```

## Intersecting alignments with annotations

First, read libraries have to be processed as described in folder
`../Preprocessing`.

> The resulting sorted, indexed alignments (`.bam` and `.bam.bai` files need to
> be available in directory `alignments/`, relative to _this_ directory. If
> not, the script called in the next step needs to be modified accordingly.

Uniquely mapping alignments are then intersected with the prepared annotations:

```bash
bash ./docs/03.count_exonic_intronic_reads.sh
```

The resulting count tables, one for each library, are then merged into a single
table:

```bash
Rscript ./docs/04.summarize_exonic_intronic_reads.R
```

## Determining differentially expressed genes

Finally, `edgeR` is run on the merged count table to determine differentially
expressed genes across four different sample group comparisons:

```bash
Rscript ./docs/05.dgea.R
```

## Results

All results are going to be written to directory `results/`, relative to _this_
directory. It contains two subdirectories `counts/` and `dgea/`. The most
important files produced during execution of the commands listed above are:

- `annotations/hsa.GRCh38_90.gene_annotations.exonic_intronic.gtf.gz`: Fully
  processed set of gene annotations
- `results/counts/merged.all_experiments.counts`: Count table produced after
  intersecting alignments of all libraries with the processed gene annotations
- `results/dgea/dgea.count_table.tab`: The final count table further filtered
  for expression; used as basis for all differential gene expression analyses
- `results/dgea/dgea.1_wildtype.vs.2_*.tab`: Tables listing differentially
  expressed genes for each of four different comparisons

## Checkpoints

Note that the fully processed annotations as well as the complete set of
results from the differential gene expression analyses are available in
the `annotations/` and `results/dgea/` directories respectively for your
convenience. You can use these files as checkpoints or to compare final
results. The precise file paths, relative to _this_ directory, are as follows:

- `annotations/processed_annotations.tar.gz`
- `results/results/dgea/dgea_results.tar.gz`
