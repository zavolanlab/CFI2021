# RNA-Seq data preprocessing

All the FASTQ files containing RNA-Seq data were processed according to the templates below.

## 1. Reverse-complement the reads

PAQR requires that single-stranded data must have reads in sense direction.

(FASTX Toolkit 0.0.14)

```
gunzip -c {path_to_the_fastq_input_file} | fastx_reverse_complement -z -o {path_to_the_fastq_output_file}
```

## 2. Adapter removal

Adapter on the original reads: `TGGAATTCTCGGGTGCCAAGG`  
Reverse-complement: `CCTTGGCACCCGAGAATTCCA`

(cutadapt version 1.16)

```
cutadapt \
-b CCTTGGCACCCGAGAATTCCA \
-e 0.1 \
--times 1 \
--trim-n \
--cores=8 \
--minimum-length 10 \
-o {path_for_the_output_file} \
{path_to_the_input_file} \
1> stdout.log \
2> stderr.log
```

## 3. Poly(A)-tails removal

(cutadapt version 1.16)

```
cutadapt \
-b AAAAAAAAAAAAAAAAAAAA \
-b TTTTTTTTTTTTTTTTTTTT \
-e 0.1 \
-O 1 \
--times 1 \
--trim-n \
--cores=8 \
--minimum-length 10 \
-o {path_for_the_output_file} \
{path_to_the_input_file} \
1> stdout.log \
2> stderr.log
```

## 4. Read mapping

Genomic resources used:  
`Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa`  
`Homo_sapiens.GRCh38.90.chr.gtf`

(STAR aligner version 2.7.1a)

```
STAR \
--runMode alignReads \
--twopassMode Basic \
--outSAMunmapped None \
--outSAMattributes All \
--outReadsUnmapped None \
--outFilterType BySJout \
--alignEndsType Local \
--outFilterMismatchNoverLmax 0.1 \
--outFilterScoreMinOverLread 0.66 \
--outFilterMatchNminOverLread 0.66 \
--outFilterMultimapNmax 10 \
--outFilterMultimapScoreRange 0 \
--runThreadN 8 \
--genomeDir ../index \
--sjdbGTFfile Homo_sapiens.GRCh38.90.chr.gtf \
--readFilesIn {path_to_the_input_file} \
--readFilesCommand zcat \
--outFileNamePrefix {prefix_for_the_output_file} \
--outSAMtype BAM Unsorted \
--quantMode TranscriptomeSAM \
1> stdout.log \
2> stderr.log
```

## 5. Sorting and Indexing the alignments

(SAMtools version 1.10)

```
samtools sort \
-@ 8 \
{path_to_the_input_file} \
1> {path_for_the_output_file} \
2> stderr.log
```

```
samtools index \
-@ 8 \
{path_to_the_input_file} \
1> {path_for_the_output_file} \
2> stderr.log
```
