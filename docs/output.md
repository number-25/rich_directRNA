# number-25/LongTranscriptomics: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

Directories corresponding to the stages listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [FASTQ quality control and summary stats](#FASTQ-quality-control)
  - [NANOQ](#NANOQ)
  - [SEQUALI](#SEQUALI)
- [Reference genome mapping](#Reference-genome-mapping)
  - [minimap2](#minimap2)
  - [samtools](#samtools-sort-index)
- [Create bigWig coverage files](#Create-files-to-visualise-mapping)
  - [bedtools](#bedtools)
  - [bedGraphToBigWig](#bedGraphToBigWig)
- [Extensive QC of alignments](#Alignment-quality-control)
  - [samtools](#samtools-flagstat)
  - [cramino](#cramino)
  - [alfred](#alfred)
  - [ngs-bits](#ngs-bits)
- [Transcriptome reconstruction](#Transcriptome-reconstruction)
  - [FLAIR](#FLAIR)
  - [bambu](#bambu)
  - [IsoQuant](#IsoQuant)
  - [StringTie](#StringTie)
<!-- 7. Fusion gene detection [`JAFFA`](github.com/Oshlack/JAFFA) -->
- [Transcriptome assessment](#Transcriptome-assessment)
  - [gffutils](#gffutils)
- [Transcript quantification](#Transcript-quantification)
  - [TranSigner](#TranSigner)
  - [oarfish](#oarfish)
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## FASTQ quality control

### NANOQ

<details markdown="1">
<summary>Output files</summary>

- `fastq_qc/nanoq/`
  - `*_nanoq.json`: `json` formatted file containing quality metrics.
  - `*_nanoq.stats`: basic NANOQ report containing quality metrics.
  - `*_nanoq_stats.verbose`: verbose NANOQ report containing quality metrics.

</details>

[NANOQ](https://github.com/esteinig/nanoq) provides general quality statistics
about the nanopore sequence reads. It outputs the statistics in both verbose and
minimal reports, which can be formatted in `json` format.

```
Nanoq Read Summary
====================

Number of reads:      100000
Number of bases:      400398234
N50 read length:      5154
Longest read:         44888
Shortest read:        5
Mean read length:     4003
Median read length:   3256
Mean read quality:    NaN
Median read quality:  NaN


Read length thresholds (bp)

> 200       99104             99.1%
> 500       96406             96.4%
> 1000      90837             90.8%
> 2000      73579             73.6%
> 5000      25515             25.5%
> 10000     4987              05.0%
> 30000     47                00.0%
> 50000     0                 00.0%
> 100000    0                 00.0%
> 1000000   0                 00.0%


Top ranking read lengths (bp)

1. 44888
2. 40044
3. 37441
4. 36543
5. 35630
```

### SEQUALI

<details markdown="1">
<summary>Output files</summary>

- `fastq_qc/sequali/`
  - `*_sequali.json`: `json` formatted file containing quality metrics.
  - `*_sequali.html`: `html` formatted containing quality metrics.

</details>

[SEQUALI](https://github.com/rhpvorderman/sequali) provides general quality statistics
about the sequence reads, along with several other features including,
overrepresentation analysis and duplication rate estimation. It outputs the
statistics in both ??

## Reference genome mapping

### minimap2

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[minimap2](https://github.com/lh3/minimap2)  is perhaps the most popular
long-read sequence aligner. In general, it aligns the sequence reads to the reference
genome/transcriptome provided by the user. Taken directly from the developers
> Minimap2 is a versatile sequence alignment program that aligns DNA or mRNA
  sequences against a large reference database. Typical use cases include: (1)
  mapping PacBio or Oxford Nanopore genomic reads to the human genome; (2)
  finding overlaps between long reads with error rate up to ~15%; (3)
  splice-aware alignment of PacBio Iso-Seq or Nanopore cDNA or Direct RNA reads
  against a reference genome; (4) aligning Illumina single- or paired-end reads;
  (5) assembly-to-assembly alignment; (6) full-genome alignment between two
  closely related species with divergence below ~15%.

### samtools sort index

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[samtools](http://www.htslib.org/doc/#manual-pages) is a multipurpose
toolkit for working with SAM/BAM files. It is used to sort the output from
minimap2 (SAM format) and output it in compressed BAM format, and then index
this file.

## Create files to visualise mapping

### bedtools

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[bedtools](https://github.com/arq5x/bedtools2) is a multipurpose
toolkit for working with tab separated genomic formats such as GTF/GFF/BED, but
also SAM/BAM/CRAM files. Here it is used convert the mapped BAM file to BEDGRAPH
format, in preparation for conversion to BigWig.

### bedGraphToBigWig

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[bedGraphToBigWig](https://hgdownload.soe.ucsc.edu/admin/exe/) is a specific
tool that is part of a broad UCSC software suite. It has one specific
function that can be guessed from it's very name. You guessed it, it converts a
bedgraph to a BigWig file, that's it. Once created, the BigWig files can be
loaded into a genome browser such as IGV, allowing the mapping to be visualised
in a lightweight way.

## Alignment quality control

### samtools flagstat

[samtools](http://www.htslib.org/doc/#manual-pages) flagstats provides summary
statistics on the mapped BAM file. Specifically, it counts the number of
alignments for each FLAG type.

### cramino

[cramino](https://github.com/wdecoster/cramino) is a tool for quick quality assessment of cram and bam files, intended for long read sequencing.

```
File name       example.cram
Number of reads 14108020
% from total reads  83.45
Yield [Gb]      139.91
N50     17447
Median length   6743.00
Mean length     9917
Median identity 94.27
Mean identity   92.53
Path    alignment/example.cram
Creation time   09/09/2022 10:53:36
```

### alfred

[alfred](https://www.gear-genomics.com/docs/alfred/cli/) computes various
alignment metrics and summary statistics by read group.

### ngs-bits

[ngs-bits
mappingQC](https://github.com/imgag/ngs-bits/blob/master/doc/tools/MappingQC/index.md)
provides one more technique for quality control of the mapped BAM files. It's
advantage is that it has an output that is compatible with
[MultiQC](https://github.com/MultiQC/MultiQC/blob/main/docs/markdown/modules/ngsbits.md).


## Transcriptome reconstruction

### FLAIR

[FLAIR](https://github.com/BrooksLabUCSC/flair) **F**ull **L**ength
**A**lternative **I**soform analysis of **R**NA is used for the correction,
isoform definition, and alternative splicing analysis of noisy reads. FLAIR has
primarily been used for nanopore cDNA, native RNA, and PacBio sequencing reads.
FLAIR is able to be used with and without read correction, making it amenable to
sensitive sample types, such as those coming from cancer where errors may
instead be putative variants which should not be corrected.

![FLAIR - example schematic](images/flair_workflow_compartmentalized.png)

### bambu

### IsoQuant

### StringTie



<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
