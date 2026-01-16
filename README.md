[![nf-core CI](https://github.com/number-25/LongTranscriptomics/actions/workflows/ci.yml/badge.svg?branch=dev)](https://github.com/number-25/LongTranscriptomics/actions/workflows/ci.yml)
[![nf-core linting comment](https://github.com/number-25/LongTranscriptomics/actions/workflows/linting_comment.yml/badge.svg)](https://github.com/number-25/LongTranscriptomics/actions/workflows/linting_comment.yml)
[![GitHub Actions Linting Status](https://github.com/number-25/LongTranscriptomics/actions/workflows/linting.yml/badge.svg)](https://github.com/number-25/LongTranscriptomics/actions/workflows/linting.yml)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.17082314-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.17082314-1073c8)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.10.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/number-25/rich_directrna)

## Introduction

**LongTranscriptomics** is a nextflow pipeline that is used for the processing of direct RNA nanopore sequencing data, providing multiple transcript reconstruction options, and quantification with the use of a reference genome, and transcriptome annotation.

<!-- Additionally, it performs post transcriptome reconstruction assessment, and recovery. -->

The pipeline accepts sequencing data from both directRNA and cDNA Oxford
Nanopore Technologies (ONT) libraries. It is recommended to provide raw FASTQ
files to the pipeline, however, it will also accept already mapped sequencing
reads in BAM format. These are provided to the samplesheet as input.

The general flow of the pipeline is as follows;

<!-- nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->

1. QC of FASTQ input files ( [`NANOQ`](https://github.com/esteinig/nanoq), [`SEQUALI`](https://github.com/rhpvorderman/sequali) )
2. Mapping to reference genome ( [`minimap2`](github.com/lh3/minimap2) )
3. Sort and index alignments ( [`samtools`](https://sourceforge.net/projects/samtools/files/samtools/) )
4. Create bigWig coverage files ( [`bedtools`](https://github.com/arq5x/bedtools2/), [`bedGraphToBigWig`](https://hgdownload.soe.ucsc.edu/admin/exe/) )
5. Extensive QC of alignments
   1. [`samtools`](https://sourceforge.net/projects/samtools/files/samtools/)
   2. [`cramino`](https://github.com/wdecoster/cramino)
   3. [`alfred`](https://www.gear-genomics.com/docs/alfred/)
   4. [`ngs-bits`](https://github.com/imgag/ngs-bits/tree/master)
6. Multiple transcriptome reconstruction options, with read correction options.
   1. [`FLAIR`](https://github.com/BrooksLabUCSC/flair) - allows read correction
   2. [`bambu`](http://github.com/GoekeLab/bambu) - very minor read correction
   3. [`IsoQuant`](https://ablab.github.io/IsoQuant/) - allows read correction
   4. [`StringTie`](https://github.com/skovaka/stringtie2)
   <!-- 7. Fusion gene detection [`JAFFA`](github.com/Oshlack/JAFFA) -->
7. Transcriptome assessment ( [`gffcompare`](https://ccb.jhu.edu/software/stringtie/gff.shtml) )
8. Transcript quantification ( [`oarfish`](https://github.com/COMBINE-lab/oarfish) )
     <!-- ( [`TranSigner`](https://github.com/haydenji0731/TranSigner),
   Small test datasets for the pipeline are included in the [assets directory](https://github.com/number-25/LongTranscriptomics/assets/test_data). -->

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,replicate,reads
CONTROL1,1,data/long_reads_1.fastq.gz
CONTROL1,2,data/long_reads_2.fastq.gz
```

Each row represents a fastq file. Replicate refers to a technical replicate, biological replicates should be named uniquely. Be sure to pay attention to sample naming, in
order to avoid duplication and file overwriting. The replicate field is optional, the other two are mandatory.

The basic reference files required to run the pipeline are 1) a genome in fasta format, and 2) a transcriptome annotation in gtf format. It is advised that the files be gzipped, but it is fine if they are not.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

Once you've created the samplesheet, acquired a genome fasta and a transcriptome annotation, you can run the pipeline using:

```bash
nextflow run . \
   -profile <docker/singularity/.../institute> \
   -c <CONFIG FILE> \
   --input <SAMPLESHEET> \
   --outdir <OUTDIR> \
   --genome_fasta <GENOME FASTA> \
   --annotation_gtf <ANNOTATION GTF>
```

To run a minimal, quick test dataset, use:

```bash
mkdir testing_dir

nextflow run . \
    -profile test,singularity \
    --outdir testing_dir \
    -c conf/test.config`
```

For additional documentation on usage of the workflow, and details on outputs, please refer to the usage [documentation](./docs/).

## Credits

number-25/LongTranscriptomics was originally written by Dean Bašić.

## Citations

If you use number-25/LongTranscriptomics for your analysis, please cite it using the following doi: [10.5281/zenodo.17082314-1073c8](https://doi.org/10.5281/zenodo.17082314-1073c8).

<!-- nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
