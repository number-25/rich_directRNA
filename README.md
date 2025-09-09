[![nf-core CI](https://github.com/number-25/rich_directRNA/actions/workflows/ci.yml/badge.svg?branch=dev)](https://github.com/number-25/rich_directRNA/actions/workflows/ci.yml)
[![nf-core linting comment](https://github.com/number-25/rich_directRNA/actions/workflows/linting_comment.yml/badge.svg)](https://github.com/number-25/rich_directRNA/actions/workflows/linting_comment.yml)
[![GitHub Actions Linting Status](https://github.com/number-25/rich_directRNA/actions/workflows/linting.yml/badge.svg)](https://github.com/number-25/rich_directRNA/actions/workflows/linting.yml)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.17082314-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.10.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/number-25/rich_directrna)

## Introduction

**rich_directRNA** is a bioinformatics pipeline that is still in the works… It
is a nextflow pipeline that is used for the processing of direct RNA nanopore sequencing data, providing multiple transcript reconstruction, and quantification
options with the use of a reference genome, and transcriptome annotation. Additionally, it performs post transcriptome reconstruction assessment, and recovery.

The pipeline currently _only_ accepts sequencing data from directRNA Oxford
Nanopore Technologies (ONT) libraries. It is recommended to provide raw FASTQ
files to the pipeline, however, it will also accept already mapped sequencing
reads in BAM format. These are provided to the samplesheet as input.

<!-- nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->

1. QC of FASTQ input files ( [`NANOQ`](https://github.com/esteinig/nanoq), [`SEQUALI`](https://github.com/rhpvorderman/sequali) )
2. Mapping to reference genome ( [`minimap2`](github.com/lh3/minimap2) )
3. Sort and index alignments ( [`samtools`](https://sourceforge.net/projects/samtools/files/samtools/) )
4. Create bigWig coverage files ( [`bedtools`](https://github.com/arq5x/bedtools2/), [bedGraphToBigWig](https://hgdownload.soe.ucsc.edu/admin/exe/) )
5. Extensive QC of alignments
   1. [`samtools`](https://sourceforge.net/projects/samtools/files/samtools/)
   2. [`cramino`](https://github.com/wdecoster/cramino)
   3. [`alfred`](https://www.gear-genomics.com/docs/alfred/)
   4. [`ngs-bits`](https://github.com/imgag/ngs-bits/tree/master)
6. Multiple transcriptome reconstruction options, with read correction options.
   1. [`FLAIR`](github.com/BrooksLabUCSC/flair) - allows read correction
   2. [`bambu`](github.com/GoekeLab/bambu) - very minor read correction
   3. [`IsoQuant`](https://ablab.github.io/IsoQuant/) - allows read correction
   4. [`StringTie`](https://github.com/skovaka/stringtie2)
7. Fusion gene detection [`JAFFA`](github.com/Oshlack/JAFFA)
8. Transcriptome assessment [`gffutils`](https://ccb.jhu.edu/software/stringtie/gff.shtml)
9. Transcript quantification ( [`TranSigner`](https://github.com/haydenji0731/TranSigner), [oarfish](https://github.com/COMBINE-lab/oarfish) )

Small test datasets for the pipeline are included in the [assets directory](https://github.com/number-25/rich_directRNA/assets/test_data).

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

<!--  nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

Download the reference genome, transcriptome etc. to be used by the pipeline.
We will use the hg38 analysis set and Ensemble 112, corresponding to GENCODE
release 46.
Navigate and enter into the `assets` directory and execute `bash download_references.sh`.
Ensure that you have gunzip and rsync installed on your system.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

-->

Now, you can run the pipeline using:

<!-- nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run . \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

To run a minimal, quick test dataset, use:

```bash
mkdir testing_dir

nextflow run . \
    -profile test,singularity \
    --outdir testing_dir \
    -c conf/test.config`
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Credits

number-25/rich_directRNA was originally written by Dean Bašić.

We thank the following people for their extensive assistance in the development of this pipeline:

<!--  nf-core: If applicable, make list of people who have also contributed -->

## Citations

<!--  nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use number-25/rich_directRNA for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
