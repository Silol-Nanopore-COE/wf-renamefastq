<!-- [![GitHub Actions CI Status](https://github.com/wf/renamefastq/actions/workflows/ci.yml/badge.svg)](https://github.com/wf/renamefastq/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/wf/renamefastq/actions/workflows/linting.yml/badge.svg)](https://github.com/wf/renamefastq/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/wf/renamefastq) -->

<img src="docs/images/silol-logo.png" width="300">

# wf-renamefastq: A workflow for renaming and demultiplexing FASTQ files

## Introduction

This workflow is used for renaming and/or demultiplexing the FASTQ file obtained from Oxford Nanopore Sequencer.

The workflow will perform the following:
1. Check a comma-separated sample sheet (*.csv) if it is provided.
2. Demultiplex FASTQ files when `--demultiplex` is specified. ([Dorado](https://software-docs.nanoporetech.com/dorado/latest))
3. Renaming FASTQ files and/or concatenating mulitple FASTQ file into a single FASTQ file. ([FASTCAT](https://github.com/epi2me-labs/fastcat))
4. Filtering reads by a Phred quality score
 (Q-score) when `--quality_filter` is specified. ([Seqkit](https://bioinf.shenwei.me/seqkit/usage))
5. Computing simple statistics of demultiplexed / filtered reads. ([Seqkit](https://bioinf.shenwei.me/seqkit/usage))

## Compute requirements

Recommended requirements:
- CPUs = 16
- Memory = 64 GB

Minimum requirements:
- CPUs = 8
- Memory = 16 GB

## Input example

This workflow takes FASTQ files as input and accepts one of three cases:

1. The path to a single FASTQ
2. The path to a top-level directory containing FASTQ files
3. The path to a directory containing one level of sub-directories which in turn contain FASTQ files

```
(1)                     (2)                 (3)    
input_reads.fastq   ─── input_directory  ─── input_directory
                        ├── reads0.fastq     ├── barcode01
                        └── reads1.fastq     │   ├── reads0.fastq
                                             │   └── reads1.fastq
                                             ├── barcode02
                                             │   ├── reads0.fastq
                                             │   ├── reads1.fastq
                                             │   └── reads2.fastq
                                             └── barcode03
                                              └── reads0.fastq
```

For the first and second cases (1 and 2), a sample name can be supplied with the `--sample` parameter.

```
--sample sample_01
```

For the last case (3), the data is assumed to be multiplexed with the names of the sub-directories as barcodes. A comma-separated sample sheet can be provided with the `--sample_sheet` parameter.

```
--sample_sheet /path/to/samplesheet.csv
```

The input samplesheet used for mapping barcode to sample aliases must contain two columns named `barcode` and `alias` as shown in the example below.

```csv
barcode,alias
barcode01,A01
barcode02,A02
barcode03,A03
```
> [!IMPORTANT]
Barcodes must match with the names of sub-directories in the input FASTQ directory.

## Usage

### Renaming FASTQ file

- A single FASTQ file or top-level directory (case 1 and 2):

```bash
nextflow run wf-renamefastq \
   --fastq <path to a FASTQ file or directory> \
   --sample <sample name> \
   --outdir <OUTDIR> \
   -profile <docker>
```

- A directory containing sub-directories (case 3):

```bash
nextflow run wf-renamefastq \
   --fastq <path to a FASTQ directory> \
   --sample_sheet <sample name> \
   --outdir <OUTDIR> \
   -profile <docker>
```

### Demultiplexing FASTQ file

The input can be both a single FASTQ file and a directory containing FASTQ files.

```bash
nextflow run wf-renamefastq \
   --fastq <path to a FASTQ file or directory> \
   --sample_sheet <sample name> \
   --demultiplex \
   --outdir <OUTDIR> \
   -profile <docker>
```

<!-- > [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files). -->

## Credits

wf-renamefastq was originally written by Arissara Tubtimyoy and Piroon Jenjaroenpun.

<!-- We thank the following people for their extensive assistance in the development of this pipeline: -->

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

<!-- ## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations -->

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use wf/renamefastq for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

<!-- An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file. -->

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
