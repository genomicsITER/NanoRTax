# NanoRTax

**Real-time analysis pipeline for nanopore 16S rRNA data**.

[![Docker](https://img.shields.io/docker/automated/hecrp/nanortax.svg)](https://hub.docker.com/r/hecrp/nanortax)

## Introduction

NanoRTax is a taxonomic and diversity analysis pipeline built originally for Nanopore 16S rRNA data with real-time analysis support in mind. It combines state-of-the-art classifiers such as Kraken2, Centrifuge and BLAST with downstream analysis steps to provide a comprehensive output. NanoRTax retrieves the final output files in the same structure/format for every classifier which enables more comprehensive tool/database comparison and better benchmarking capabilities. Additionally, NanoRTax includes a web application (./viz_webapp/) for visualizing complete or in-progress pipeline runs. 


The NanoRTax pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with conda environments and docker containers making installation trivial and results highly reproducible.

## Quick Start

i. Install [`nextflow`](https://nf-co.re/usage/installation)

ii. Install either [`Docker`](https://docs.docker.com/engine/installation/) for full pipeline reproducibility or use [`Conda`](https://conda.io/miniconda.html)

iii. Download the pipeline and example databases and test it on a minimal dataset with a single command

```bash
#BLAST database
mkdir db db/taxdb
wget https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz && tar -xzvf 16S_ribosomal_RNA.tar.gz -C db
wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz && tar -xzvf taxdb.tar.gz -C db/taxdb
#Kraken2 RDP database
wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/16S_RDP11.5_20200326.tgz && tar -xzvf 16S_RDP11.5_20200326.tgz -C db
#Centrifuge P_COMPRESSED database (more information: https://ccb.jhu.edu/software/centrifuge/manual.shtml#database-download-and-index-building)
wget https://genome-idx.s3.amazonaws.com/centrifuge/p_compressed_2018_4_15.tar.gz && tar -xzvf p_compressed_2018_4_15.tar.gz -C db
```

```bash
nextflow run main.nf -profile test,<docker/conda>
```

iv. Start running your own analysis!

We provide an example configuration profile with the default parameters for running the pipeline (conf/default.config) and it is a good starting point to easily customize your NanoRTax workflow. This configuration is loaded by specifying "default" in the profiles list of pipeline command. 

a. Run classification on a single FASTQ file
```bash
nextflow run main.nf -profile <default,docker/conda> --reads '/seq_path/sample.fastq'
```
b. Run classification on an entire sequencing run directory. NanoRTax will detect the barcode directories and analyze all samples:
```bash
nextflow run main.nf -profile <default,docker/conda> --reads '/seq_path/fastq_pass/**/*.fastq'
```

c. Real-time mode.

```bash
nextflow run main.nf -profile <default,docker/conda> --reads_rt '/seq_path/fastq_pass/**/*.fastq'
```

Similar to the normal mode but using --reads_rt for input. Partial results stored at output directory and are also accesible (.csv files and webapp visualization). In this mode, the workflow will run endlessly, so it needs to be stopped manually by Ctrl+C once all consumed FASTQ files are completely processed.

Note: This mode is intended to work with non-bulk FASTQ files (ie: 500 reads per file) in order to provide a fluid real-time analysis of generated reads. This aspect can be configured before starting the experiment via MinKNOW sequencing software.

v. Visualize partial/complete outputs

NanoRTax comes with a Python Dash web application which provides interactive visualization of partial and complete results.
```bash
cd viz_webapp && python dashboard.py
```
Start the web application server with the command above and access the interface with a web browser (http://127.0.0.1:8050/ by default).

See [usage docs](docs/usage.md) for all of the available options when running the pipeline.


## Documentation

The NanoRTax pipeline comes with documentation about the pipeline, found in the `docs/` directory:

[Running the pipeline](docs/usage.md)
[Output and how to interpret the results](docs/output.md)

## Credits

NanoRTax was originally written by Laura Ciuffreda, Héctor Rodríguez Pérez and Carlos Flores.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

