![](./docs/images/nanortax_logo.png)

**Real-time analysis pipeline for nanopore 16S rRNA data**.



<img src="./docs/images/nanortax_workflow.png" style="zoom: 40%;" />

## Introduction

NanoRTax is a taxonomic and diversity analysis pipeline built originally for Nanopore 16S rRNA data with real-time analysis support in mind. It combines state-of-the-art classifiers such as Kraken2, Centrifuge and BLAST with downstream analysis steps to provide a framework for the analysis of in-progress sequencing runs. NanoRTax retrieves the final output files in the same structure/format for every classifier which enables more comprehensive tool/database comparison and better benchmarking capabilities. Additionally, NanoRTax includes a web application (./viz_webapp/) for visualizing complete or partial pipeline outputs. 


The NanoRTax pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with conda environments and docker containers making installation trivial and results highly reproducible.

## Quick Start

i. Install [`nextflow`](https://nf-co.re/usage/installation)

ii. Install either [`Docker`](https://docs.docker.com/engine/installation/) for full pipeline reproducibility or use [`Conda`](https://conda.io/miniconda.html)

iii. Download the pipeline and example databases and test it on a minimal dataset with a single command

```bash
#UPDATE:Taxonomic data necessary for taxonkit
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzvf taxdump.tar.gz -C db/
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

Similar to the normal mode but using --reads_rt for input. Partial results are stored in an output directory and are accessible as .csv file and web app visualization files. In this mode, the workflow will run endlessly, so it needs to be stopped manually by Ctrl+C once all  FASTQ files are completely processed.

Note: This mode is intended to work with non-bulk FASTQ files (ie: 500 reads per file) in order to provide a fluid real-time analysis of generated reads. This aspect can be configured before starting the experiment via MinKNOW sequencing software.

v. Visualize partial/complete outputs using NanoRTax web application (./viz_webapp)

Before running the web application, make sure to have the necessary dependencies installed or use the provided viz_webapp/environment.yml file to build a conda environment (recommended):

```bash
conda env create -f environment.yml
conda activate nanortax_webapp
```
Start the **web application** server with the command below and access the interface with a web browser (http://127.0.0.1:8050/ by default).

```bash
cd viz_webapp && python dashboard.py
```

See [usage docs](docs/usage.md) for all of the available options when running the pipeline.


## Documentation

The NanoRTax pipeline comes with documentation about the pipeline, found in the `docs/` directory:

[Running the pipeline](docs/usage.md)
[Output and how to interpret the results](docs/output.md)

## Credits

Rodríguez-Pérez H, Ciuffreda L, Flores C. NanoRTax, a real-time pipeline for taxonomic and diversity analysis of nanopore 16S rRNA amplicon sequencing data. Comput Struct Biotechnol J. 2022;20:5350-5354. doi: https://doi.org/10.1016/j.csbj.2022.09.024

This work was supported by Instituto de Salud Carlos III [PI14/00844, PI17/00610, and FI18/00230] and co-financed by the European Regional Development Funds, “A way of making Europe” from the European Union; Ministerio de Ciencia e Innovación [RTC-2017–6471-1, AEI/FEDER, UE]; Cabildo Insular de Tenerife [CGIEU0000219140]; Fundación Canaria Instituto de Investigación Sanitaria de Canarias [PIFUN48/18]; and by the agreement with Instituto Tecnológico y de Energías Renovables (ITER) to strengthen scientific and technological education, training, research, development and innovation in Genomics, Personalized Medicine and Biotechnology [OA17/008].

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

