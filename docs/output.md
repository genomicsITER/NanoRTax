# NanoRTax: Output

This page describes the output produced by the pipeline. Final report and diversity files are stored in the ./results directory or a custom directory specified by --out

## Pipeline overview

**Default output directory: `results/`**

Inside output directory, sample results are organized based on barcode/sample directory name. The following files can be found there:

* `classifier_diveristy_taxlevel.csv`
  * Diversity information at every step (individual read file) of the analysis for the specific classifier tool and taxonomic level. Includes Shannon index, Simpson index and taxa count.
* `classifier_report_taxlevel.csv`
  * Simple report containing tax name (retrieved from Taxonkit) and read count for the specific classifier tool and taxonomic level. 
* `classifier_report_full.csv`
  * Simple report containing the complete taxonomic information extracted with taxonkit for every read in the dataset. Also available in OTU format.
* `qc_report.csv`
  * Complete QC report file.

Note that diversity and classification reports are created for all included classifiers in the workflow and for class, order, family, genus and species level.

## Web Application

All the above output files are also consumed by the web application, providing an interactive visualization mode for results inspection

Get the conda environment for the webapp and start the server:
```bash
cd viz_webapp && python dashboard.py
```
Access the interface with a web browser (http://127.0.0.1:8050/ by default).
