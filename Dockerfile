FROM nfcore/base:1.9
LABEL authors="Laura Ciuffreda, Héctor Rodríguez Pérez" \
      description="Docker image containing all software requirements for the NanoRTax pipeline"
USER root

# Install the conda environment
COPY environment.yml /
RUN mkdir /.taxonkit && wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz && tar -xzvf taxdump.tar.gz -C /.taxonkit
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nanortax/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nanortax > nanortax.yml
