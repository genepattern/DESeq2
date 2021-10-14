FROM jupyter/datascience-notebook:r-4.0.3
MAINTAINER Anthony S. Castanza <acastanza@ucsd.edu>

ENV LANG=C LC_ALL=C
USER root

RUN mkdir /build

# install R dependencies
RUN R -e 'install.packages("foreign", repos = "https://cloud.r-project.org/")'
RUN R -e 'install.packages("Hmisc", repos = "https://cloud.r-project.org/")'
RUN R -e 'install.packages("getopt", repos = "https://cloud.r-project.org/")'
RUN R -e 'install.packages("optparse", repos = "https://cloud.r-project.org/")'
RUN R -e 'install.packages("BiocManager", repos = "https://cloud.r-project.org/")'
RUN R -e 'BiocManager::install()'
RUN R -e 'BiocManager::install("tximport")'
RUN R -e 'BiocManager::install("DESeq2")'
RUN R -e 'BiocManager::install("GenomicFeatures")'
RUN R -e 'BiocManager::install("rhdf5")'

# copy module files
COPY module/* /build/
RUN chmod a+x /build/tximport.deseq2.R

RUN R -e "sessionInfo()"
RUN rm -rf /tmp/downloaded_packages/

CMD ["Rscript", "--version"]

# build using this:
# docker build -t genepattern/deseq2:0.1.
