# BASE IMAGE
FROM bioconductor/bioconductor_docker:3.17

USER root
RUN mkdir /mnt/genepatt
WORKDIR /mnt/genepatt
COPY src/* /mnt/genepatt


RUN R -e "install.packages('getopt', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('optparse', repos='http://cran.us.r-project.org')"
RUN R -e "BiocManager::install('DESeq2')"