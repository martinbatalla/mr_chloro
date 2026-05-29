FROM continuumio/miniconda3:latest

RUN conda install -c bioconda -c conda-forge bwa-mem2 samtools blast -y && conda clean -a -y

COPY std_modified.sh /usr/local/bin/std_modified.sh

RUN chmod +x /usr/local/bin/std_modified.sh