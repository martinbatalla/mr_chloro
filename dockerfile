FROM continuumio/miniconda3:latest

RUN conda install -c bioconda -c conda-forge bwa-mem2 samtools blast -y && conda clean -a -y