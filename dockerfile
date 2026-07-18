FROM continuumio/miniconda3:latest

RUN apt-get update && apt-get install -y libgomp1 && apt-get clean

RUN conda install -c bioconda -c conda-forge bwa-mem2 samtools blast -y && conda clean -a -y