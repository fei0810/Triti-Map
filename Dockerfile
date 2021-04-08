FROM continuumio/miniconda3
RUN conda create -y -c conda-forge -c bioconda -n tritimap tritimap \
&& echo "source activate tritimap" > ~/.bashrc
ENV PATH /opt/conda/envs/tritimap/bin:$PATH
