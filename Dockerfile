FROM jupyter/scipy-notebook:latest


RUN conda install --quiet --yes --channel https://conda.anaconda.org/conda-forge rdkit && conda clean -tipsy
RUN pip3 install sygma 
