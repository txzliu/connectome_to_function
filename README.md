# connectome_to_function
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5716323.svg)](https://doi.org/10.5281/zenodo.5716323)

This repo contains Python code to interface with NEURON and simulate neurons from electron microscopy data, companion to "Connectomic features underlying diverse synaptic connection strengths and subcellular computation."

To run, use a virtual environment (i.e. Anaconda) from which you can run NEURON 7.7.2 and Python 3.7. Further, make sure to input your own user token to access the neuprint API (to load neural morphologies and synapse locations). An account can easily be created at the neuprint website (https://neuprint.janelia.org/). 

Requisite packages are imported in `nrn_cell.py`, and can be installed using a package manager like pip (i.e. `pip install numpy`). The `requirements_reference.txt` file provides package versions, but I do NOT recommend simultaneously installing all packages in the txt file with pip/conda, as it may introduce compatibility issues. Instead, use it as a reference in case version issues arise here with the import statements below. 

Please reach out at any point if questions come up :) 
