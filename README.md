# connectome_to_function
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5716323.svg)](https://doi.org/10.5281/zenodo.5716323)

This repo contains Python code to interface with NEURON and simulate neurons (i.e. mEPSPs, uEPSPs, input resistance measurements) from electron microscopy data (i.e. neural morphology SWC files, can be Drosophila neurons from neuprint). This is a companion to **Liu, Davoudian, Lizbinski, Jeanne 2021**: *Connectomic features underlying diverse synaptic connection strengths and subcellular computation*. 

### How to use the code

In general, this code contains functionality to: 
- from Python, instantiate a neural morphology file (in the widely-used [SWC format](http://www.neuronland.org/NLMorphologyConverter/MorphologyFormats/SWC/Spec.html)) as a cell object in the NEURON interface 
  - can also use synapse locations to instantiate synapses onto the cell object
  - for more help, check out the tutorials on [using Python with NEURON](https://neuron.yale.edu/neuron/docs/scripting-neuron-basics) 
  - notably, the morphology files are not restricted to being from Drosophila only! Any SWC file, i.e. from [online databases](http://neuromorpho.org/index.jsp) should be importable. In the supplementary figures, we also used related code to simulate mammalian Purkinje and pyramidal cells. 
- once the object is in NEURON, the user can manipulate it in the built in graphical user interface (i.e. look at the morphology) and also conduct functional simulations, such as by activating assortments of synapses and measuring the depolarization at any location on the neuron 
  - the user can also conduct manipulations to the neuron, such as by changing morphological features (i.e. lengthening or shortening neurites) or biophysical properties (i.e. making certain compartments isopotential by setting their axial resistances to close to zero), and understand how those manipulations affect functional properties
- for the present work, the simulations are purely passive, but active properties can also be added in NEURON

The notebook `generate attrs per LHN synapse.ipynb` will walk through the simulation of mEPSPs (activating a single synapse at a time, and also probing other synapse properties, such as local input resistances, distances from morphology landmarks, and more) onto the list of lateral horn neurons (higher-order Drosophila olfactory neurons) given in `21-05-07_LHN_SWCPointNo_to_NEURON.csv`. This data is used for Figures 5 and 6 of the paper. 

### Installation

To run, use a virtual environment (i.e. Anaconda) from which you can run NEURON 7.7.2 and Python 3.7. Further, make sure to input your own user token to access the neuprint API (to load neural morphologies and synapse locations). An account can easily be created at the neuprint website (https://neuprint.janelia.org/). 

Requisite packages are imported in `nrn_cell.py`, and can be installed using a package manager like pip (i.e. `pip install numpy`). The `requirements_reference.txt` file provides package versions, but I do NOT recommend simultaneously installing all packages in the txt file with pip/conda, as it may introduce compatibility issues. Instead, use it as a reference in case version issues arise here with the import statements below. 

---

Please reach out at any point if questions or comments come up to txliu@stanford.edu and james.jeanne@yale.edu. 
