# Code examples for pycoalescence usage
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/thompsonsed/extinction_debt_eco_let/master?filepath=%2Fhome%2Fextinction_debt_eco_let%2Fexample_simulation.ipynb)
[![Documentation status](https://readthedocs.org/projects/pip/badge/?version=latest&style=flat)](https://pycoalescence.readthedocs.io)


This repository contains examples of the pycoalescence package described [here](https://pycoalescence.readthedocs.io) 
in support of Thompson *et al* (2018, in prep.). Primarily, the binder jupyter notebook can be run through a browser to
demonstrate the simulation process used in the paper for both spatially explicit neutral models and calculations of the 
effective connectivity parameter. Additional simulation data and collections of maps used to generate the results can
be found [here](doi::10.5281/zenodo.1493944).

### File overview

- **example\_simulation.ipynb**: run simulations through MyBinder as an example of how models were performed.

- **figure\_generation.R**: generate the figures shown in the main text using simulation results stored in csv files in
  *results* folder.
    
- **plot\_colours.R**: controls the colours and labels for the plots.

- **preston.R**: contains functions for generating the analytical solutions from the Preston function.



### Technical details

[MyBinder](https://mybinder.org/) is used to generate an interactive [jupyter notebook](http://jupyter.org/) containing
the code for running the example simulations and analyses.

The docker image used to generate the binder notebook is available 
[here](https://hub.docker.com/r/thompsonsed/pycoalescence-circleci-0.0.1/) or installable using
 `docker pull thompsonsed/pycoalescence-circleci-0.0.1` 
