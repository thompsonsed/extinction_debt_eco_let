# Code examples for pycoalescence usage
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/thompsonsed/extinction_debt_eco_let/master?filepath=%2Fhome%2Fextinction_debt_eco_let%2Fexample_simulation.ipynb)
[![Documentation status](https://readthedocs.org/projects/pip/badge/?version=latest&style=flat)](https://pycoalescence.readthedocs.io)


This repository contains examples of the pycoalescence package described [here](https://pycoalescence.readthedocs.io) 
in support of Thompson *et al* (2019). Primarily, the binder jupyter notebook can be run through a browser to demonstrate the simulation process used in the paper for both spatially explicit neutral models and calculations of the 
effective connectivity parameter. Examples of map files can be found [here](10.6084/m9.figshare.c.4660199).

### File and folder overview

- **results**: folder contains the outputs of simulations used for all analyses.

- **code_examples**: folder contains scripts and jupyter notebooks used for generating the landscapes found [here](10.6084/m9.figshare.c.4660199), performing the dispersal simulations on HPC systems and performing coalescence simulations on HPC systems.

- **figures**: folder containing all figures present in the main text.

- **example\_simulation.ipynb**: run simulations through MyBinder as an example of how models were performed.

- **figure\_generation.R**: generate the figures shown in the main text using simulation results stored in csv files in
  *results* folder.
    
- **plot\_colours.R**: controls the colours and labels for the plots.

- **preston.R**: contains functions for generating the analytical solutions from the Preston function.



## Features
#### Generation of map files

Maps are split into three types:

- Clustered maps (also called "contrived maps") involve near-circular islands of habitat see *ContrivedMapGeneration.ipynb* for examples.
- Random maps were generated using random noise (see "fragment_generation.py")
- Real maps were extracted from data (see "fragment_generation.py")

#### Coalescence simulations

Coalescence simulations were performed using the *pycoalescence* package (detailed [here](https://pycoalescence.readthedocs.io/en/develop/#)) and run on HPC systems at Imperial College London. Examples of code can be found in *example\_simulation.ipynb* and in *code_examples/coalescence\_simulations*.

### Technical details

[MyBinder](https://mybinder.org/) is used to generate an interactive [jupyter notebook](http://jupyter.org/) containing
the code for running the example simulations and analyses.

The docker image used to generate the binder notebook is available 
[here](https://hub.docker.com/r/thompsonsed/pycoalescence-circleci-0.0.1/) or installable using
 `docker pull thompsonsed/pycoalescence-circleci-0.0.1` 


