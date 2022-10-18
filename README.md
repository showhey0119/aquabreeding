# Aquabreeding

Aquabreeding software is a python tool for simulating aquaculture breeding.

Aquabreeding generates a founder population using coalescent simulation implemented in [msprime](https://tskit.dev/msprime/docs/stable/intro.html).  A progeny population is produced by crossing the founder individuals, and the progenies with larger breeding/phenotypic values are selected as the founder of the next generation.  Phenotypic values, breeding values, variance components, and inbreeding coefficient are calculated in each generation.

This is a beta version.  Please do not use it for any publications.


## Requirements
- macOS (No guarantee that it will work on unix)
- python3  
  (e.g., brew install python3)
- python libraries (automatically installed)
    - numpy
    - scipy
    - numba  
    - [msprime](https://tskit.dev/msprime/docs/stable/intro.html)  
    - pandas (for jupyter notebook)  
    - matplotlib (for jupyter notebook)  
    - seaborn (for jupyter notebook)  
- cmake  
  (e.g., brew install cmake)
- C++ library
    - eigen  
      (e.g., brew install eigen)


## Installation
`git clone --recursive https://github.com/showhey0119/aquabreeding`  
`git submodule update --init --recursive`  
`cd ./aquabreeding`  
`pip3 install .`  


## Quick startup
see [aquabreeding\_tutorial.ipynb](https://github.com/showhey0119/aquabreeding/blob/master/aquabreeding_tutorial.ipynb)

## History
- Version 0.7.0 is released  
-- Sex is implemented  
-- Mating design can be flexibly set  
-- ABLUP and GBLUP are implemented  

- Version 0.6.0 is released
