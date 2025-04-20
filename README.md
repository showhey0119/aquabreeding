# Aquabreeding

Aquabreeding software is a python tool for simulating aquaculture breeding.

Aquabreeding generates a founder population using coalescent simulation implemented in [msprime](https://tskit.dev/msprime/docs/stable/intro.html).  A progeny population is produced by crossing the founder individuals, and the progenies with larger breeding/phenotypic values are selected as parents of the next generation.  Phenotypic values, breeding values, variance components, and inbreeding coefficient are calculated in each generation.

This is a beta version.  Please do not use it for any publications.


## Last update
- 12/29/2023


## Note
- Aquabreeding supports python3.11.  


## Requirements
- macOS, Linux
- python3 (3.11)  
  (e.g., brew install python3)
- python libraries (automatically installed)
    - numpy
    - scipy
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
`cd ./aquabreeding`  
`git submodule update --init --recursive`  
`pip3 install .`  

If the installation fails, try
`pip3 install -e .`  

Trouble shooting: the cpp library sometimes cannot be built.  Updating pybind11 manually may solve this problem.  
`git submodule deinit -f pybind11`  
`git rm -f pybind11`  
`rm -rf .git/modules/pybind11`  
`git submodule add -b master https://github.com/pybind/pybind11.git pybind11`  
`pip3 install .`  

## Quick startup
see [aquabreeding\_tutorial.ipynb](https://github.com/showhey0119/aquabreeding/blob/master/aquabreeding_tutorial.ipynb)


## History
- Version 0.9.0 is released (12/29/2023)  
    - Aquabreeding supports python3.11  

- Version 0.8.1 is released (05/19/2023)  
    - Bug fix  

- Version 0.8.0 is released (03/22/2023)  
    - Breeding parents are sampled either from a standar Wright-Fisher population or structured population.

- Version 0.7.2 is released (10/26/2022)  
    - New mating design that minimizes inbreeding  

- Version 0.7.1 is released (10/19/2022)  
    - Family selection is implemented  

- Version 0.7.0 is released (10/18/2022)  
    - Sex is implemented  
    - Mating design can be flexibly set  
    - ABLUP and GBLUP are implemented  

- Version 0.6.0 is released
