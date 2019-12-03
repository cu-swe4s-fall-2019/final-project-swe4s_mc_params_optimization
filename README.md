ParamBayes
==============================
[//]: # (Badges)
[![Travis Build Status](https://travis-ci.org/REPLACE_WITH_OWNER_ACCOUNT/ParamBayes.svg?branch=master)](https://travis-ci.org/REPLACE_WITH_OWNER_ACCOUNT/ParamBayes)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/ParamBayes/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/ParamBayes/branch/master)

Bayesian MCMC parameterization for CSCI 7000 (swe4s) class project

## Installation

This package uses the `conda` package. Installation instructions can be found at [conda's website](https://docs.anaconda.com/anaconda/install/). 

Create a dedicated python environment for parambayes using the following command: `conda create --yes -n parambayes`

Activate parambayes environment command: `conda activate parambayes`

Conda install the following packages:

* matplotlib
* pandas
* scipy
* numpy 
* tqdm
* pyyaml


## Usage

1. Navigate to [parambayes](https://github.com/cu-swe4s-fall-2019/final-project-swe4s_mc_params_optimization/tree/master/parambayes).
2. Run `python create_runfile.py`. This creates a runfile in [runfiles](https://github.com/cu-swe4s-fall-2019/final-project-swe4s_mc_params_optimization/tree/master/parambayes/runfiles) folder. 
3. Finally, run `python mcmc_run_yaml.py -f <full path to runfile>`:
    1. This creates an `output` file folder with relevant graphs, as well as a copy of the runfile. 
    2. Multiple runs using different runfiles will be stored in the output folder. 


### Copyright

Copyright (c) 2019, Owen Madin, Lenny Fobe, Ryan Morelock


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
