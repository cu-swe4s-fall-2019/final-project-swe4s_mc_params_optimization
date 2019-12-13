ParamBayes
==============================
[//]: # (Badges)
[![Build Status](https://travis-ci.com/cu-swe4s-fall-2019/final-project-swe4s_mc_params_optimization.svg?branch=master)](https://travis-ci.com/cu-swe4s-fall-2019/final-project-swe4s_mc_params_optimization)
[![codecov](https://codecov.io/gh/cu-swe4s-fall-2019/final-project-swe4s_mc_params_optimization/branch/master/graph/badge.svg)](https://codecov.io/gh/cu-swe4s-fall-2019/final-project-swe4s_mc_params_optimization/branch/master)

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
* pathlib
* pycodestyle

## Usage

### Description

The Parambayes software package uses a Bayesian inference scheme to predict force field parameters for molecular dynamics (MD) simulations. This data-driven paradigm uses a Monte Carlo algorithm to explore the parameter space (energy, size, bond distance and quadrupole moment) of a subset of select two-centered molecules (O2, N2, Br2, C2H2, C2H4, C2H6, C2F4 and F2.) The parameter distributions generated by the Monte Carlo algorithm should best reproduce the experimentally observed physical properties (density, saturation pressure, and surface tension) of the molecule of interest. The data used in this optimization procedure comes from NIST (National Institute of Standards and Technology.)  

### Parameter Optimization

1. Navigate to [parambayes](https://github.com/cu-swe4s-fall-2019/final-project-swe4s_mc_params_optimization/tree/master/parambayes).
2. Run `python create_runfile.py`. This creates a runfile in [runfiles](https://github.com/cu-swe4s-fall-2019/final-project-swe4s_mc_params_optimization/tree/master/parambayes/runfiles) folder. 
3. Finally, run `python mcmc_run_yaml.py -f <full path to runfile>`:
    1. This creates an `output` file folder with relevant graphs, as well as a copy of the runfile. 
    2. Multiple runs using different runfiles will be stored in the output folder.

### Creating a Runfile

The .yaml runfile format allows for fast and easy input parameter modification. An example runfile exists [here](https://github.com/cu-swe4s-fall-2019/final-project-swe4s_mc_params_optimization/blob/master/parambayes/runfiles/C2H2_rhol%2BPsat_test_2019-11-20.yml). Input parameters include:
* compound
* label
* number_data_points
* priors
    * Function and shape parameters for energy, size, bond distance and quadrupole moment prior distributions
* properties
    * Experimentally-tabulated values (density, saturation pressure, surface tension) to use as data
* save_traj
* steps
    * number of Monte Carlo steps
* trange 
    * range of temperature to consider

### Manifest

`parambayes` contains the following files, which serve the following functions:

`parambayes.py` contains the `MCMC_Simulation` and `MCMC_Prior` classes.  
`MCMC_Simulation` creates and runs the MCMC simulations, while `MCMC_Prior` sets up the prior functions that are used in calculating the posterior.

`mcmc_run_yaml.py` contains the driver which takes a runfile and creates and runs an MCMC simulation.

`create_runfile.py` contains an application to create runfiles for use in `mcmc_run_yaml.py`

`output_models.py` contains an implementation of the surrogate models that interfaces with the `MCMC_Simulation` object.

`data_import.py` is used by `MCMC_Simulation` to import the data and parameters from files.

`plotting.py` creates the triangle plots shown here.

`utility.py` contains various auxiliary functions.

# Science

## Background

### Force fields
Molecular Dynamics (the atomistic-level simulation of molecules) has been used to study the behavior of fluids, proteins, drug binding, and other phenomena since the late 1970's.  Crucial to molecular dynamics are *force fields*, which are energy functions that encode the physical interactions of atoms into simple physical forms.  Within a specific functional form, force fields are parameterized (i.e. constants are chosen for particular chemical species so that the potential matches reality).  Force fields are generally broken into two components: 

![Nonbond](output/figures/nonbond.png)

While bonded interactions are typically fit to quantum mechanical simulations like DFT, non-bonded interactions have a less obvious connection to data and are generally fit to physical property data, making their parameterization a more difficult problem.

### Bayesian Parameter Inference

An attractive paradigm for parameterizing non-bonded interactions is Bayesian parameter sampling, which encodes prior knowledge about the distribution through a prior distribution and fitting "constraints" through a likelihood (or energy) function, combining these together into the posterior distribution.  This posterior distribution gives us the dependence of the data fit to the parameters of the model, and allows us to make estimates of good parameter sets.

### MCMC Simulation

Part of the issue with Bayesian Inference is that, due to the likelihood functions that are chose, posterior distributions are often non-analytical and non-trivial to compute.  In this case, we turn to Markov Chain Monte Carlo sampling to estimate the posterior distribution.  By proposing random moves and accepting/rejecting them with a certain criteria, we draw samples from the posterior distribution. In the Bayesian paradigm (and with a careful choice of move proposals) this criteria turns out to be:

![metropolis](output/figures/metropolis.png)

Collecting enough of these samples gives a reasonable approximation of the posterior distribution.

## 2CLJQ model

### Force field

The 2 center Lennard-Jones + Quadrupole (2CLJQ) model is chosen for this investigation for several reasons:
1. It is a simple (4-parameter) model that focuses on non-bonded interactions (only bonded parameter is the bond length between LJ sites & the model is rigid).
2. Analytical correlations for physical properties (liquid density, saturation pressure, surface tension) are available from Stoll and Werth, so no actual simulation is required to evaluate these properties as a function of many parameter sets.
3.  A previous optimzation study for this set of fluids was done by Stobener, which provides a convenient benchmark for our Bayesian approach.

![2CLJQ](output/figures/2CLJQ.png)

### Experimental Data

Experimental data for the chosen compounds is available for our collaborators in the temperature range 55-95% of critical temperature.

### Uncertainty models

Because the computed properties come from analytical correlations rather than direct simulations, one must account for uncertainy inherent in the correlations (due to their inability to fit every temperature and property combination perfectly.  These uncertainty models are described in `uncertainty_models.py`

## Experiments
In this study, we attempt to reparameterize the 2CLJQ model by Bayesian parameter sampling from the analytical correlations provided by Stoll and Werth.

### Compounds and criteria chosen
Parameterization was performed in the following cases, to match the criteria of Stobener: a 2 criteria case (liquid density and saturation pressure) and a 3 criteria case (liquid density, saturation pressure, surface tension).

Compounds were chosen based on inclusion in the Stobener study and availability of experimental data.  
Included in this study: O2, N2, Br2, F2, C2H2, C2H4, C2H6, C2F4.
Note that C2F4 is not included in the 3-criteria case due to a lack of experimental surface tension data

### Simulation parameters

In all cases, simulations were run for 2,000,000 MCMC steps.  Running all 15 cases 4 at a time in parallel, this took roughly 3 hours.  In each case, 10 data points between 55-95% of critical temperature were chosen (unless data was not available for this entire range, in which case 10 data points were chosen from the available range).  The likelihood function was computed by forming an normal distribution around each experimental value with a shape parameter determined by the model uncertainty and then attempting to reproduce the experimental value with the analytical correlations and parameters.  The probability of this value being produced under the normal distribution is then multiplied by the corresponding probability for all other experimental values.  In this way, the closer the predictions are to the experimental values, the higher the likelihood.

Priors were chosen to be simple exponential priors, which encode that the values must be positive, and likely to be in a region based on the rate parameter.

### Results

As a metric of evaluating our parameter sets, we use average percent deviations from experimental values (i.e. for each property, compute the deviation from experiment for each temperature and average them), and compare these deviations to those from the Stobener parameter sets. 

In the following table, we include two parameter sets from our simulations:
1. (PM) A maximum a posterior (MAP) parameter set - essentially the parameter set that was sampled the most during the simulation.
2. (PP) A pareto set based on minimizing the *total* average percent deviation over all properties - i.e., it chooses the parameter set that minimizes the *combined* deviations. 
These are compared to:
3. (SP) A selected Pareto set from the Stobener paper.

For the 2-criteria case:

![2_criteria_devs](output/figures/results/2_criteria_devs.png)

For the 3 criteria case:

![3_criteria_devs](output/figures/results/3_criteria_devs.png)

In tabular form:

![2_criteria_table](output/figures/results/2_criteria_table.png)

![3_criteria_table](output/figures/results/3_criteria_table.png)

For the most part, the Parambayes MAP (PM) parameter sets end up doing worse than the Stobener Pareto (SP) set, but the Parambayes Pareto (PP) set does slightly better.  This indicates that there is influence from the prior that is shifting the posterior distributions away from the desired values.  This could be mitigated by using a more informed prior, perhaps by fitting a gamma or normal distribution to a shorter MCMC run.

### Individual Molecules/Cases

### O2 rhol+Psat
![O2_rhol+Psat_params](output/figures/results/O2_rhol+Psat_triangle_plot_params.png)

![O2_rhol+Psat_devs](output/figures/results/O2_rhol+Psat_triangle_plot_percent_dev_trace.png)

The PP param set is an improvement over the SP param set, while the PM set is slightly worse.

### O2 All
![O2_All_params](output/figures/results/O2_All_triangle_plot_params.png)

![O2_All_devs](output/figures/results/O2_All_triangle_plot_percent_dev_trace.png)


While both Parambayes sets are slightly worse than the SP, the Pareto set is comparable.


### N2 rhol+Psat
![N2_rhol+Psat_params](output/figures/results/N2_rhol+Psat_triangle_plot_params.png)

![N2_rhol+Psat_devs](output/figures/results/N2_rhol+Psat_triangle_plot_percent_dev_trace.png)


The PP param set is an improvement over the SP param set, while the PM set is slightly worse.

### N2 All
![N2_All_params](output/figures/results/N2_All_triangle_plot_params.png)

![N2_All_devs](output/figures/results/N2_All_triangle_plot_percent_dev_trace.png)

The PP param set is an improvement over the SP param set, while the PM set is slightly worse. However, all are pretty similar. Parambayes prioritizes lower deviations in the surface tension.


### Br2 rhol+Psat
![Br2_rhol+Psat_params](output/figures/results/Br2_rhol+Psat_triangle_plot_params.png)

![Br2_rhol+Psat_devs](output/figures/results/Br2_rhol+Psat_triangle_plot_percent_dev_trace.png)


The PP param set is an improvement over the SP param set, while the PM set is slightly worse.

### Br2 All
![Br2_All_params](output/figures/results/Br2_All_triangle_plot_params.png)

![Br2_All_devs](output/figures/results/Br2_All_triangle_plot_percent_dev_trace.png)


Both Parambayes sets are worse than the SP.  Parambayes prioritizes lower deviations in surface tension, to its detriment.

### F2 rhol+Psat
![F2_rhol+Psat_params](output/figures/results/F2_rhol+Psat_triangle_plot_params.png)

![F2_rhol+Psat_devs](output/figures/results/F2_rhol+Psat_triangle_plot_percent_dev_trace.png)


The PP param set is an improvement over the SP param set, while the PM set is slightly worse.

### F2 All
![F2_All_params](output/figures/results/F2_All_triangle_plot_params.png)

![F2_All_devs](output/figures/results/F2_All_triangle_plot_percent_dev_trace.png)

The PP param set is a significant improvement over the SP param set, whereas the PM set is significantly worse.

### C2H2 rhol+Psat
![C2H2_rhol+Psat_params](output/figures/results/C2H2_rhol+Psat_triangle_plot_params.png)

![C2H2_rhol+Psat_devs](output/figures/results/C2H2_rhol+Psat_triangle_plot_percent_dev_trace.png)


The PP param set is an improvement over the SP param set, while the PM set is slightly worse.

### C2H2 All
![C2H2_All_params](output/figures/results/C2H2_All_triangle_plot_params.png)

![C2H2_All_devs](output/figures/results/C2H2_All_triangle_plot_percent_dev_trace.png)

Both the PP and PM param sets are an improvement over the SP set.

### C2H4 rhol+Psat
![C2H4_rhol+Psat_params](output/figures/results/C2H4_rhol+Psat_triangle_plot_params.png)

![C2H4_rhol+Psat_devs](output/figures/results/C2H4_rhol+Psat_triangle_plot_percent_dev_trace.png)


The PP param set is an improvement over the SP param set, while the PM set essentially the same.

### C2H4 All
![C2H4_All_params](output/figures/results/C2H4_All_triangle_plot_params.png)

![C2H4_All_devs](output/figures/results/C2H4_All_triangle_plot_percent_dev_trace.png)


While both Parambayes sets are slightly worse than the SP, the Pareto set is comparable.

### C2H6 rhol+Psat
![C2H6_rhol+Psat_params](output/figures/results/C2H6_rhol+Psat_triangle_plot_params.png)

![C2H6_rhol+Psat_devs](output/figures/results/C2H6_rhol+Psat_triangle_plot_percent_dev_trace.png)


Both the PP and PM param sets are better than the SP set.

### C2H6 All
![C2H6_All_params](output/figures/results/C2H6_All_triangle_plot_params.png)

![C2H6_All_devs](output/figures/results/C2H6_All_triangle_plot_percent_dev_trace.png)

Both Parambayes sets are worse than the SP.  Parambayes prioritizes lower deviations in density.

### C2F4 rhol+Psat
![C2F4_rhol+Psat_params](output/figures/results/C2F4_rhol+Psat_triangle_plot_params.png)

![C2F4_rhol+Psat_devs](output/figures/results/C2F4_rhol+Psat_triangle_plot_percent_dev_trace.png)


The PP param set is an improvement over the SP param set, while the PM set is slightly worse.


### Advantages of Bayesian parameter sampling

In addition to providing us with parameter sets, we get additional information about the dependence of the model performance on the parameters. For example, in essentially all cases, we see that the parameters sigma, epsilon, and L are highly correlated.  This indicates that they should be co-optimized in the future.  We also see that in many cases, the Q (Quadrupole parameter) is driven towards zero, while in a few cases, we see that the quadrupole is definitively non-zero.  This indicates that the quadrupole parameter is not useful in some cases, but is essential in others.

## Reproducibility

In order to reproduce these MCMC simulations and outputs, after installing the package, run the script `gen_basic_runfiles.sh`.  It will produce the necessary runfiles for all the simulations done here, and then immediately run them.

### Copyright

Copyright (c) 2019, Owen Madin, Lenny Fobe, Ryan Morelock

#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
