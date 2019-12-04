#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 18:58:08 2019

@author: owenmadin
"""

import emcee
from parambayes import MCMC_Simulation,MCMC_Prior
from LennardJones_2Center_correlations import LennardJones_2C
import numpy as np
import plotting


simulation_params = {}

#BASIC PARAMS
simulation_params['compound'] = 'N2'
# Compound to use (C2H2,C2H4,C2H6,C2F4,O2,N2,Br2,F2)
simulation_params['properties'] = 'rhol+Psat'
# Which properties to simulate ('rhol','rhol+Psat','All')
simulation_params['trange'] = [0.55,0.95]
#Temperature range (fraction of Tc)
simulation_params['steps'] = 2000000
#Number of MCMC steps
simulation_params['swap_freq'] = 0.1
#Frequency of model swaps
simulation_params['number_data_points'] = 10


#CUSTOM SIMULATION OPTIONS
simulation_params['priors'] = {'epsilon': ['gamma', [134.665,0,0.2683]],
        'sigma': ['gamma', [2646.618, 0, 0.0001246]],
        'L': ['gamma', [53.1725, 0, 0.002059]],
        'Q': ['exponential', [0,1]]}

prior = MCMC_Prior(simulation_params['priors'])
prior.epsilon_prior()
prior.sigma_prior()
prior.L_prior()
prior.Q_prior()

mcmc_simulator = MCMC_Simulation(simulation_params['compound'], 
                                 simulation_params['trange'],
                                 simulation_params['properties'],
                                 simulation_params['number_data_points'],
                                 simulation_params['steps'])

mcmc_simulator.prepare_data()
logp = mcmc_simulator.calc_posterior

compound_2CLJ = LennardJones_2C(mcmc_simulator.M_w)

n_walkers = 10

initial_pos = []
for i in range(n_walkers):    
    mcmc_simulator.set_initial_state(prior, compound_2CLJ)
    initial_pos.append(mcmc_simulator.initial_values)
initial_pos = np.asarray(initial_pos)

sampler = emcee.EnsembleSampler(n_walkers,4,mcmc_simulator.calc_posterior_emcee,args=[compound_2CLJ,prior])
sampler.run_mcmc(initial_pos,10000,skip_initial_state_check=True)
A = sampler.flatchain
plotting.create_param_triangle_plot_4D(A,'Test',mcmc_simulator.lit_params,simulation_params['properties'],simulation_params['compound'],100000)

