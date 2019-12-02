#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 14:11:34 2019

@author: owenmadin
"""

import yaml
from datetime import date
import os

simulation_params = {}

#BASIC PARAMS
simulation_params['compound'] = 'C2H2'
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
#Options: exponential [loc (should always be 0), scale]
#         gamma [alpha,beta,loc,scale]
#See scipy.stats.expon and scipy.stats.gengamma

#RECORDING OPTIONS
simulation_params['save_traj'] = False
#Saves trajectories
simulation_params['label'] = 'test'
#Label for output files
today = str(date.today())

if os.path.exists('runfiles') is False:
    os.mkdir('runfiles')

filename = 'runfiles/'+simulation_params['compound'] + '_'+simulation_params['properties']+'_'+simulation_params['label']+'_'+today+'.yml' 


with open(filename,'w') as outfile:
    yaml.dump(simulation_params,outfile,default_flow_style=False)
