#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 22:40:13 2019

@author: owenmadin
"""

import unittest
from parambayes import MCMC_Simulation

def params():
    simulation_params = {}

    #BASIC PARAMS
    simulation_params['compound'] = 'C2H4'
    # Compound to use (C2H2,C2H4,C2H6,C2F4,O2,N2,Br2,F2)
    simulation_params['properties'] = 'All'
    # Which properties to simulate ('rhol','rhol+Psat','All')
    simulation_params['trange'] = [0.55,0.95]
    #Temperature range (fraction of Tc)
    simulation_params['steps'] = 2000000
    #Number of MCMC steps
    simulation_params['number_data_points'] = 10
    
    
    #CUSTOM SIMULATION OPTIONS
    simulation_params['priors'] = {'epsilon': ['exponential', [0,400]],
            'sigma': ['exponential', [0,1]],
            'L': ['exponential', [0,1]],
            'Q': ['exponential', [0,1]]}
    return simulation_params

simulation_params = params()

class TestInit(unittest.TestCase):
    def test_AllInputsExist(self):
        self.assertRaises(ValueError,MCMC_Simulation,None,simulation_params['trange'],simulation_params['properties'],simulation_params['number_data_points'],simulation_params['steps'])
        self.assertRaises(ValueError,MCMC_Simulation,simulation_params['compound'],None,simulation_params['properties'],simulation_params['number_data_points'],simulation_params['steps'])
        self.assertRaises(ValueError,MCMC_Simulation,simulation_params['compound'],simulation_params['trange'],None,simulation_params['number_data_points'],simulation_params['steps'])
        self.assertRaises(ValueError,MCMC_Simulation,simulation_params['compound'],simulation_params['trange'],simulation_params['properties'],None,simulation_params['steps'])
        self.assertRaises(ValueError,MCMC_Simulation,simulation_params['compound'],simulation_params['trange'],simulation_params['properties'],simulation_params['number_data_points'],None)
    def test_correct_compounds(self):
        self.assertRaises(ValueError,MCMC_Simulation,'NH4',simulation_params['trange'],simulation_params['properties'],simulation_params['number_data_points'],simulation_params['steps'])
    def test_correct_properties(self):
        self.assertRaises(ValueError,MCMC_Simulation,simulation_params['compound'],simulation_params['trange'],'Enthalpy',simulation_params['number_data_points'],simulation_params['steps'])
        self.assertRaises(ValueError,MCMC_Simulation,simulation_params['compound'],simulation_params['trange'],'rhol+PSat',simulation_params['number_data_points'],simulation_params['steps'])
    def test_trange(self):
        self.assertRaises(TypeError,MCMC_Simulation,simulation_params['compound'],['int',[1,1,1]],simulation_params['properties'],simulation_params['number_data_points'],simulation_params['steps'])
        self.assertRaises(ValueError,MCMC_Simulation,simulation_params['compound'],[-1,0.5],simulation_params['properties'],simulation_params['number_data_points'],simulation_params['steps'])
        self.assertRaises(ValueError,MCMC_Simulation,simulation_params['compound'],[0.5,2],simulation_params['properties'],simulation_params['number_data_points'],simulation_params['steps'])
        self.assertRaises(ValueError,MCMC_Simulation,simulation_params['compound'],[0.6,0.5],simulation_params['properties'],simulation_params['number_data_points'],simulation_params['steps'])
        self.assertRaises(ValueError,MCMC_Simulation,simulation_params['compound'],[0.5,0.5],simulation_params['properties'],simulation_params['number_data_points'],simulation_params['steps'])
    def test_num_data_points(self):
        self.assertRaises(TypeError,MCMC_Simulation,simulation_params['compound'],simulation_params['trange'],simulation_params['properties'],'int',simulation_params['steps'])
        self.assertRaises(ValueError,MCMC_Simulation,simulation_params['compound'],simulation_params['trange'],simulation_params['properties'],-10,simulation_params['steps'])
    def test_steps(self):
        self.assertRaises(TypeError,MCMC_Simulation,simulation_params['compound'],simulation_params['trange'],simulation_params['properties'],simulation_params['number_data_points'],100.3)
        self.assertRaises(ValueError,MCMC_Simulation,simulation_params['compound'],simulation_params['trange'],simulation_params['properties'],simulation_params['number_data_points'],-10)
        
    def test_tune(self):
        self.assertRaises(TypeError,MCMC_Simulation,simulation_params['compound'],simulation_params['trange'],simulation_params['properties'],simulation_params['number_data_points'],simulation_params['steps'],tune_for='blop')
        self.assertRaises(ValueError,MCMC_Simulation,simulation_params['compound'],simulation_params['trange'],simulation_params['properties'],simulation_params['number_data_points'],simulation_params['steps'],tune_for=-3)

        self.assertRaises(TypeError,MCMC_Simulation,simulation_params['compound'],simulation_params['trange'],simulation_params['properties'],simulation_params['number_data_points'],simulation_params['steps'],tune_freq='blop')
        self.assertRaises(ValueError,MCMC_Simulation,simulation_params['compound'],simulation_params['trange'],simulation_params['properties'],simulation_params['number_data_points'],simulation_params['steps'],tune_freq=-3)