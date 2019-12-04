#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 22:40:13 2019

@author: owenmadin
"""

import unittest
from parambayes import MCMC_Simulation,MCMC_Prior
from LennardJones_2Center_correlations import LennardJones_2C

import math



def params():
    simulation_params = {}

    # BASIC PARAMS
    simulation_params['compound'] = 'C2H4'
    # Compound to use (C2H2,C2H4,C2H6,C2F4,O2,N2,Br2,F2)
    simulation_params['properties'] = 'All'
    # Which properties to simulate ('rhol','rhol+Psat','All')
    simulation_params['trange'] = [0.55, 0.95]
    # Temperature range (fraction of Tc)
    simulation_params['steps'] = 2000000
    # Number of MCMC steps
    simulation_params['number_data_points'] = 10

    # CUSTOM SIMULATION OPTIONS
    simulation_params['priors'] = {'epsilon': ['exponential', [0, 400]],
                                   'sigma': ['exponential', [0, 1]],
                                   'L': ['exponential', [0, 1]],
                                   'Q': ['exponential', [0, 1]]}
    return simulation_params


simulation_params = params()


class TestInit(unittest.TestCase):
    def test_AllInputsExist(self):
        self.assertRaises(
            ValueError,
            MCMC_Simulation,
            None,
            simulation_params['trange'],
            simulation_params['properties'],
            simulation_params['number_data_points'],
            simulation_params['steps'])
        self.assertRaises(
            ValueError,
            MCMC_Simulation,
            simulation_params['compound'],
            None,
            simulation_params['properties'],
            simulation_params['number_data_points'],
            simulation_params['steps'])
        self.assertRaises(
            ValueError,
            MCMC_Simulation,
            simulation_params['compound'],
            simulation_params['trange'],
            None,
            simulation_params['number_data_points'],
            simulation_params['steps'])
        self.assertRaises(
            ValueError,
            MCMC_Simulation,
            simulation_params['compound'],
            simulation_params['trange'],
            simulation_params['properties'],
            None,
            simulation_params['steps'])
        self.assertRaises(
            ValueError,
            MCMC_Simulation,
            simulation_params['compound'],
            simulation_params['trange'],
            simulation_params['properties'],
            simulation_params['number_data_points'],
            None)

    def test_correct_compounds(self):
        self.assertRaises(
            ValueError,
            MCMC_Simulation,
            'NH4',
            simulation_params['trange'],
            simulation_params['properties'],
            simulation_params['number_data_points'],
            simulation_params['steps'])

    def test_correct_properties(self):
        self.assertRaises(
            ValueError,
            MCMC_Simulation,
            simulation_params['compound'],
            simulation_params['trange'],
            'Enthalpy',
            simulation_params['number_data_points'],
            simulation_params['steps'])
        self.assertRaises(
            ValueError,
            MCMC_Simulation,
            simulation_params['compound'],
            simulation_params['trange'],
            'rhol+PSat',
            simulation_params['number_data_points'],
            simulation_params['steps'])

    def test_trange(self):
        self.assertRaises(TypeError, MCMC_Simulation, simulation_params['compound'], ['int', [
                          1, 1, 1]], simulation_params['properties'], simulation_params['number_data_points'], simulation_params['steps'])
        self.assertRaises(ValueError,
                          MCMC_Simulation,
                          simulation_params['compound'],
                          [-1,
                           0.5],
                          simulation_params['properties'],
                          simulation_params['number_data_points'],
                          simulation_params['steps'])
        self.assertRaises(ValueError, MCMC_Simulation, simulation_params['compound'], [
                          0.5, 2], simulation_params['properties'], simulation_params['number_data_points'], simulation_params['steps'])
        self.assertRaises(ValueError,
                          MCMC_Simulation,
                          simulation_params['compound'],
                          [0.6,
                           0.5],
                          simulation_params['properties'],
                          simulation_params['number_data_points'],
                          simulation_params['steps'])
        self.assertRaises(ValueError,
                          MCMC_Simulation,
                          simulation_params['compound'],
                          [0.5,
                           0.5],
                          simulation_params['properties'],
                          simulation_params['number_data_points'],
                          simulation_params['steps'])

    def test_num_data_points(self):
        self.assertRaises(
            TypeError,
            MCMC_Simulation,
            simulation_params['compound'],
            simulation_params['trange'],
            simulation_params['properties'],
            'int',
            simulation_params['steps'])
        self.assertRaises(
            ValueError,
            MCMC_Simulation,
            simulation_params['compound'],
            simulation_params['trange'],
            simulation_params['properties'],
            -10,
            simulation_params['steps'])

    def test_steps(self):
        self.assertRaises(
            TypeError,
            MCMC_Simulation,
            simulation_params['compound'],
            simulation_params['trange'],
            simulation_params['properties'],
            simulation_params['number_data_points'],
            100.3)
        self.assertRaises(
            ValueError,
            MCMC_Simulation,
            simulation_params['compound'],
            simulation_params['trange'],
            simulation_params['properties'],
            simulation_params['number_data_points'],
            -10)

    def test_tune(self):
        self.assertRaises(
            TypeError,
            MCMC_Simulation,
            simulation_params['compound'],
            simulation_params['trange'],
            simulation_params['properties'],
            simulation_params['number_data_points'],
            simulation_params['steps'],
            tune_for='blop')
        self.assertRaises(
            ValueError,
            MCMC_Simulation,
            simulation_params['compound'],
            simulation_params['trange'],
            simulation_params['properties'],
            simulation_params['number_data_points'],
            simulation_params['steps'],
            tune_for=-3)

        self.assertRaises(
            TypeError,
            MCMC_Simulation,
            simulation_params['compound'],
            simulation_params['trange'],
            simulation_params['properties'],
            simulation_params['number_data_points'],
            simulation_params['steps'],
            tune_freq='blop')
        self.assertRaises(
            ValueError,
            MCMC_Simulation,
            simulation_params['compound'],
            simulation_params['trange'],
            simulation_params['properties'],
            simulation_params['number_data_points'],
            simulation_params['steps'],
            tune_freq=-3)
        
class TestGetAttributes(unittest.TestCase):
    def test_returns(self):
        mcmc = MCMC_Simulation(simulation_params['compound'],
            simulation_params['trange'],
            simulation_params['properties'],
            simulation_params['number_data_points'],
            simulation_params['steps'])
        self.assertIs(dict,type(mcmc.get_attributes()))
    def test_values_exists(self):
        mcmc = MCMC_Simulation(simulation_params['compound'],
            simulation_params['trange'],
            simulation_params['properties'],
            simulation_params['number_data_points'],
            simulation_params['steps'])
        mcmc_dict = mcmc.get_attributes()
        for i in mcmc_dict:
            self.assertIsNotNone(mcmc_dict[i])
            
class TestPrepareData(unittest.TestCase):
    pass

class TestCalcPosterior(unittest.TestCase):
    def setup():
        mcmc_simulator = MCMC_Simulation(simulation_params['compound'],
            simulation_params['trange'],
            simulation_params['properties'],
            simulation_params['number_data_points'],
            simulation_params['steps'])
        mcmc_simulator.prepare_data()
        compound_2CLJ = LennardJones_2C(mcmc_simulator.M_w)
        prior = MCMC_Prior(simulation_params['priors'])
        prior.epsilon_prior()
        prior.sigma_prior()
        prior.L_prior()
        prior.Q_prior()
        return mcmc_simulator,compound_2CLJ,prior
    
    def test_correct_inputs(self):
        mcmc_simulator,compound_2CLJ,prior=TestCalcPosterior.setup()
        chain_values = [70,0.3,0.1,0.3]
        self.assertRaises(TypeError,mcmc_simulator.calc_posterior,prior,compound_2CLJ,'string')
        self.assertRaises(TypeError,mcmc_simulator.calc_posterior,prior,compound_2CLJ,1.1)
        self.assertRaises(IndexError,mcmc_simulator.calc_posterior,prior,compound_2CLJ,chain_values[1:])
        self.assertRaises(TypeError,mcmc_simulator.calc_posterior,'prior',compound_2CLJ,chain_values)
        self.assertRaises(TypeError,mcmc_simulator.calc_posterior,prior,'C2H6',chain_values)
        chain_values = [70,0.3,0.1,'blep']
        self.assertRaises(TypeError,mcmc_simulator.calc_posterior,prior,compound_2CLJ,chain_values)
    def test_nan_handle(self):
        mcmc_simulator,compound_2CLJ,prior=TestCalcPosterior.setup()
        chain_values = [70,0.3,math.nan,0.1]
        self.assertWarns(UserWarning,mcmc_simulator.calc_posterior,prior,compound_2CLJ,chain_values)
    

    
    
