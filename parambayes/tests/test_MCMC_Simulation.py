#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 22:40:13 2019

@author: owenmadin
"""

import unittest
from parambayes import MCMC_Simulation,MCMC_Prior
from parambayes.LennardJones_2Center_correlations import LennardJones_2C
import numpy as np
import copy
import os
import math
import shutil



def params():
    simulation_params = {}

    # BASIC PARAMS
    simulation_params['compound'] = 'C2H4'
    # Compound to use (C2H2,C2H4,C2H6,C2F4,O2,N2,Br2,F2)
    simulation_params['properties'] = 'All'
    # Which properties to simulate ('rhol','rhol+Psat','All')
    simulation_params['trange'] = [0.55, 0.95]
    # Temperature range (fraction of Tc)
    simulation_params['steps'] = 200
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
        prior.make_priors()
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
        
class TestSetInitialState(unittest.TestCase):
    def setup():
        mcmc_simulator = MCMC_Simulation(simulation_params['compound'],
            simulation_params['trange'],
            simulation_params['properties'],
            simulation_params['number_data_points'],
            simulation_params['steps'])
        mcmc_simulator.prepare_data()
        compound_2CLJ = LennardJones_2C(mcmc_simulator.M_w)
        prior = MCMC_Prior(simulation_params['priors'])
        prior.make_priors()
        return mcmc_simulator,compound_2CLJ,prior
    def test_inputs(self):
        mcmc_simulator,compound_2CLJ,prior = TestSetInitialState.setup()
        self.assertRaises(TypeError,mcmc_simulator.set_initial_state,'text',compound_2CLJ)
        self.assertRaises(TypeError,mcmc_simulator.set_initial_state,prior,'text')
        self.assertRaises(TypeError,mcmc_simulator.set_initial_state,prior,compound_2CLJ,initial_position = 'text')
        self.assertRaises(IndexError,mcmc_simulator.set_initial_state,prior,compound_2CLJ,initial_position = [1,1,1])
        self.assertRaises(TypeError,mcmc_simulator.set_initial_state,prior,compound_2CLJ,initial_position = [1,1,1,'text'])
    def test_outputs_sanity(self):
        mcmc_simulator,compound_2CLJ,prior = TestSetInitialState.setup()
        mcmc_simulator.set_initial_state(prior,compound_2CLJ)
        self.assertEqual(mcmc_simulator.n_params,4)
        self.assertIs(np.float64,type(mcmc_simulator.initial_logp))
        self.assertIs(np.ndarray,type(mcmc_simulator.initial_values))
        self.assertIs(tuple,type(mcmc_simulator.initial_percent_deviation))
        self.assertIs(np.ndarray,type(mcmc_simulator.new_lit_devs))
    def test_logp_output(self):
        mcmc_simulator,compound_2CLJ,prior = TestSetInitialState.setup()
        self.assertRaises(ValueError,mcmc_simulator.set_initial_state,prior,compound_2CLJ,initial_position = [-1,-1,-1,-1])
class TestMCMCOuterLoop(unittest.TestCase):
    def setup():
        mcmc_simulator = MCMC_Simulation(simulation_params['compound'],
            simulation_params['trange'],
            simulation_params['properties'],
            simulation_params['number_data_points'],
            simulation_params['steps'],
            tune_for=100,
            tune_freq=50)
        mcmc_simulator.prepare_data()
        compound_2CLJ = LennardJones_2C(mcmc_simulator.M_w)
        prior = MCMC_Prior(simulation_params['priors'])
        prior.make_priors()
        mcmc_simulator.set_initial_state(prior,compound_2CLJ)
        return mcmc_simulator,compound_2CLJ,prior
    def test_inputs(self):
        mcmc_simulator,compound_2CLJ,prior = TestMCMCOuterLoop.setup()
        self.assertRaises(TypeError,mcmc_simulator.MCMC_Outerloop,prior,'text')
        self.assertRaises(TypeError,mcmc_simulator.MCMC_Outerloop,'text',compound_2CLJ)
    def test_outputs(self):
         mcmc_simulator,compound_2CLJ,prior = TestMCMCOuterLoop.setup()
         mcmc_simulator.MCMC_Outerloop(prior,compound_2CLJ)
         self.assertEqual(len(mcmc_simulator.logp_trace),mcmc_simulator.steps+1)
         self.assertEqual(len(mcmc_simulator.trace),mcmc_simulator.steps+1)
         self.assertEqual(len(mcmc_simulator.percent_dev_trace),mcmc_simulator.steps+1)
         self.assertEqual(len(mcmc_simulator.logp_trace_tuned),mcmc_simulator.steps-mcmc_simulator.tune_for)
         self.assertEqual(len(mcmc_simulator.trace_tuned),mcmc_simulator.steps-mcmc_simulator.tune_for)
         self.assertEqual(len(mcmc_simulator.percent_dev_trace_tuned),mcmc_simulator.steps-mcmc_simulator.tune_for)
         self.assertGreaterEqual(mcmc_simulator.move_proposals,mcmc_simulator.move_acceptances)
                  

class TestAcceptReject(unittest.TestCase):
    def test_inputs(self):
        mcmc_simulator = MCMC_Simulation(simulation_params['compound'],
            simulation_params['trange'],
            simulation_params['properties'],
            simulation_params['number_data_points'],
            simulation_params['steps'])
        self.assertRaises(TypeError,mcmc_simulator.accept_reject,'text')
        acceptance = mcmc_simulator.accept_reject(-1.05)
        self.assertIs(bool,type(acceptance))

class TestParameterProposal(unittest.TestCase):
    def setup():
        mcmc_simulator = MCMC_Simulation(simulation_params['compound'],
            simulation_params['trange'],
            simulation_params['properties'],
            simulation_params['number_data_points'],
            simulation_params['steps'])
        mcmc_simulator.prepare_data()
        compound_2CLJ = LennardJones_2C(mcmc_simulator.M_w)
        prior = MCMC_Prior(simulation_params['priors'])
        prior.make_priors()
        mcmc_simulator.set_initial_state(prior,compound_2CLJ)
        return mcmc_simulator,compound_2CLJ,prior
    def test_sanity_check(self):
        mcmc_simulator,compound_2CLJ,prior = TestParameterProposal.setup()
        params = [1,1,1,1]
        proposed_params = copy.deepcopy(params)
        new_params, proposed_log_prob = mcmc_simulator.parameter_proposal(prior,proposed_params,compound_2CLJ)
        self.assertNotEqual(new_params,params)

class TestMCMCSteps(unittest.TestCase):
    def setup():
        mcmc_simulator = MCMC_Simulation(simulation_params['compound'],
                                         simulation_params['trange'],
                                         simulation_params['properties'],
                                         simulation_params['number_data_points'],
                                         simulation_params['steps'])
        mcmc_simulator.prepare_data()
        compound_2CLJ = LennardJones_2C(mcmc_simulator.M_w)
        prior = MCMC_Prior(simulation_params['priors'])
        prior.make_priors()
        mcmc_simulator.set_initial_state(prior,compound_2CLJ)
        mcmc_simulator.trace = [mcmc_simulator.initial_values]
        mcmc_simulator.logp_trace = [mcmc_simulator.initial_logp]
        mcmc_simulator.percent_dev_trace = [mcmc_simulator.initial_percent_deviation]
        mcmc_simulator.move_proposals = 0
        mcmc_simulator.move_acceptances = 0
        mcmc_simulator.current_params = mcmc_simulator.trace[0].copy()
        mcmc_simulator.current_model = int(mcmc_simulator.current_params[0])
        mcmc_simulator.current_log_prob = mcmc_simulator.logp_trace[0].copy()
        return mcmc_simulator,compound_2CLJ,prior
    def test_inputs(self):
        mcmc_simulator,compound_2CLJ,prior = TestMCMCSteps.setup()

        self.assertRaises(TypeError,mcmc_simulator.MCMC_Steps,prior,'text')
        self.assertRaises(TypeError,mcmc_simulator.MCMC_Steps,'text',compound_2CLJ)
    def test_outputs(self):
        mcmc_simulator,compound_2CLJ,prior = TestMCMCSteps.setup()
        for i in range(100):
            move_proposals = copy.deepcopy(mcmc_simulator.move_proposals)
            move_acceptances = copy.deepcopy(mcmc_simulator.move_acceptances)
            new_params,new_logp,acceptance = mcmc_simulator.MCMC_Steps(prior,compound_2CLJ)
            self.assertEqual(move_proposals+1,mcmc_simulator.move_proposals)
            self.assertIs(bool,type(acceptance))
            
            if acceptance is True:
                self.assertEqual(move_acceptances+1,mcmc_simulator.move_acceptances)
                self.assertNotEqual(sum(new_params),sum(mcmc_simulator.current_params))
                self.assertNotEqual(new_logp,mcmc_simulator.current_log_prob)
                
            if acceptance is False:
                self.assertEqual(move_acceptances,mcmc_simulator.move_acceptances)
                self.assertEqual(new_params.all(),mcmc_simulator.current_params.all())
                self.assertEqual(new_logp,mcmc_simulator.current_log_prob)
                        
class TestTuneRJMC(unittest.TestCase):     
    def setup():
        mcmc_simulator = MCMC_Simulation(simulation_params['compound'],
                                         simulation_params['trange'],
                                         simulation_params['properties'],
                                         simulation_params['number_data_points'],
                                         simulation_params['steps'])
        mcmc_simulator.prepare_data()
        compound_2CLJ = LennardJones_2C(mcmc_simulator.M_w)
        prior = MCMC_Prior(simulation_params['priors'])
        prior.make_priors()
        mcmc_simulator.set_initial_state(prior,compound_2CLJ)
        mcmc_simulator.move_proposals = 0
        mcmc_simulator.move_acceptances = 0
        return mcmc_simulator,compound_2CLJ,prior

    def test_tune_down(self):
        mcmc_simulator,compound_2CLJ,prior = TestTuneRJMC.setup()
        mcmc_simulator.move_proposals = 10000
        mcmc_simulator.move_acceptances = 9000
        prop_sd_orig = copy.deepcopy(mcmc_simulator.prop_sd)
        mcmc_simulator.Tune_MCMC()
        self.assertGreater(sum(mcmc_simulator.prop_sd),sum(prop_sd_orig))
    def test_tune_up(self):
        mcmc_simulator,compound_2CLJ,prior = TestTuneRJMC.setup()
        mcmc_simulator.move_proposals = 10000
        mcmc_simulator.move_acceptances = 1000
        prop_sd_orig = copy.deepcopy(mcmc_simulator.prop_sd)
        mcmc_simulator.Tune_MCMC()
        self.assertLess(sum(mcmc_simulator.prop_sd),sum(prop_sd_orig))
    def test_errors(self):
        mcmc_simulator,compound_2CLJ,prior = TestTuneRJMC.setup()
        mcmc_simulator.move_proposals = 1000
        mcmc_simulator.move_acceptances = 1001
        self.assertRaises(ValueError,mcmc_simulator.Tune_MCMC)
        
        mcmc_simulator.move_proposals = -1000
        mcmc_simulator.move_acceptances = -1002
        self.assertRaises(ValueError,mcmc_simulator.Tune_MCMC)
        
        mcmc_simulator.move_proposals = 1000
        mcmc_simulator.move_acceptances = -100
        self.assertRaises(ValueError,mcmc_simulator.Tune_MCMC)
    
class TestWriteOutput(unittest.TestCase):
    def setup():
        mcmc_simulator = MCMC_Simulation(simulation_params['compound'],
                                         simulation_params['trange'],
                                         simulation_params['properties'],
                                         simulation_params['number_data_points'],
                                         simulation_params['steps'])
        mcmc_simulator.prepare_data()
        compound_2CLJ = LennardJones_2C(mcmc_simulator.M_w)
        prior = MCMC_Prior(simulation_params['priors'])
        prior.make_priors()
        mcmc_simulator.set_initial_state(prior,compound_2CLJ)
        mcmc_simulator.MCMC_Outerloop(prior,compound_2CLJ)  
        return mcmc_simulator,prior,compound_2CLJ
    
    def test_bad_input(self):
        mcmc_simulator,prior,compound_2CLJ = TestWriteOutput.setup()
        self.assertRaises(TypeError,mcmc_simulator.write_output,'wrong_prior')
        self.assertRaises(TypeError,mcmc_simulator.write_output,simulation_params['priors'],tag=1)
        self.assertRaises(TypeError,mcmc_simulator.write_output,simulation_params['priors'],save_traj='int')
    
    def test_output_created(self):
        mcmc_simulator,prior,compound_2CLJ = TestWriteOutput.setup()
        path = mcmc_simulator.write_output(simulation_params['priors'],output_path = True,save_traj = True)
        self.assertTrue(os.path.isdir(path))
        self.assertTrue(os.path.isdir(path+'/figures'))
        self.assertTrue(os.path.isfile(path+'/figures/logp_trace.png'))
        self.assertTrue(os.path.isfile(path+'/figures/triangle_plot_params.png'))
        self.assertTrue(os.path.isfile(path+'/figures/triangle_plot_percent_dev_trace.png'))
        self.assertTrue(os.path.isfile(path+'/datapoints.pkl'))
        self.assertTrue(os.path.isfile(path+'/metadata.pkl'))
        self.assertTrue(os.path.isdir(path+'/trace'))
        self.assertTrue(os.path.isfile(path+'/trace/trace.npy'))        
        self.assertTrue(os.path.isfile(path+'/trace/logp_trace.npy'))    
        self.assertTrue(os.path.isfile(path+'/trace/percent_dev_trace_tuned.npy'))
        shutil.rmtree(path)

class TestFindMaxima(unittest.TestCase):
    def setup():
        mcmc_simulator = MCMC_Simulation(simulation_params['compound'],
                                         simulation_params['trange'],
                                         simulation_params['properties'],
                                         simulation_params['number_data_points'],
                                         simulation_params['steps'])
        mcmc_simulator.prepare_data()
        compound_2CLJ = LennardJones_2C(mcmc_simulator.M_w)
        prior = MCMC_Prior(simulation_params['priors'])
        prior.make_priors()
        mcmc_simulator.set_initial_state(prior,compound_2CLJ)
        mcmc_simulator.MCMC_Outerloop(prior,compound_2CLJ)  
        return mcmc_simulator,prior,compound_2CLJ
    
    def test_input(self):
        mcmc_simulator,prior,compound_2CLJ = TestFindMaxima.setup()
        self.assertRaises(TypeError,mcmc_simulator.find_maxima,'not_array')
    def test_return_values(self):
        mcmc_simulator,prior,compound_2CLJ = TestFindMaxima.setup()
        mcmc_simulator.find_maxima(mcmc_simulator.trace)
        self.assertIsNotNone(mcmc_simulator.max_values)
        
        
    
    
