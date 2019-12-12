#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 14:09:21 2019

@author: owenmadin
"""

import yaml
from parambayes.LennardJones_2Center_correlations import LennardJones_2C
from datetime import date
from parambayes import MCMC_Simulation, MCMC_Prior
import argparse


def parse_input_yaml(filepath):
    print('Loading simulation params from ' + filepath + '...')
    with open(filepath) as yfile:
        simulation_params = yaml.load(yfile)  # ,Loader=yaml.FullLoader)
    return simulation_params


def basic(simulation_params):

    print(simulation_params['priors'])
    mcmc_prior = MCMC_Prior(simulation_params['priors'])
    mcmc_prior.make_priors()

    mcmc_simulator = MCMC_Simulation(simulation_params['compound'],
                                     simulation_params['trange'],
                                     simulation_params['properties'],
                                     simulation_params['number_data_points'],
                                     simulation_params['steps'])

    mcmc_simulator.prepare_data()

    print('Simulation Attributes:', mcmc_simulator.get_attributes())

    compound_2CLJ = LennardJones_2C(mcmc_simulator.M_w)
    mcmc_simulator.set_initial_state(mcmc_prior, compound_2CLJ)

    mcmc_simulator.MCMC_Outerloop(mcmc_prior, compound_2CLJ)

    mcmc_simulator.write_output(
        simulation_params['priors'],
        tag=simulation_params['label'],
        save_traj=simulation_params['save_traj'])

    path = '../output/' + simulation_params['compound'] + '/' + \
        simulation_params['properties'] + '/' + \
        simulation_params['compound'] + \
        '_' + simulation_params['properties'] + '_' + \
        str(simulation_params['steps']) + '_' + simulation_params['label'] + \
        '_' + str(date.today()) + '/runfile.yaml'

    with open(path, 'w') as outfile:
        yaml.dump(simulation_params, outfile, default_flow_style=False)
    return mcmc_simulator


def main():
    parser = argparse.ArgumentParser(description='Find YAML file')

    parser.add_argument('--filepath', '-f',
                        type=str,
                        help='',
                        required=True)

    args = parser.parse_args()
    filepath = args.filepath
    print('Parsing simulation params')
    simulation_params = parse_input_yaml(filepath)
    mcmc_simulator = basic(simulation_params)

    print('Finished!')


if __name__ == '__main__':
    mcmc_simulator = main()
