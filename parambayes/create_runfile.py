#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 14:11:34 2019

@author: owenmadin
"""

import yaml
from datetime import date
import os
import argparse
import sys
import json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--compound",
                        help="Name of the compound to retrieve data for",
                        choices=["C2H2", "C2H4", "C2H6", "C2F4", "O2", "N2",
                                 "Br2", "F2"],
                        default="N2",
                        type=str
                        )
    parser.add_argument("--properties",
                        help="properties to sample from",
                        choices=["All", "rhol", "rhol+Psat"],
                        default="All",
                        type=str,
                        )
    parser.add_argument("--trange",
                        help="Reduced temperature range to pull samples from",
                        nargs=2,
                        type=float,
                        default=[0.55, 0.95]
                        )
    parser.add_argument("--n_steps",
                        help="Number of steps to run the MCMC simulation for",
                        default=1000000,
                        type=int,
                        )
    parser.add_argument("--n_data_points",
                        help="number of thermoproperty data points to use",
                        default=10,
                        type=int,
                        )
    parser.add_argument("--priors_JSON",
                        default='{"epsilon":["gamma",[134.665,0,0.2683]],' +
                                '"sigma":["gamma",[2646.618,0,0.0001246]],' +
                                '"L":["gamma",[53.1725,0,0.002059]],' +
                                '"Q":["exponential",[0,1]]}',
                        help="JSON formatted prior dictorionary." +
                             "Example of a JSON input is shown below:\n" +
                             '{"epsilon":["gamma",[134.665,0,0.2683]],' +
                                '"sigma":["gamma",[2646.618,0,0.0001246]],' +
                                '"L":["gamma",[53.1725,0,0.002059]],' +
                                '"Q":["exponential",[0,1]]}',
                        type=str,
                        )
    parser.add_argument("--date",
                        default="True",
                        choices=["True", "False"],
                        help="append date to the end of output label",
                        type=str
                        )
    parser.add_argument("--label",
                        default="prod",
                        help="file label for data output",
                        type=str,
                        )
    parser.add_argument("--save_traj",
                        default="False",
                        choices=["True", "False"],
                        help="save trajectory to an output file",
                        type=str,
                        )
    args = parser.parse_args()

    simulation_params = {}

    # BASIC PARAMS
    simulation_params["compound"] = args.compound
    # Compound to use (C2H2,C2H4,C2H6,C2F4,O2,N2,Br2,F2)
    simulation_params["properties"] = args.properties
    # Which properties to simulate ("rhol","rhol+Psat","All")

    if all([g >= 0 and g <= 100 for g in args.trange]):
        simulation_params["trange"] = args.trange
    else:
        print("create_runfile.py: trange provided is unphysical!",
              file=sys.stderr)
        sys.exit(1)
    # Temperature range (fraction of Tc)
    if args.n_steps > 0:
        simulation_params["steps"] = args.n_steps
    else:
        print("create_runfile.py: --n_steps must be a positive integer!",
              file=sys.stderr)
        sys.exit(1)
    # Number of MCMC steps
    if args.n_data_points > 0:
        simulation_params["number_data_points"] = args.n_data_points
    else:
        print("create_runfile.py: --n_data_points must be a positive integer!",
              file=sys.stderr)
        sys.exit(1)
<<<<<<< HEAD
    # [[72.21694778065564, 0, 0.5602885008739708],
    #  [1373.1569234248782, 0, 0.0002187094597278935],
    #  [24.898977020356103, 0, 0.0037408599155028944],
    #  [0.07855623224190525, 0, 0.8478830261261804]]

    # [[134.6652855956637, 0, 0.26832299241910723],
    #  [2646.618017963573, 0, 0.0001245996481485222],
    #  [53.172495977171614, 0, 0.002059236619566919],
    #  [1.8099969497105[134.6652855956637, 0, 0.26832299241910723],
    #  [484, 0, 0.06640261374979685],
    #  [134.6652855956637, 0, 0.26832299241910723]]
=======
>>>>>>> e17ba2aeadf2df3174a33d52eb7577e44dfee6ec

    # CUSTOM SIMULATION OPTIONS
    try:
        simulation_params["priors"] = json.loads(args.priors_JSON)
    except json.decoder.JSONDecodeError:
        print("create_runfile.py: Unable to load custom json input!",
              file=sys.stderr)
        sys.exit(1)

    for key in simulation_params["priors"].keys():
        if simulation_params["priors"][key][0] not in ["exponential", "gamma"]:
            print("create_runfile.py:", simulation_params["priors"][key][0],
                  "is an invalide prior distribution.",
                  "Please select from exponential or gamma.", file=sys.stderr)
            sys.exit(1)
        if simulation_params["priors"][key][0] == "exponential":
            if len(simulation_params["priors"][key][1]) != 2:
                print("create_runfile.py: exponential distributions must",
                      "have 2 inputs parameters", file=sys.stderr)
                sys.exit(1)
        if simulation_params["priors"][key][0] == "gamma":
            if len(simulation_params["priors"][key][1]) != 3:
                print("create_runfile.py: gamma distributions must",
                      "have 3 inputs parameters", file=sys.stderr)
                sys.exit(1)

    # Options: exponential [loc (should always be 0), scale]
    #         gamma [alpha,beta,loc,scale]
    # See scipy.stats.expon and scipy.stats.gengamma

    # RECORDING OPTIONS
    simulation_params["save_traj"] = bool(args.save_traj)
    # Saves trajectories
    simulation_params["label"] = args.label
    # Label for output files
    today = ""
    if bool(args.date):
        today = "_" + str(date.today())


    if os.path.exists("runfiles") is False:
        os.mkdir("runfiles")

    filename = "runfiles/" + simulation_params["compound"] + "_"
    filename += simulation_params["properties"] + "_"
    filename += simulation_params["label"]
    filename += today + ".yml"

    with open(filename, "w") as outfile:
        yaml.dump(simulation_params, outfile, default_flow_style=False)

    print("create_runfile.py: a run file has been written to " + filename)
    print("Please edit this file to your simulation specifications.")


if __name__ == "__main__":
    main()

