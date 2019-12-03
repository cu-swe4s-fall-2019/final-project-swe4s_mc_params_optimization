#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 21:26:41 2019

@author: owenmadin
"""
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd
import yaml
from uncertainty_models import uncertainty_models


def filter_thermo_data(thermo_data, T_min, T_max, n_points):
    for name in thermo_data:
        df = thermo_data[name]

        df = df[df.values[:, 0] > T_min]
        df = df[df.values[:, 0] < T_max]
        if int(np.floor(df.shape[0] / (n_points - 1))) == 0:
            slicer = 1
        else:
            slicer = int(np.floor(df.shape[0] / (n_points - 1)))
        # print(slicer)
        df = df[::slicer]
        thermo_data[name] = df

    return thermo_data


def import_literature_values(criteria, compound):
    df = pd.read_csv(
        'data/Pareto_Hasse_' +
        criteria +
        '_criteria.txt',
        delimiter=' ',
        skiprows=2,
        usecols=[
            0,
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8])

    df = df[df.Substance == compound]
    df1 = df.iloc[:, 1:5]
    df2 = df.iloc[:, 5:9]
    df1 = df1[['epsilon', 'sigma', 'L', 'Q']]
    lit_values = np.asarray(df1)
    lit_values[:, 1:] /= 10

    return lit_values, np.asarray(df2)


def calculate_uncertainties(thermo_data, T_c):
    u_dict = {}
    for name in thermo_data:

        # Extract data from our data arrays
        data = np.asarray(thermo_data[name])
        T = data[:, 0]
        values = data[:, 1]
        u_exp = data[:, 2]

        pu_corr = uncertainty_models(T, T_c, name)
        u_corr = pu_corr * values

        u_tot = np.sqrt(u_corr**2 + u_exp**2)
        u_dict[name] = u_tot
    return u_dict


def parse_data_ffs(compound):
    fname = "data/lit_forcefields/" + compound + ".yaml"
    with open(fname) as yfile:
        yfile = yaml.load(yfile)  # ,Loader=yaml.FullLoader)
    ff_params = []
    params = ['eps_lit', 'sig_lit', 'Lbond_lit', 'Q_lit']
    for name in params:
        ff_params.append(yfile["force_field_params"][name])

    ff_params_ref = np.transpose(np.asarray(ff_params))
    ff_params_ref[:, 1:] = ff_params_ref[:, 1:] / 10

    Tc_lit = np.loadtxt('data/TRC_data/' + compound + '/Tc.txt', skiprows=1)
    M_w = np.loadtxt('data/TRC_data/' + compound + '/Mw.txt', skiprows=1)

    df = pd.read_csv('data/NIST_bondlengths/NIST_bondlengths.txt',
                     delimiter='\t')
    df = df[df.Compound == compound]
    NIST_bondlength = np.asarray(df)

    data = ['rhoL', 'Pv', 'SurfTens']
    data_dict = {}
    for name in data:
        df = pd.read_csv('data/TRC_data/' + compound + '/' +
                         name + '.txt', sep='\t')
        df = df.dropna()
        data_dict[name] = df
    return ff_params_ref, Tc_lit, M_w, data_dict, NIST_bondlength[0][1] / 10
