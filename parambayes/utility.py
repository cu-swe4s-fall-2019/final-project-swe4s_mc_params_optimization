#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 21:44:21 2019

@author: owenmadin
"""
import numpy as np
from parambayes.LennardJones_2Center_correlations import LennardJones_2C

# Create functions that return properties for a given model, eps, sig


def rhol_hat_models(compound_2CLJ, Temp, eps, sig, L, Q):
    '''
    L_nm=L/10
    sig_nm=sig/10
    Q_nm=Q/10
    '''
    if not isinstance(compound_2CLJ, LennardJones_2C):
        raise TypeError("utility.py: compound_2CLJ must be a LennardJones_2Center_correlation.LennardJones2C object")
    if all([not isinstance(ins, (float, int, np.float, np.int, np.ndarray)) for ins in [Temp, eps, sig, L, Q]]):
        raise TypeError("utility.py: Temp, eps, sig, L, and Q must be an int/float/np.int/np.float/np.array!")
    rhol_hat = compound_2CLJ.rhol_hat_2CLJQ(Temp, eps, sig, L, Q)
    if all([isinstance(ins, (int, float, np.int, np.float)) for ins in [Temp, eps, sig, L, Q]]):
        raise TypeError("utility.py: At least one of Temp, eps, sig, L or Q must be a Numpy array")
    return rhol_hat  # [kg/m3]


def Psat_hat_models(compound_2CLJ, Temp, eps, sig, L, Q):

    Psat_hat = compound_2CLJ.Psat_hat_2CLJQ(Temp, eps, sig, L, Q)

    return Psat_hat  # [kPa]


def SurfTens_hat_models(compound_2CLJ, Temp, eps, sig, L, Q):
    '''
    L_nm=L/10
    sig_nm=sig/10
    Q_nm=Q/10
    '''

    SurfTens_hat = compound_2CLJ.ST_hat_2CLJQ(Temp, eps, sig, L, Q)

    return SurfTens_hat


def T_c_hat_models(compound_2CLJ, eps, sig, L, Q):
    '''
    L_nm=L/10
    sig_nm=sig/10
    Q_nm=Q/10
    '''

    T_c_hat = compound_2CLJ.T_c_hat_2CLJQ(eps, sig, L, Q)

    return T_c_hat


def computePercentDeviations(
        compound_2CLJ,
        temp_values_rhol,
        temp_values_psat,
        temp_values_surftens,
        parameter_values,
        rhol_data,
        psat_data,
        surftens_data,
        T_c_data,
        rhol_hat_models,
        Psat_hat_models,
        SurfTens_hat_models,
        T_c_hat_models):

    rhol_model = rhol_hat_models(compound_2CLJ, temp_values_rhol, *parameter_values)
    psat_model = Psat_hat_models(compound_2CLJ, temp_values_psat, *parameter_values)
    if len(surftens_data) != 0:
        surftens_model = SurfTens_hat_models(compound_2CLJ, temp_values_surftens, *parameter_values)
        surftens_deviation_vector = ((surftens_data - surftens_model) / surftens_data)**2
        surftens_mean_relative_deviation = np.sqrt(
            np.sum(surftens_deviation_vector) / np.size(surftens_deviation_vector)) * 100
    else:
        surftens_mean_relative_deviation = 0
    T_c_model = T_c_hat_models(compound_2CLJ, *parameter_values)

    rhol_deviation_vector = ((rhol_data - rhol_model) / rhol_data)**2
    psat_deviation_vector = ((psat_data - psat_model) / psat_data)**2

    T_c_relative_deviation = (T_c_data - T_c_model) * 100 / T_c_data

    rhol_mean_relative_deviation = np.sqrt(np.sum(rhol_deviation_vector) / np.size(rhol_deviation_vector)) * 100
    psat_mean_relative_deviation = np.sqrt(np.sum(psat_deviation_vector) / np.size(psat_deviation_vector)) * 100

    return rhol_mean_relative_deviation, psat_mean_relative_deviation, surftens_mean_relative_deviation, T_c_relative_deviation
