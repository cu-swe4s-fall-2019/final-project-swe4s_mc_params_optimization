#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 21:33:24 2019

@author: owenmadin
"""
import numpy as np


def uncertainty_models(T, T_c, thermo_property):
    Tr = T / T_c
    u = np.zeros(np.size(Tr))

    # Linear models for uncertainties in the 2CLJQ correlation we are using,
    # determined from Messerly analysis of figure from Stobener, Stoll, Werth

    # Starts at 0.3% for low values and ramps up to 1% for large values
    if thermo_property == 'rhoL':
        for i in range(np.size(Tr)):
            if Tr[i] < 0.9:
                u[i] = 0.3
            elif 0.9 <= Tr[i] <= 0.95:
                u[i] = 0.3 + (1 - 0.3) * (Tr[i] - 0.9) / (0.95 - 0.9)
            else:
                u[i] = 1.0

    # Starts at 20% for low values and ramps down to 2% for large values
    if thermo_property == 'Pv':
        for i in range(np.size(Tr)):
            if Tr[i] <= 0.55:
                u[i] = 20
            elif 0.55 <= Tr[i] <= 0.7:
                u[i] = 20 + (2 - 20) * (Tr[i] - 0.55) / (0.7 - 0.55)
            else:
                u[i] = 2.0

    # Starts at 4% for low values and ramps up to 12% for higher values
    if thermo_property == 'SurfTens':
        for i in range(np.size(Tr)):
            if Tr[i] <= 0.75:
                u[i] = 4
            elif 0.75 <= Tr[i] <= 0.95:
                u[i] = 4 + (12 - 4) * (Tr[i] - 0.75) / (0.95 - 0.75)
            else:
                u[i] = 12.0
    u /= 100
    return u
