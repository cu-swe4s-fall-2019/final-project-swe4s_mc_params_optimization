"""
parambayes.py
Bayesian MCMC parameterization for CSCI 7000 (swe4s) class project

Handles the primary functions
"""
import numpy as np

class MCMC_Simulation():
    """ Builds an object that runs an RJMC simulation based on the parameters the user gives to it


    Attributes:


    """

    def __init__(self, compound, T_range, properties, n_points, steps, 
                 tune_freq=100, tune_for=10000):
        """Initializes the basic state of the simulator object.

        Parameters
        ----------
        data_physical_state: dict
            Passes information about the physical state of the simulation, i.e.
            which compound, what temperature range, which properties and how
            many data points
        steps: int
            Number of steps which the simulation should run for.
        swap_freq: float
            Percentage of times the simulation tries to jump between models
        biasing_factor: array
            Applies a biasing factor to a certain model
        prior: class
            Initializing priors for RJMC sampling
        """

        self.compound = compound
        self.T_range = T_range
        self.properties = properties
        self.n_points = n_points
        self.steps = steps
        self.tune_for = tune_for
        self.tune_freq = tune_freq

    def get_attributes(self):
        """Return attributes of RJMC system
        """

        return {
            'compound': self.compound,
            'properties': self.properties,
            'T_range': self.T_range,
            'n_points': self.n_points,
            'steps': self.steps,
            'swap_freq': self.swap_freq,
            'biasing_factor': self.biasing_factor}
        
    def prepare_data(self):
        """From input parameters, pull appropriate experimental data and
        uncertainty information.
        """

        self.ff_params_ref, self.Tc_lit, self.M_w, thermo_data, self.NIST_bondlength = parse_data_ffs(self.compound)
        # Retrieve force field literature values, constants, and thermo data

        self.T_min = self.T_range[0] * self.Tc_lit[0]
        self.T_max = self.T_range[1] * self.Tc_lit[0]

        # Select temperature range of data points to select, and how many
        # temperatures within that range to use data at.

        thermo_data = filter_thermo_data(thermo_data, self.T_min, self.T_max, self.n_points)

        # Filter data to selected conditions.

        uncertainties = calculate_uncertainties(thermo_data, self.Tc_lit[0])
        # Calculate uncertainties for each data point, based on combination of
        # experimental uncertainty and correlation uncertainty

        self.thermo_data_rhoL = np.asarray(thermo_data['rhoL'])
        self.thermo_data_Pv = np.asarray(thermo_data['Pv'])
        self.thermo_data_SurfTens = np.asarray(thermo_data['SurfTens'])
        # Convert dictionaries to numpy arrays

        # RJMC stuff

        # Calculate the estimated standard deviation
        sd_rhol = uncertainties['rhoL'] / 2.
        sd_Psat = uncertainties['Pv'] / 2.
        sd_SurfTens = uncertainties['SurfTens'] / 2

        # Calculate the precision in each property
        self.t_rhol = np.sqrt(1. / sd_rhol)
        self.t_Psat = np.sqrt(1. / sd_Psat)
        self.t_SurfTens = np.sqrt(1. / sd_SurfTens)
        
        if self.properties == 'rhol+Psat':    
            self.lit_params, self.lit_devs = import_literature_values('two', self.compound)
        elif self.properties == 'All':
            self.lit_params, self.lit_devs = import_literature_values('three',self.compound)