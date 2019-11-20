"""
parambayes.py
Bayesian MCMC parameterization for CSCI 7000 (swe4s) class project

Handles the primary functions
"""
import numpy as np
from data_import import filter_thermo_data,import_literature_values,parse_data_ffs,calculate_uncertainties
from scipy.stats import distributions
from LennardJones_correlations import LennardJones
from LennardJones_2Center_correlations import LennardJones_2C

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
            
    def calc_posterior(self, prior, compound_2CLJ, chain_values):
        # def calc_posterior(model,eps,sig,L,Q,biasing_factor_UA=0,biasing_factor_AUA=0,biasing_factor_AUA_Q=0):

        dnorm = distributions.norm.logpdf

        logp = 0
        
        '''
        if chain_values[1] or chain_values[2] or chain_values[3] <= 0:
            #disallow values below 0 as nonphysical
            #print('Reject negative value')
            logp = -1*np.inf
        '''   
        
        logp += prior.sigma_prior_function.logpdf(chain_values[2], *prior.sigma_prior_values)
        logp += prior.epsilon_prior_function.logpdf(chain_values[1], *prior.epsilon_prior_values)
        # Create priors for parameters common to all models
        
        
        if chain_values[0] == 2:
            chain_values[4] = 0
            logp += self.biasing_factor[2]
            # Ensure Q=0 for UA model

        elif chain_values[0] == 0:
            chain_values[4] = 0
            logp += prior.L_prior_function.logpdf(chain_values[3], *prior.L_prior_values)
            logp += self.biasing_factor[0]
            # Add prior over L for AUA model

        elif chain_values[0] == 1:
            logp += prior.Q_prior_function.logpdf(chain_values[4], *prior.Q_prior_values)
            logp += prior.L_prior_function.logpdf(chain_values[3], *prior.L_prior_values)
            logp += self.biasing_factor[1]
            # Add priors for Q and L for AUA+Q model

        rhol_hat = rhol_hat_models(compound_2CLJ, self.thermo_data_rhoL[:, 0], *chain_values)  # [kg/m3]
        Psat_hat = Psat_hat_models(compound_2CLJ, self.thermo_data_Pv[:, 0], *chain_values)  # [kPa]
        SurfTens_hat = SurfTens_hat_models(compound_2CLJ, self.thermo_data_SurfTens[:, 0], *chain_values)
        # Compute properties at temperatures from experimental data

        # Data likelihood: Compute likelihood based on gaussian penalty function
        if self.properties == 'rhol':
            logp += sum(dnorm(self.thermo_data_rhoL[:, 1], rhol_hat, self.t_rhol**-2.))
            #logp += sum(dnorm(rhol_data,rhol_hat,t_rhol**-2.))
        elif self.properties == 'Psat':
            logp += sum(dnorm(self.thermo_data_Pv[:, 1], Psat_hat, self.t_Psat**-2.))
        elif self.properties == 'rhol+Psat':
            logp += sum(dnorm(self.thermo_data_rhoL[:, 1], rhol_hat, self.t_rhol**-2.))
            logp += sum(dnorm(self.thermo_data_Pv[:, 1], Psat_hat, self.t_Psat**-2.))
        elif self.properties == 'All':
            logp += sum(dnorm(self.thermo_data_rhoL[:, 1], rhol_hat, self.t_rhol**-2.))
            logp += sum(dnorm(self.thermo_data_Pv[:, 1], Psat_hat, self.t_Psat**-2.))
            logp += sum(dnorm(self.thermo_data_SurfTens[:, 1], SurfTens_hat, self.t_SurfTens**-2))

        return logp