"""
parambayes.py
Bayesian MCMC parameterization for CSCI 7000 (swe4s) class project

Handles the primary functions
"""
import numpy as np
from data_import import filter_thermo_data,import_literature_values,parse_data_ffs,calculate_uncertainties
from scipy.stats import distributions
import math
import os 
from plotting import create_param_triangle_plot_4D,create_percent_dev_triangle_plot
from utility import rhol_hat_models,Psat_hat_models,SurfTens_hat_models,T_c_hat_models,computePercentDeviations
from datetime import date,datetime
import pickle
import matplotlib.pyplot as plt
from tqdm import tqdm
from shutil import rmtree

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
        """Return attributes of MCMC system
        """

        return {
            'compound': self.compound,
            'properties': self.properties,
            'T_range': self.T_range,
            'n_points': self.n_points,
            'steps': self.steps}
        
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
        
        logp += prior.sigma_prior_function.logpdf(chain_values[1], *prior.sigma_prior_values)
        logp += prior.epsilon_prior_function.logpdf(chain_values[0], *prior.epsilon_prior_values)
        logp += prior.Q_prior_function.logpdf(chain_values[3], *prior.Q_prior_values)
        logp += prior.L_prior_function.logpdf(chain_values[2], *prior.L_prior_values)
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
    
    def set_initial_state(self, prior, compound_2CLJ, initial_position=None):
        initial_logp = math.nan
        while math.isnan(initial_logp):
            initial_values = np.empty(4)

            rnorm = np.random.normal

            initial_values[0] = rnorm(self.ff_params_ref[0][0], self.ff_params_ref[0][0] / 10)
            initial_values[1] = rnorm(self.ff_params_ref[0][1], self.ff_params_ref[0][1] / 10)
            initial_values[2] = rnorm(self.ff_params_ref[0][2], self.ff_params_ref[0][2] / 10)
            initial_values[3] = rnorm(self.ff_params_ref[0][3], self.ff_params_ref[0][3] / 10)

            if initial_position is not None:
                initial_values = initial_position
            print('Markov Chain initialized at values:', initial_values)
            print('==============================')
            self.n_params = len(initial_values)
            self.prop_sd = np.asarray(initial_values) / 100
            initial_logp = self.calc_posterior(prior, compound_2CLJ, initial_values)
            if math.isnan(initial_logp):
                print('Nan detected! Finding new values')

        print('Initial log posterior:', initial_logp)
        print('==============================')
        self.initial_values = initial_values
        self.initial_logp = initial_logp
        self.initial_percent_deviation = computePercentDeviations(compound_2CLJ,
                                                                  self.thermo_data_rhoL[:, 0],
                                                                  self.thermo_data_Pv[:, 0],
                                                                  self.thermo_data_SurfTens[:, 0],
                                                                  self.initial_values,
                                                                  self.thermo_data_rhoL[:, 1],
                                                                  self.thermo_data_Pv[:, 1],
                                                                  self.thermo_data_SurfTens[:, 1],
                                                                  self.Tc_lit[0],
                                                                  rhol_hat_models,
                                                                  Psat_hat_models,
                                                                  SurfTens_hat_models,
                                                                  T_c_hat_models)
        self.new_lit_devs = []
        for i in range(len(self.lit_params[:,0])):
            self.new_lit_devs.append(computePercentDeviations(compound_2CLJ,
                                                                      self.thermo_data_rhoL[:, 0],
                                                                      self.thermo_data_Pv[:, 0],
                                                                      self.thermo_data_SurfTens[:, 0],
                                                                      self.lit_params[i],
                                                                      self.thermo_data_rhoL[:, 1],
                                                                      self.thermo_data_Pv[:, 1],
                                                                      self.thermo_data_SurfTens[:, 1],
                                                                      self.Tc_lit[0],
                                                                      rhol_hat_models,
                                                                      Psat_hat_models,
                                                                      SurfTens_hat_models,
                                                                      T_c_hat_models))
        self.new_lit_devs = np.asarray(self.new_lit_devs)
        
    def MCMC_Outerloop(self, prior, compound_2CLJ):
        self.trace = [self.initial_values]
        self.logp_trace = [self.initial_logp]
        self.percent_dev_trace = [self.initial_percent_deviation]

        print('Initializing Simulation...')
        print('Tuning Proposals...')
        print('==============================')
        for i in tqdm(range(self.steps)):
            if not i % 50000:
                #print('Iteration ' + str(i)), print('Log Posterior:', self.logp_trace[i])
                pass
            self.current_params = self.trace[i].copy()
            self.current_model = int(self.current_params[0])
            self.current_log_prob = self.logp_trace[i].copy()
            self.move_proposals = 0
            self.move_acceptances = 0

            new_params, new_log_prob, acceptance = self.MCMC_Steps(prior, compound_2CLJ)

            #self.move_proposals[int(self.current_params[0]), int(new_params[0])] += 1

            # accept_vector[i]=1
            self.logp_trace.append(new_log_prob)
            self.trace.append(new_params)
            self.percent_dev_trace.append(computePercentDeviations(compound_2CLJ,
                                                                   self.thermo_data_rhoL[:, 0],
                                                                   self.thermo_data_Pv[:, 0],
                                                                   self.thermo_data_SurfTens[:, 0],
                                                                   self.trace[i + 1],
                                                                   self.thermo_data_rhoL[:, 1],
                                                                   self.thermo_data_Pv[:, 1],
                                                                   self.thermo_data_SurfTens[:, 1],
                                                                   self.Tc_lit[0],
                                                                   rhol_hat_models,
                                                                   Psat_hat_models,
                                                                   SurfTens_hat_models,
                                                                   T_c_hat_models))

            if (not (i + 1) % self.tune_freq) and (i < self.tune_for):
                self.Tune_MCMC()
            if i == self.tune_for:
                self.move_proposals = 0
                self.move_acceptances = 0

                #print('Tuning complete!')
                #print('==============================')
        self.trace = np.asarray(self.trace)
        self.logp_trace = np.asarray(self.logp_trace)
        self.percent_dev_trace = np.asarray(self.percent_dev_trace)
        self.trace_tuned = self.trace[self.tune_for + 1:]
        self.logp_trace_tuned = self.logp_trace[self.tune_for + 1:]
        self.percent_dev_trace_tuned = self.percent_dev_trace[self.tune_for + 1:]
        print('Simulation Done!')
        print('==============================')
        
    def MCMC_Steps(self, prior, compound_2CLJ):
        proposed_params = self.current_params.copy()


        proposed_params, proposed_log_prob = self.parameter_proposal(prior, proposed_params, compound_2CLJ)
        self.move_proposals += 1
        alpha = (proposed_log_prob - self.current_log_prob)

        
        acceptance = self.accept_reject(alpha)
        if acceptance == 'True':
            new_log_prob = proposed_log_prob
            new_params = proposed_params
            self.move_acceptances += 1            

        elif acceptance == 'False':
            new_log_prob = self.current_log_prob
            new_params = self.current_params
   

        return new_params, new_log_prob, acceptance

    def accept_reject(self, alpha):
        urv = np.random.random()
        # Metropolis-Hastings accept/reject criteria
        if np.log(urv) < alpha:
            acceptance = 'True'
        else:
            acceptance = 'False'
        return acceptance

    
    def parameter_proposal(self, prior, proposed_params, compound_2CLJ):

        rnorm = np.random.normal
        modified_param = int(np.random.randint(0,4))

        proposed_params[modified_param] = rnorm(proposed_params[modified_param], self.prop_sd[modified_param])
        proposed_log_prob = self.calc_posterior(prior, compound_2CLJ, proposed_params)

        return proposed_params, proposed_log_prob

    def Tune_MCMC(self):
        # print(np.sum(self.move_proposals))
        acceptance_rate = np.sum(self.move_acceptances) / np.sum(self.move_proposals)
        # print(acceptance_rate)
        if acceptance_rate < 0.2:
            self.prop_sd *= 0.9
            # print('Yes')
        elif acceptance_rate > 0.5:
            self.prop_sd *= 1.1
            # print('No')
                
    def write_output(self, prior_dict, tag=None, save_traj=False):

	

        # Ask if output exists
        if os.path.isdir('../output') is False:
            os.mkdir('../output')
        if os.path.isdir('../output/' + self.compound) is False:
            os.mkdir('../output/' + self.compound)
        if os.path.isdir('../output/' + self.compound + '/' + self.properties) is False:
            os.mkdir('../output/' + self.compound + '/' + self.properties)

        path = '../output/' + self.compound + '/' + self.properties + '/' + self.compound + \
            '_' + self.properties + '_' + str(self.steps) + '_' + tag + '_' + str(date.today())

        if os.path.isdir(path):
            print('Directory Exists, overwriting')
            rmtree(path)
            os.mkdir(path)
        else:
            os.mkdir(path)

        os.mkdir(path + '/figures')

        print('Creating figures...')
        print('==============================')
        plt.plot(self.logp_trace_tuned)
        plt.savefig(path + '/figures/logp_trace.png')
        plt.close()


        create_param_triangle_plot_4D(
            self.trace_tuned,
            'triangle_plot_trace_model_0',
            self.lit_params,
            self.properties,
            self.compound,
            self.steps,
            file_loc=path +
            '/figures/')
        create_percent_dev_triangle_plot(
            self.percent_dev_trace_tuned,
            'triangle_plot_percent_dev_trace',
            self.new_lit_devs,
            self.properties,
            self.compound,
            self.steps,
            file_loc=path + '/figures/')

        print('Writing metadata...')
        print('==============================')
        self.write_datapoints(path)

        self.write_metadata(path, prior_dict)

        if save_traj:
            print('Saving Trajectories')
            print('==============================')
            self.write_traces(path)

    def write_datapoints(self, path):

        datapoints = {'Density Temperatures': self.thermo_data_rhoL[:, 0],
                      'Density Values': self.thermo_data_rhoL[:, 1],
                      'Density Measurement Uncertainties': self.thermo_data_rhoL[:, 2],
                      'Saturation Pressure Temperatures': self.thermo_data_Pv[:, 0],
                      'Saturation Pressure Values': self.thermo_data_Pv[:, 1],
                      'Saturation Pressure Uncertainties': self.thermo_data_Pv[:, 2],
                      'Surface Tension Temperatures': self.thermo_data_SurfTens[:, 0],
                      'Surface Tension Values': self.thermo_data_SurfTens[:, 1],
                      'Surface Tension Uncertainties': self.thermo_data_SurfTens[:, 2],
                      'Literature Critical Temperature': self.Tc_lit[0]}


        filename = path + '/datapoints.pkl'

        with open(filename, 'wb') as f:
            pickle.dump(datapoints,f)

    def write_metadata(self, path, prior_dict):

        metadata = self.get_attributes()
        filename = path + '/metadata.pkl'

        with open(filename, 'wb') as f:
            pickle.dump(metadata,f)



    def write_traces(self, path):
        if os.path.isdir(path + '/trace') == False:
            os.mkdir(path + '/trace')
        np.save(path + '/trace/trace.npy', self.trace_tuned)
        np.save(path + '/trace/logp_trace.npy', self.logp_trace_tuned)
        np.save(path + '/trace/percent_dev_trace_tuned.npy', self.percent_dev_trace_tuned)
        
class MCMC_Prior():
    """ Sets up a prior based on the user-specified prior types and parameters

    Attributes
    """

    def __init__(self, prior_dict):
        self.prior_dict = prior_dict

        self.dnorm = distributions.norm
        self.dgamma = distributions.gamma
        self.dgengamma = distributions.gengamma
        self.duni = distributions.uniform
        self.dlogit = distributions.logistic
        self.dexp = distributions.expon

    def epsilon_prior(self):
        eps_prior_type, eps_prior_vals = self.prior_dict['epsilon']

        if eps_prior_type == 'exponential':
            self.epsilon_prior_function = self.dexp
            self.epsilon_prior_values = [eps_prior_vals[0],eps_prior_vals[1]]
        elif eps_prior_type == 'gamma':
            self.epsilon_prior_function = self.dgamma
            self.epsilon_prior_values = [eps_prior_vals[0], eps_prior_vals[1],eps_prior_vals[2]]

    def sigma_prior(self):
        sig_prior_type, sig_prior_vals = self.prior_dict['sigma']

        if sig_prior_type == 'exponential':
            self.sigma_prior_function = self.dexp
            self.sigma_prior_values = [sig_prior_vals[0],sig_prior_vals[1]]
        elif sig_prior_type == 'gamma':
            self.sigma_prior_function = self.dgamma
            self.sigma_prior_values = [sig_prior_vals[0], sig_prior_vals[1],sig_prior_vals[2]]

    def L_prior(self):
        L_prior_type, L_prior_vals = self.prior_dict['L']

        if L_prior_type == 'exponential':
            self.L_prior_function = self.dexp
            self.L_prior_values = [L_prior_vals[0],L_prior_vals[1]]
        elif L_prior_type == 'gamma':
            self.L_prior_function = self.dgamma
            self.L_prior_values = [L_prior_vals[0], L_prior_vals[1], L_prior_vals[2]]

    def Q_prior(self):
        Q_prior_type, Q_prior_vals = self.prior_dict['Q']

        if Q_prior_type == 'exponential':
            self.Q_prior_function = self.dexp
            self.Q_prior_values = [Q_prior_vals[0], Q_prior_vals[1]]
        elif Q_prior_type == 'gamma':
            self.Q_prior_function = self.dgamma
            self.Q_prior_values = [Q_prior_vals[0], Q_prior_vals[1], Q_prior_vals[2]]


def main():
    pass
