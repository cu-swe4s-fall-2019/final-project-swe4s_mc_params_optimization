"""
parambayes.py
Bayesian MCMC parameterization for CSCI 7000 (swe4s) class project

Handles the primary functions
"""
import numpy as np
from parambayes.data_import import (filter_thermo_data,
                                    import_literature_values,
                                    parse_data_ffs,
                                    calculate_uncertainties)
from scipy.stats import distributions
import math
import os
from parambayes.plotting import (create_param_triangle_plot_4D,
                                 create_percent_dev_triangle_plot)
from parambayes.utility import (rhol_hat_models, Psat_hat_models,
                                SurfTens_hat_models, T_c_hat_models,
                                computePercentDeviations)
from parambayes.LennardJones_2Center_correlations import LennardJones_2C
from datetime import date, datetime
import pickle
import matplotlib.pyplot as plt
from tqdm import tqdm
from shutil import rmtree
import warnings


class MCMC_Simulation():
    """ Builds an object that runs an RJMC simulation
    based on the parameters the user gives to it


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
        if (compound is None or T_range is None or properties
                is None or n_points is None or steps is None):
            raise ValueError('MCMC_Simulation: A vital ' +
                             'simulation parameter is missing')
        if compound not in ['C2H6', 'C2H2', 'C2H4', 'C2F4',
                            'C2Cl4', 'O2', 'N2', 'Br2', 'F2']:
            raise ValueError("MCMC_Simulation: Compound is not available. \
                             Available compounds are: 'C2H6','C2H2','C2H4'," +
                             "'C2F4','C2Cl4','O2','N2','Br2','F2'")
        if properties not in ['rhol+Psat', 'All', 'rhol', 'Psat']:
            raise ValueError(
                "MCMC_Simulation: Properties Not Implemented. " +
                "Currently Available: 'rhol+Psat','All','rhol','Psat'")

        for T in T_range:
            if not isinstance(T, (float, int)):
                raise TypeError("MCMC_Simulation: T_range must be " +
                                "of type 'int' or 'float'")
            if T < 0 or T > 1:
                raise ValueError("MCMC_Simulation: T_range must " +
                                 "be between 0 and 1 (fraction of T_c)")
        if T_range[0] >= T_range[1]:
            raise ValueError("MCMC_Simulation: First number in " +
                             "T range must be less than second")
        if not isinstance(n_points, int):
            raise TypeError("MCMC_Simulation.__init__: " +
                            "n_points must be of type 'int'")
        if n_points <= 0:
            raise ValueError("MCMC_Simulation.__init__: " +
                             "n_points must be of positive integer")
        if not isinstance(steps, int):
            raise TypeError("MCMC_Simulation.__init__: " +
                            "steps must be of type 'int'")
        if steps <= 0:
            raise ValueError("MCMC_Simulation.__init__: " +
                             "steps must be of positive integer")
        if not isinstance(tune_freq, int):
            raise TypeError("MCMC_Simulation.__init__: " +
                            "steps must be of type 'int'")
        if tune_freq <= 0:
            raise ValueError("MCMC_Simulation.__init__: " +
                             "steps must be of positive integer")
        if not isinstance(tune_for, int):
            raise TypeError("MCMC_Simulation.__init__: " +
                            "steps must be of type 'int'")
        if tune_for <= 0:
            raise ValueError("MCMC_Simulation.__init__: " +
                             "steps must be of positive integer")

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

        self.ff_params_ref, self.Tc_lit, self.M_w, thermo_data, self.NIST_bondlength = parse_data_ffs(self.compound)  # noqa

        if not isinstance(self.ff_params_ref, np.ndarray):
            raise TypeError('MCMC_Simulation.prepare_data: expected ' +
                            'to receive ff_params_ref as numpy array')
        if np.shape(self.ff_params_ref)[1] != 4:
            raise ValueError('MCMC_Simulation.prepare_data: ' +
                             'expected ff_params_ref to have 4 columns')
        if not isinstance(self.Tc_lit, np.ndarray):
            raise TypeError('MCMC_Simulation.prepare_data: ' +
                            'expected to receive Tc_lit as list')
        if self.Tc_lit[0] <= 0:
            raise ValueError('MCMC_Simulation.prepare_data: ' +
                             'Tc_lit must be positive')
        if self.M_w <= 0:
            raise ValueError('MCMC_simulation.prepare_data: ' +
                             'M_w must be positive')
        if not isinstance(thermo_data, dict):
            raise TypeError('MCMC_simulation.prepare_data: ' +
                            'Expected to receive thermo_data as dict')
        if self.NIST_bondlength <= 0:
            raise ValueError('MCMC_simulation.prepare_data: ' +
                             'NIST_bondlength should be positive')
            # Retrieve force field literature values,
            # constants, and thermo data

        self.T_min = self.T_range[0] * self.Tc_lit[0]
        self.T_max = self.T_range[1] * self.Tc_lit[0]

        # Select temperature range of data points to select, and how many
        # temperatures within that range to use data at.

        thermo_data = filter_thermo_data(
            thermo_data, self.T_min, self.T_max, self.n_points)

        # Filter data to selected conditions.

        uncertainties = calculate_uncertainties(thermo_data, self.Tc_lit[0])
        # Calculate uncertainties for each data point, based on combination of
        # experimental uncertainty and correlation uncertainty
        self.thermo_data_rhoL = np.asarray(thermo_data['rhoL'])
        self.thermo_data_Pv = np.asarray(thermo_data['Pv'])
        self.thermo_data_SurfTens = np.asarray(thermo_data['SurfTens'])
        # Convert dictionaries to numpy arrays

        # Calculate the estimated standard deviation
        sd_rhol = uncertainties['rhoL'] / 2.
        sd_Psat = uncertainties['Pv'] / 2.
        sd_SurfTens = uncertainties['SurfTens'] / 2

        # Calculate the precision in each property
        self.t_rhol = np.sqrt(1. / sd_rhol)
        self.t_Psat = np.sqrt(1. / sd_Psat)
        self.t_SurfTens = np.sqrt(1. / sd_SurfTens)

        if self.properties == 'rhol+Psat':
            self.lit_params, self.lit_devs = import_literature_values(
                'two', self.compound)
        elif self.properties == 'All':
            self.lit_params, self.lit_devs = import_literature_values(
                'three', self.compound)
        else:
            print('Warning: no reference data available for ' +
                  'FF params; comparison is not possible')
            self.lit_params, self.lit_devs = (import_literature_values
                                              ('three', self.compound))

    def calc_posterior(self, mcmc_prior, compound_2CLJ, chain_values):
        # def calc_posterior(model,eps,sig,L,Q,biasing_factor_UA=0,
        # biasing_factor_AUA=0,biasing_factor_AUA_Q=0):

        if not isinstance(chain_values, (list, np.ndarray)):
            raise TypeError('MCMC_Simulation.calc_posterior: ' +
                            'chain values must be list or ndarray')
        if np.shape(chain_values)[0] != 4:
            raise IndexError('MCMC_Simulation.calc_posterior: ' +
                             'chain values must have length 4')
        for value in chain_values:
            if not isinstance(value, (int, float)):
                raise TypeError('MCMC_Simulation.calc_posterior: ' +
                                'chain values must all be floats or integers')
        if not isinstance(compound_2CLJ, LennardJones_2C):
            raise TypeError('MCMC_Simulation.calc_posterior: ' +
                            'compound_2CLJ must be instance of ' +
                            'LennardJones_2C class')
        if not isinstance(mcmc_prior, MCMC_Prior):
            raise TypeError('MCMC_Simulation.calc_posterior: ' +
                            'prior must be instance of MCMC_Prior class')

        dnorm = distributions.norm.logpdf

        logp = 0

        logp += (mcmc_prior.priors['sigma']['function'].
                 logpdf(chain_values[1],
                        *mcmc_prior.priors['sigma']['values']))
        logp += mcmc_prior.priors['epsilon']['function'].logpdf(
            chain_values[0], *mcmc_prior.priors['epsilon']['values'])
        logp += (mcmc_prior.priors['Q']['function'].
                 logpdf(chain_values[3],
                        *mcmc_prior.priors['Q']['values']))
        logp += (mcmc_prior.priors['L']['function'].
                 logpdf(chain_values[2],
                        *mcmc_prior.priors['L']['values']))
        # Add priors for Q and L for AUA+Q model

        rhol_hat = rhol_hat_models(compound_2CLJ,
                                   self.thermo_data_rhoL[:, 0],
                                   *chain_values)  # [kg/m3]
        Psat_hat = Psat_hat_models(compound_2CLJ,
                                   self.thermo_data_Pv[:, 0],
                                   *chain_values)  # [kPa]
        SurfTens_hat = SurfTens_hat_models(compound_2CLJ,
                                           self.thermo_data_SurfTens[:, 0],
                                           *chain_values)
        # Compute properties at temperatures from experimental data

        # Data likelihood: Compute likelihood
        # based on gaussian penalty function
        if self.properties == 'rhol':
            logp += sum(dnorm(self.thermo_data_rhoL[:, 1],
                              rhol_hat, self.t_rhol**-2.))
            # logp += sum(dnorm(rhol_data,rhol_hat,t_rhol**-2.))
        elif self.properties == 'Psat':
            logp += sum(dnorm(self.thermo_data_Pv[:, 1],
                              Psat_hat, self.t_Psat**-2.))
        elif self.properties == 'rhol+Psat':
            logp += sum(dnorm(self.thermo_data_rhoL[:, 1],
                              rhol_hat, self.t_rhol**-2.))
            logp += sum(dnorm(self.thermo_data_Pv[:, 1],
                              Psat_hat, self.t_Psat**-2.))
        elif self.properties == 'All':
            logp += sum(dnorm(self.thermo_data_rhoL[:, 1],
                              rhol_hat, self.t_rhol**-2.))
            logp += sum(dnorm(self.thermo_data_Pv[:, 1],
                              Psat_hat, self.t_Psat**-2.))
            logp += sum(dnorm(self.thermo_data_SurfTens[:, 1],
                              SurfTens_hat, self.t_SurfTens**-2))

        if np.isnan(logp):
            logp = -1 * math.inf
            warnings.warn('Warning: Proposed move returned ' +
                          'logp of NaN. Setting logp to -inf')
        return logp

    def set_initial_state(self, mcmc_prior,
                          compound_2CLJ, initial_position=None):
        if not isinstance(mcmc_prior, MCMC_Prior):
            raise TypeError('MCMC_Simulation.set_initial_state: ' +
                            'prior must be an instance of MCMC_Prior object')
        if not isinstance(compound_2CLJ, LennardJones_2C):
            raise TypeError(
                'MCMC_Simulation.set_initial_state: compound_2CLJ ' +
                'must be an instance of LennardJones_2C object')
        if initial_position is not None:
            if not isinstance(initial_position, (list, np.ndarray)):
                raise TypeError(
                    'MCMC_Simulation.set_initial_state: ' +
                    'User supplied initial position must be list or ndarray')
            if np.shape(initial_position)[0] != 4:
                raise IndexError('MCMC_Simulation.set_initial_state: user ' +
                                 'supplied initial_position must have len 4')
            for value in initial_position:
                if not isinstance(value, (float, int)):
                    raise TypeError(
                        'MCMC_Simulation.set_intial_state: user supplied ' +
                        'initial position must be list of floats or ints')

        initial_logp = math.nan
        while math.isnan(initial_logp):
            initial_values = np.empty(4)

            rnorm = np.random.normal

            initial_values[0] = rnorm(self.ff_params_ref[0][0],
                                      self.ff_params_ref[0][0] / 100)
            initial_values[1] = rnorm(self.ff_params_ref[0][1],
                                      self.ff_params_ref[0][1] / 100)
            initial_values[2] = rnorm(self.ff_params_ref[0][2],
                                      self.ff_params_ref[0][2] / 100)
            initial_values[3] = rnorm(self.ff_params_ref[0][3],
                                      self.ff_params_ref[0][3] / 100)

            if initial_position is not None:
                initial_values = initial_position
            print('Markov Chain initialized at values:', initial_values)
            print('==============================')
            self.n_params = len(initial_values)
            self.prop_sd = np.asarray(initial_values) / 100
            initial_logp = self.calc_posterior(
                mcmc_prior, compound_2CLJ, initial_values)
            if math.isnan(initial_logp):
                print('Nan detected! Finding new values')

        print('Initial log posterior:', initial_logp)
        print('==============================')
        self.initial_values = initial_values
        self.initial_logp = initial_logp
        if self.initial_logp == -1 * math.inf:
            raise ValueError(
                'MCMC_Simulation.set_initial_state: initial state has ' +
                '0 probability.  This is probably due to user-provided values')
        self.initial_percent_deviation = computePercentDeviations(
            compound_2CLJ,
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
        for i in range(len(self.lit_params[:, 0])):
            self.new_lit_devs.append(computePercentDeviations
                                     (compound_2CLJ,
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

        if not isinstance(prior, MCMC_Prior):
            raise TypeError('MCMC_Simulation.set_initial_state: ' +
                            'prior must be an instance of MCMC_Prior object')
        if not isinstance(compound_2CLJ, LennardJones_2C):
            raise TypeError(
                'MCMC_Simulation.set_initial_state: ' +
                'compound_2CLJ must be an instance of LennardJones_2C object')
        self.trace = [self.initial_values]
        self.logp_trace = [self.initial_logp]
        self.percent_dev_trace = [self.initial_percent_deviation]

        print('Initializing Simulation...')
        print('Tuning Proposals...')
        print('==============================')
        for i in tqdm(range(self.steps)):
            if not i % 50000:
                # print('Iteration ' + str(i)),
                # print('Log Posterior:', self.logp_trace[i])
                pass
            self.current_params = self.trace[i].copy()
            self.current_model = int(self.current_params[0])
            self.current_log_prob = self.logp_trace[i].copy()
            self.move_proposals = 0
            self.move_acceptances = 0

            new_params, new_log_prob, acceptance = \
                self.MCMC_Steps(prior, compound_2CLJ)

            # self.move_proposals[int(self.current_params[0]),
            # int(new_params[0])] += 1

            # accept_vector[i]=1
            self.logp_trace.append(new_log_prob)
            self.trace.append(new_params)
            self.percent_dev_trace.append(
                computePercentDeviations(compound_2CLJ,
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

                # print('Tuning complete!')
                # print('==============================')
        self.trace = np.asarray(self.trace)
        self.logp_trace = np.asarray(self.logp_trace)
        self.percent_dev_trace = np.asarray(self.percent_dev_trace)
        self.trace_tuned = self.trace[self.tune_for + 1:]
        self.logp_trace_tuned = self.logp_trace[self.tune_for + 1:]
        self.percent_dev_trace_tuned = (self.percent_dev_trace[
                                        self.tune_for + 1:])
        print('Simulation Done!')
        print('==============================')

    def MCMC_Steps(self, prior, compound_2CLJ):

        if not isinstance(prior, MCMC_Prior):
            raise TypeError('MCMC_Simulation.set_initial_state: ' +
                            'prior must be an instance of MCMC_Prior object')
        if not isinstance(compound_2CLJ, LennardJones_2C):
            raise TypeError(
                'MCMC_Simulation.set_initial_state: ' +
                'compound_2CLJ must be an instance of LennardJones_2C object')
        proposed_params = self.current_params.copy()

        proposed_params, proposed_log_prob = (self.
                                              parameter_proposal(
                                                  prior, proposed_params,
                                                  compound_2CLJ))
        self.move_proposals += 1
        alpha = (proposed_log_prob - self.current_log_prob)

        acceptance = self.accept_reject(alpha)

        if alpha > 0:
            assert acceptance is True

        if acceptance is True:
            new_log_prob = proposed_log_prob
            new_params = proposed_params
            self.move_acceptances += 1

        elif acceptance is False:
            new_log_prob = self.current_log_prob
            new_params = self.current_params

        return new_params, new_log_prob, acceptance

    def accept_reject(self, alpha):
        if not isinstance(alpha, (float, np.float)):
            raise TypeError('MCMC_Simulation.accept_reject: ' +
                            'alpha must be float')
        urv = np.random.random()
        # Metropolis-Hastings accept/reject criteria
        if np.log(urv) < alpha:
            acceptance = True
        else:
            acceptance = False
        return acceptance

    def parameter_proposal(self, prior, proposed_params, compound_2CLJ):

        rnorm = np.random.normal
        modified_param = int(np.random.randint(0, 4))

        proposed_params[modified_param] = (rnorm(proposed_params[
            modified_param],
            self.prop_sd[
            modified_param]))
        proposed_log_prob = (self.calc_posterior(
            prior, compound_2CLJ, proposed_params))

        return proposed_params, proposed_log_prob

    def Tune_MCMC(self):
        # print(np.sum(self.move_proposals))
        if self.move_acceptances > self.move_proposals:
            raise ValueError(
                'MCMC_Simulation.Tune_MCMC: ' +
                'Somehow more moves accepted then proposed. ' +
                'Counter is likely broken')
        if self.move_acceptances < 0:
            raise ValueError('MCMC_Simulation.Tune_MCMC: ' +
                             'move_acceptance < 0.  Counter is broken')
        if self.move_proposals < 0:
            raise ValueError('MCMC_Simulation.Tune_MCMC: ' +
                             'move_proposals < 0. Counter is broken')
        acceptance_rate = (np.sum(self.move_acceptances) /
                           np.sum(self.move_proposals))
        # print(acceptance_rate)
        if acceptance_rate < 0.2:
            self.prop_sd *= 0.9
            # print('Yes')
        elif acceptance_rate > 0.5:
            self.prop_sd *= 1.1
            # print('No')

    def write_output(self, prior_dict, tag='',
                     save_traj=False, output_path=False):

        if not isinstance(prior_dict, dict):
            raise TypeError(
                'MCMC_Simulation.write_output: ' +
                'prior dict must be dictionary of priors ' +
                '(generate with MCMC_Priors)')
        if tag is not None:
            if not isinstance(tag, str):
                raise TypeError('MCMC_Simulation.write_output: ' +
                                'tag must be None or str')
        if not isinstance(save_traj, bool):
            raise TypeError('MCMC_Simulation.write_output: ' +
                            'save_traj must be bool')

        # Ask if output exists
        if os.path.isdir('../output') is False:
            os.mkdir('../output')
        if os.path.isdir('../output/' + self.compound) is False:
            os.mkdir('../output/' + self.compound)
        if (os.path.isdir('../output/' + self.compound
                          + '/' + self.properties) is False):
            os.mkdir('../output/' + self.compound + '/' + self.properties)

        path = '../output/' + self.compound + '/' + \
            self.properties + '/' + self.compound + \
            '_' + self.properties + '_' + str(self.steps) + \
            '_' + tag + '_' + str(date.today())

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
            'triangle_plot_params',
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
        

        if self.max_values is not None:    
            np.savetxt(path+'/max_values.csv', (self.max_values,
                                                self.max_percent_dev,
                                                self.ff_params_ref[0],
                                                self.ff_percent_dev),
                                                delimiter=',')
        

        if save_traj is True:
            print('Saving Trajectories')
            print('==============================')
            self.write_traces(path)
        if output_path is True:
            return path

    def write_datapoints(self, path):

        datapoints = {'Density Temperatures': self.thermo_data_rhoL[:, 0],
                      'Density Values': self.thermo_data_rhoL[:, 1],
                      'Density Measurement Uncertainties': (self.
                                                            thermo_data_rhoL[
                                                                :, 2]),
                      'Saturation Pressure Temperatures': (self.
                                                           thermo_data_Pv[
                                                               :, 0]),
                      'Saturation Pressure Values': self.thermo_data_Pv[:, 1],
                      'Saturation Pressure Uncertainties': (self.
                                                            thermo_data_Pv[
                                                                :, 2]),
                      'Surface Tension Temperatures': (self.
                                                       thermo_data_SurfTens[
                                                           :, 0]),
                      'Surface Tension Values': (self.
                                                 thermo_data_SurfTens[:, 1]),
                      'Surface Tension Uncertainties': (self.
                                                        thermo_data_SurfTens[
                                                            :, 2]),
                      'Literature Critical Temperature': self.Tc_lit[0]}

        filename = path + '/datapoints.pkl'

        with open(filename, 'wb') as f:
            pickle.dump(datapoints, f)

    def write_metadata(self, path, prior_dict):

        metadata = self.get_attributes()
        filename = path + '/metadata.pkl'

        with open(filename, 'wb') as f:
            pickle.dump(metadata, f)

    def write_traces(self, path):
        if os.path.isdir(path + '/trace') is False:
            os.mkdir(path + '/trace')
        np.save(path + '/trace/trace.npy', self.trace_tuned)
        np.save(path + '/trace/logp_trace.npy', self.logp_trace_tuned)
        np.save(path + '/trace/percent_dev_trace_tuned.npy',
                self.percent_dev_trace_tuned)

    def find_maxima(self, trace, compound_2CLJ):
        if not isinstance(trace, (np.ndarray, list)):
            raise TypeError('MCMC_Simulation.find_maxima: ' +
                            'trace must be np.ndarray or list')
        num_bins = 20
        hist = np.histogramdd(trace, bins=num_bins, density=True)
        val = hist[0].max()
        for i in range(num_bins):
            for j in range(num_bins):
                for k in range(num_bins):
                    for l in range(num_bins):
                        if hist[0][i][j][k][l] == val:
                            key = [i, j, k, l]
                            break
        max_values = []
        for index in range(len(key)):
            low = hist[1][index][key[index]]
            high = hist[1][index][key[index]]
            max_values.append((low + high) / 2)
        self.max_values = max_values
        self.max_percent_dev = computePercentDeviations(
            compound_2CLJ,
            self.thermo_data_rhoL[:, 0],
            self.thermo_data_Pv[:, 0],
            self.thermo_data_SurfTens[:, 0],
            self.max_values,
            self.thermo_data_rhoL[:, 1],
            self.thermo_data_Pv[:, 1],
            self.thermo_data_SurfTens[:, 1],
            self.Tc_lit[0],
            rhol_hat_models,
            Psat_hat_models,
            SurfTens_hat_models,
            T_c_hat_models)
        self.ff_percent_dev = computePercentDeviations(
            compound_2CLJ,
            self.thermo_data_rhoL[:, 0],
            self.thermo_data_Pv[:, 0],
            self.thermo_data_SurfTens[:, 0],
            self.ff_params_ref[0],
            self.thermo_data_rhoL[:, 1],
            self.thermo_data_Pv[:, 1],
            self.thermo_data_SurfTens[:, 1],
            self.Tc_lit[0],
            rhol_hat_models,
            Psat_hat_models,
            SurfTens_hat_models,
            T_c_hat_models)
        
            


class MCMC_Prior():
    """Sets up a prior based on the user-specified prior types and parameters

    Attributes
    """

    def __init__(self, prior_dict):

        if not isinstance(prior_dict, dict):
            raise TypeError("parambayes.py:MCMC_Prior: " +
                            "prior_dict must be a dictionary!")
        for key in prior_dict.keys():
            if not isinstance(key, str):
                raise TypeError("parambayes.py:MCMC_Prior: " +
                                "prior_dict keys must be strings!")

        self.prior_dict = prior_dict

        self.dnorm = distributions.norm
        self.dgamma = distributions.gamma
        self.dgengamma = distributions.gengamma
        self.duni = distributions.uniform
        self.dlogit = distributions.logistic
        self.dexp = distributions.expon
        self.str_2_fxn_map = {
            "exponential": self.dexp,
            "gamma": self.dgamma,
            "gengamma": self.dgengamma,
            "uniform": self.duni,
            "logistic": self.dlogit,
            "normal": self.dnorm,
        }
        for method in [prior_dict[a][0] for a in self.prior_dict.keys()]:
            if method not in self.str_2_fxn_map.keys():
                raise KeyError("parambayes.py:MCMC_Prior: " +
                               method +
                               " is not implemented in MCMC_Prior.\n" +
                               "Please select from the following " +
                               "distributions:\n" +
                               ", ".join(self.str_2_fxn_map.keys()))

    def make_priors(self):
        self.priors = {}
        for prior in self.prior_dict.keys():
            self.priors[prior] = {"function": self.str_2_fxn_map[
                self.prior_dict[prior][0]],
                "values": self.prior_dict[prior][1]
            }


def main():
    pass
