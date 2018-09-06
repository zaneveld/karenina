#/usr/bin/env python
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2016, The Karenina Project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "0.0.1-dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"

from optparse import OptionParser
from optparse import OptionGroup
from scipy.stats import norm
from numpy import array
from karenina.process import Process
from karenina.fit_timeseries import make_OU_objective_fn, fit_timeseries, aic
from os.path import join
import numpy as np
import os
import pandas as pd
import seaborn as sns
import sys
import matplotlib
import matplotlib.pyplot as plt


def make_option_parser():
    """Return an optparse OptionParser object"""

    parser = OptionParser(usage = "%prog -o ./simulation_results",
        description="This script fits microbiome change over time to Ornstein-Uhlenbeck (OU) models."+
                    "Demo fitting an OU model using default parameters.",
        version= __version__)

    required_options = OptionGroup(parser, "Required options")

    required_options.add_option('-o','--output', type="string",
                                help='the output directory for the simulation results')

    parser.add_option_group(required_options)

    optional_options = OptionGroup(parser, "Optional options")

    optional_options.add_option('--n_timepoints', type=int, default=300,
                                help='Benchmark the timeseries up to n_timepoints')

    optional_options.add_option('-v', '--verbose', action="store_true", dest="verbose", default=False,
                                help='-v, allows for verbose output' +
                                     ' [default: %default]')

    parser.add_option_group(optional_options)

    return parser


def write_logfile(output_dir,log_lines):

    f = open(output+"fit_timeseries_benchmark"+str(max_timepoints)+"_log.txt","w+")
    for line in log:
        f.write(str(line+"\n"))
    f.close()


def benchmark_simulated_datasets(max_timepoints = 300, output_dir = None, verbose = False,\
    simulation_type = "Ornstein-Uhlenbeck",\
    simulation_params={"lambda": 0.12, "delta": 0.25, "mu": 0.5}):
    """Test the ability of fit_timeseries to recover correct OU model params

    :param output: location for output log
    :param max_timepoints: maximum timepoints to test. The benchmark will test every number of timepoints up to this number. So specifying max_timpoints = 50 will test all numbers of timepoints between 2 and 50.
    :param simulation_type: string describing the model for the simulation (passed to a Process object)
    :param simulation_params: dict of parameters for the simulation (these will be passed on to the Process object).
    :param verbose: verbosity
    :return: output dataframe of benchmarked data
    """

    log = []   
    final_errors = {}
    dt = 0.1
   
    #model, n_timepoints, key, simulation_input,\
    #inferred_parameters, simulation_parameters,\
    # absolute_error, nLogLik, aic = ([] for i in range(9))
    
    results = []
    
    result_columns = ["model", "n_timepoints", "simulation_input", "inferred_parameters", "expected_parameters", "absolute_error", "AIC"]
    
    for curr_timepoint in list(range(1, max_timepoints+1)):
        
        
        curr_log,curr_results = benchmark_simulated_dataset(n_timepoints=curr_timepoint,\
          verbose = verbose,simulation_type=simulation_type,simulation_params=simulation_params,log=log,dt=dt)
        
        log.extend(curr_log)
        results.extend(curr_results)


    df = pd.DataFrame(results,columns = sorted(results[0].keys()))

    for line in log:
        if verbose:
            print(line)

    if output:
        #Write logfile
        log_filepath = output+"fit_timeseries_benchmark_log"+str(max_timepoints)+".txt"
        write_logfile(log_filepath,log)
        
        #Write benchamrk results file
        output_filepath = output+"fit_timeseries_benchmark_"+str(max_timepoints)+"_results.csv"
        df.to_csv(output_filepath, index=False)

        if verbose:
            print("Output log saved to: %s" %log_filepath)
            print("Benchmark results saved to: %s" %output_filepath)

    #TODO: is it essential that we return this reduced df instead of the full one?
    df = df[["model","n_timepoints","mag_err"]]
    return df



def generate_logfile_lines_from_result_dict(results):
    """Return a list of result lines

    :param results: a dict of results (all values must be able to be converted to strings)"""
    #Ensure results appear in the log in same order each time
    log = []
    for key in sorted(list(results.keys())):
        value = results[key]
        log.append("%s:%s" %(str(key),str(value)))
    
    log.append(str("*" * 80))
    return log


def benchmark_simulated_dataset(n_timepoints,verbose = False,\
    simulation_type = "Ornstein-Uhlenbeck",\
    simulation_params={"lambda": 0.12, "delta": 0.25, "mu": 0.5},local_optimizer=['L-BFGS-B'],\
    local_optimizer_niter=[5],local_optimizer_stepsize = 0.005,log=[],dt = 0.1):


    
    # run ou_process to get history
    ou = Process(start_coord=0.20, motion=simulation_type,
                 history=None, params= simulation_params)
    for t in range(0, n_timepoints):
        ou.update(dt)
    
    #Set timepoints and position (x) values
    xs = array(ou.History)
    ts = np.arange(0, len(ou.History)) * dt
    
    fn_to_optimize = make_OU_objective_fn(xs, ts)

    results = []
    
    # Estimate correct parameters for each tested local optimizer/niter value
    for local_optimizer in local_optimizer:
        for niter in local_optimizer_niter:

            #Set up simulation parameters 
            #(Using intentionally kinda bad starting estimates for benchmark)
            start_Sigma = 0.1
            start_Lambda = 0.0
            start_Theta = np.mean(xs)
            x0 = array([start_Sigma, start_Lambda, start_Theta])
            
            #set optimizer min/max bounds for each parameter
            xmax = array([1.0, 1.0, 1.0])
            xmin = array([0.0, 0.0, -1.0])


            #Results will be stored in the simulation_result dict
            result = {}
            result['benchmark_name'] = "_".join([local_optimizer,"niter",str(niter),"t",str(n_timepoints)])
            result['xs'] = xs
            result['ts'] = ts
            result['n_timepoints'] = n_timepoints
            result['local_optimizer']=local_optimizer
            result['local_optimizer_niter']=niter
            result['local_optimizer_param_names'] = ["Sigma","Lambda","Theta"]
            result['local_optimizer_start_param_values']=[start_Sigma,start_Lambda,start_Theta]
            result["expected_parameter_values"] = correct_values


            #Look up expected values from user input simulation parameters
            exp_Sigma,exp_Lambda,exp_Mu =\
               (simulation_parameters[param] for param in ['sigma','lambda','mu'])
            
            #Fit the model to the data
            inferred_params, nLogLik = \
              fit_timeseries(fn_to_optimize, x0, xmin, xmax,\
              stepsize= local_optimizer_stepsize,
              niter=niter, local_optimizer=local_optimizer)
 
            #Report inferred parameter values
            result["nLogLik"]=f_at_global_min
            result["AIC"]=aic(len(simulation_parameters),nLogLik)

            #Add results to the log
            log.extend(log_lines_from_result_dict(result))
            #NOTE: we do this now so later log entries with final errors
            #will appear last in the log
            
            #Calculate errors
            expected = list(correct_values)
            observed = list(inferred_params)
            parameter_names = result['local_optimizer_param_names']
            model_fit_error_results = calculate_errors(observed,expected,paramter_names)
            
            result.update(model_fit_error_results)
            #Add error results to the log
            log.extend(log_lines_from_result_dict(model_fit_error_results))

            #Update results with 
            results.append(result)

    return log,results

def calculate_errors(observed,expected,parameter_names=None):
    """Return per-parameter absolute and squared errors to results
    
    :param parameter names: list of strings describing parameter names (in order)
    :param observed: an array of observed values
    :param expected: an array of expected values     
    """
    
    if not parameter_names:
    
        #If we don't have informative parameter names, call them 'parameter_0', etc
        parameter_names = ["parameter_%i" %i for i in range(len(list(observed)))]
    
    #Calculate errors for individual parameters and all parameters 
    #together.
    
    result = {}

    #Calculate absolute and squared error for each
    #parameter individually
    absolute_errors = []
    squared_errors = []
    for i,curr_parameter in enumerate(parameter_names):
        absolute_error = abs(observed[i]-expected[i])
        absolute_errors.append(absolute_error)
        squared_error = absolute_error**2
        squared_errors.append(squared_errors)
        result['%s_absolute_error'%curr_parameter] = absolute_error
        result['%s_squared_error'%curr_parameter] = squared_error
    
    result['all_parameters_absolute_error']=absolute_errors
    result['all_parameters_squared_error']=squared_errors
    result['all_parameters_total_absolute_error'] = sum(absolute_errors)
    result['all_paramters_sum_of_squared_errors'] = sum(squared_errors)
       
    return result
            
 
def vis(df, output):
    """
    Visualizes benchmarking output error for various tested timepoints

    :param df: Dataframe containing model, timepoints, and list of errors
    :param output: output directory
    """

    # Generate visualizations for each model tested (Local optimizer)
    for model in df.model.unique():
        df_t = df.loc[df['model'] == model]
        df_split = pd.DataFrame(df_t.mag_err.values.tolist(),index=df_t.index)
        df_t = pd.concat([df_t.drop("mag_err", axis=1), df_split], axis=1)
        df_t.columns = ["model","n_timepoints","sigma_err","lambda_err","theta_err"]

        # Create individual dataframes for visualization
        df_sig = df_t[["n_timepoints","sigma_err"]]
        df_lam = df_t[["n_timepoints","lambda_err"]]
        df_the = df_t[["n_timepoints", "theta_err"]]

        # Create scatterplots with regression lines for each variable
        sig_plot = sns.lmplot(x="n_timepoints", y="sigma_err", data=df_sig, fit_reg=False, palette=["firebrick"], hue="sigma_err", legend=False)
        ax = plt.gca()
        ax.set_title("Observed Sigma Error")
        sig_plot.savefig(output+"benchmark_sigma_err.png")

        lam_plot = sns.lmplot(x="n_timepoints", y="lambda_err", data=df_lam, fit_reg=False, palette=["dodgerblue"], hue="lambda_err", legend=False)
        ax = plt.gca()
        ax.set_title("Observed Lambda Error")
        lam_plot.savefig(output + "benchmark_lambda_err.png")

        the_plot = sns.lmplot(x="n_timepoints", y="theta_err", data=df_the, fit_reg=False, palette=["goldenrod"], hue="theta_err", legend=False)
        ax = plt.gca()
        ax.set_title("Observed Theta Error")
        the_plot.savefig(output + "benchmark_theta_err.png")



           


def main():
    parser = make_option_parser()
    opts, args = parser.parse_args()
    verbose = opts.verbose
    n_timepoints = opts.n_timepoints
    output = opts.output

    if verbose:
        print(opts)

    if output is None:
        raise IOError('Output must be defined.')

    if not os.path.exists(output):
        os.makedirs(output)

    df = benchmark(n_timepoints, output, verbose)
    vis(df, output)

if __name__ == "__main__":
    main()
