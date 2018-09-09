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

from matplotlib import use
use('Agg')
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
                                help='Benchmark the timeseries up to n_timepoints [default: %default]')

    optional_options.add_option('-v', '--verbose', action="store_true", dest="verbose", default=False,
                                help='-v, allows for verbose output' +
                                     ' [default: %default]')

    parser.add_option_group(optional_options)

    return parser


def write_logfile(output_fp,log_lines):
    """Write logfile lines to output filepath

    :param output_fp: a filepath (as string) for the output log file
    :param log_lines: lines of the logfile (linebreaks will be appended)
    """
    f = open(output_fp,"w+")
    for line in log_lines:
        f.write(str(line+"\n"))
    f.close()


def benchmark_simulated_datasets(max_timepoints = 300, output_dir = None,\
    simulation_type = "Ornstein-Uhlenbeck",\
    simulation_params={"lambda": 0.12, "delta": 0.25, "mu": 0.5},
    local_optimizer=['L-BFGS-B'],niter=[5],verbose=False):
    """Test the ability of fit_timeseries to recover correct OU model params

    :param output: location for output log
    :param max_timepoints: maximum timepoints to test. The benchmark will test every number of timepoints up to this number. So specifying max_timpoints = 50 will test all numbers of timepoints between 2 and 50.
    :param simulation_type: string describing the model for the simulation (passed to a Process object)
    :param simulation_params: dict of parameters for the simulation (these will be passed on to the Process object).
    :param local_optimizer: list of names of local optimizers to use
    :param niter: list of niter values to try. This parameter controls the number of independent restarts for global optimizer
    :param verbose: verbosity
    :return: output dataframe of benchmarked data
    """

    log = []   
    final_errors = {}
    dt = 1.0
   
    #model, n_timepoints, key, simulation_input,\
    #inferred_parameters, simulation_parameters,\
    # absolute_error, nLogLik, aic = ([] for i in range(9))
    
    results = []
    
    #result_columns = ["model", "n_timepoints", "simulation_input",\
    #  "inferred_parameters", "expected_parameters", "absolute_error", "AIC"]

    timepoints_to_simulate = list(range(1,max_timepoints+1))

    for curr_timepoint in timepoints_to_simulate:
        if curr_timepoint <= 2:
            continue    
        
        curr_log,curr_results = benchmark_simulated_dataset(n_timepoints=curr_timepoint,\
          verbose = verbose,simulation_type=simulation_type,simulation_params=simulation_params,log=log,dt=dt)
        
        log.extend(curr_log)
        results.extend(curr_results)


    df = pd.DataFrame(results,columns = sorted(results[0].keys()))

    for line in log:
        if verbose:
            print(line)

    if output_dir:
        #Write logfile
        log_filepath = join(output_dir,"fit_timeseries_benchmark_log_{0}.txt".format(max_timepoints))
        write_logfile(log_filepath,log)
        
        #Write benchamrk results file
        output_filepath = join(output_dir,"fit_timeseries_benchmark_{0}_results.csv".format(max_timepoints))
        df.to_csv(output_filepath, index=False)


        if verbose:
            print("Output log saved to: %s" %log_filepath)
            print("Benchmark results saved to: %s" %output_filepath)

    #TODO: is it essential that we return this reduced df instead of the full one?
    
    #df = df[["benchmark_name","n_timepoints","all_parameters_absolute_error"]]
    return df



def log_lines_from_result_dict(results):
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
                 history=[0.20,0.20], params= simulation_params)

    timepoints = list(range(0,n_timepoints))

    for t in timepoints:
        if verbose:
            print("simulating timepoint t:",t)
        
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
            result['model'] = "_".join([simulation_type,local_optimizer,"niter",str(niter)])
            result['benchmark_name'] = "_".join([local_optimizer,"niter",str(niter),"t",str(n_timepoints)])
            result['xs'] = xs
            result['ts'] = ts
            result['n_timepoints'] = n_timepoints
            result['local_optimizer']=local_optimizer
            result['local_optimizer_niter']=niter
            result['model_parameter_names'] = ["sigma","lambda","theta"]
            result['local_optimizer_start_param_values']=[start_Sigma,start_Lambda,start_Theta]
            result["expected_parameter_values"] = [simulation_params[param] for param in ['delta','lambda','mu']]

            #Look up expected values from user input simulation parameters
            exp_Sigma,exp_Lambda,exp_Mu =\
                result["expected_parameter_values"]
 
            #Fit the model to the data
            inferred_params, nLogLik = \
              fit_timeseries(fn_to_optimize, x0, xmin, xmax,\
              stepsize= local_optimizer_stepsize,
              niter=niter, local_optimizer=local_optimizer)
 
            for i,p in enumerate(result['model_parameter_names']):
                result["inferred_{}".format(p)] = inferred_params[i]
            
            #Report inferred parameter values
            result["nLogLik"]= nLogLik
            result["AIC"]=aic(len(simulation_params),nLogLik)

            #Add results to the log
            log.extend(log_lines_from_result_dict(result))
            #NOTE: we do this now so later log entries with final errors
            #will appear last in the log
            
            #Calculate errors
            expected = list(result['expected_parameter_values'])
            observed = list(inferred_params)
            parameter_names = result['model_parameter_names']
            model_fit_error_results = calculate_errors(observed,expected,parameter_names)
            
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
        squared_errors.append(squared_error)
        result['%s_absolute_error'%curr_parameter] = absolute_error
        result['%s_squared_error'%curr_parameter] = squared_error
   
    result['all_parameters_absolute_error']=absolute_errors
    result['all_parameters_squared_error']=squared_errors
    result['all_parameters_total_absolute_error'] = sum(absolute_errors)
    result['all_parameters_sum_of_squared_errors'] = sum(squared_errors)
       
    return result
            
 
def graph_absolute_errors(df, output_dir,model_name_col="model"):
    """
    Graph absolute error in simulated datasets with varying numbers of tested timepoints

    :param df: Dataframe containing model, timepoints, and list of errors
    :param output_dir: output directory
    :param model_name_col: column of the datafram that has the name of each model/benchmark
    """

    #TODO: we need to specify more clearly exactly what columns need to be in 
    #the dataframe 


    # Generate visualizations for each model tested (Local optimizer)
    for model in df[model_name_col].unique():
        df_t = df.loc[df[model_name_col] == model]
        #df_split = pd.DataFrame(df_t.mag_err.values.tolist(),index=df_t.index)
        #df_t = pd.concat([df_t.drop("mag_err", axis=1), df_split], axis=1)
        #df_t.columns = ["model","n_timepoints","sigma_err","lambda_err","theta_err"]
        curr_df = df_t
        #Assign each parameter a seaborn palette name so graphs will be different colors
        palettes_by_parameter = {"sigma":"firebrick","lambda":"dodgerblue","theta":"goldenrod"}
        
        model_parameters = df_t['model_parameter_names'] 

        for model_param in ['sigma','lambda','theta']:
            #TODO: is the below line necessary? Can't we just pass seaborn the full df?
            #curr_df = df_t[["n_timepoints",model_param]]  
            # Create and save scatterplots with regression lines for each variable
            curr_inferred_parameter_column = "inferred_{}".format(model_param)
            curr_absolute_error_column = "{}_absolute_error".format(model_param)
            curr_palette = palettes_by_parameter[model_param]
            curr_title = "Obseved Error {}".format(model_param.capitalize()) 

            curr_outfile = join(output_dir,"benchmark_{}_error.png".format(model_param))

            make_scatterplot(x_col='n_timepoints',y_col=curr_inferred_parameter_column,\
              data=curr_df,hue = curr_absolute_error_column,palette= [curr_palette],\
              outfile=curr_outfile,title = curr_title)

def make_scatterplot(x_col,y_col,data,hue,palette=["firebrick"],title='',\
        outfile='./results_scatter.png'):
        
        sig_plot = sns.lmplot(x=x_col, y=y_col, data=data,\
          fit_reg=False, palette=palette, hue=hue, legend=False)
        
        ax = plt.gca()
        ax.set_title(title)
        sig_plot.savefig(outfile)
        plt.close()

           


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

    df = benchmark_simulated_datasets(n_timepoints, output, verbose=verbose)
    graph_absolute_errors(df, output)

if __name__ == "__main__":
    main()
