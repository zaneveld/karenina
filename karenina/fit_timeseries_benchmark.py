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


def benchmark(max_tp = 300, output = None, verbose = False):
    """
    Verifies that fit_timeseries recovers OU model params
    :param output: location for output log
    :param max_tp: maximum timepoints to test
    :param verbose: verbosity
    :return output dataframe of benchmarked data
    """
    if output:
        f = open(output+"fit_timeseries_benchmark"+str(max_tp)+"_log.txt","w+")
    log = []

    # generate several normal distributions
    test_normal_data = {}
    n_obs = 1000
    dt = 0.01
    for delta in [float(i) / 100.0 for i in range(0, 100, 1)]:
        curr_data = norm.rvs(loc=0, size=n_obs, scale=delta ** (2 * dt))
        test_normal_data[delta ** (2 * dt)] = curr_data
    BasicNormalData = test_normal_data

    # generate OU process for testing
    ou_process = Process(start_coord=0.20, motion="Ornstein-Uhlenbeck",
                         history=None, params={"lambda": 0.20, "delta": 0.25, "mu": 0.0})
    # run ou_process to get history
    for t in range(1, 30):
        dt = 1
        ou_process.update(dt)
    OU = ou_process

    final_errors = {}
    dt = 1
    model, n_tp, key, sim_i, sim_o, exp_o, err_o, nLogLik, aic_o = ([] for i in range(9))
    columns = ["model", "n_timepoints", "sim_input", "sim_output", "exp_output", "error", "aic"]
    for n_timepoints in list(range(1, max_tp+1)):
        log.append(str("Building OU model for %i timepoints" % n_timepoints))
        n_tp.append(n_timepoints)
        # run ou_process to get history
        ou = Process(start_coord=0.20, motion="Ornstein-Uhlenbeck",
                     history=None, params= {"lambda": 0.12, "delta": 0.25, "mu": 0.5})
        for t in range(0, n_timepoints):
            ou.update(dt)
        log.append(str(str(n_timepoints)+", "+str(ou.History)))
        xs = array(ou.History)
        ts = np.arange(0, len(ou.History)) * dt
        log.append(str(str(xs)+", "+str(ts)+", "+str(dt)))
        fn_to_optimize = make_OU_objective_fn(xs, ts)

        # Estimate correct parameters
        for niter in [5]:
            for local_optimizer in ['L-BFGS-B']:
                log.append(str("Running optimizer: "+str(local_optimizer)))
                model.append(local_optimizer)
                # Using intentionally kinda bad estimates
                start_Sigma = 0.1
                start_Lambda = 0.0
                start_Theta = np.mean(xs)
                log.append(str("niter: "+str(niter)))
                log.append(str("start_Theta: "+str(start_Theta)))
                key.append(["sigma","lambda","theta"])
                sim_i.append([start_Sigma,start_Lambda,start_Theta])
                log.append(str("n_timepoints: "+str(n_timepoints)))
                xmax = array([1.0, 1.0, 1.0])
                xmin = array([0.0, 0.0, -1.0])
                x0 = array([start_Sigma, start_Lambda, start_Theta])

                global_min, f_at_global_min = \
                    fit_timeseries(fn_to_optimize, x0, xmin, xmax, stepsize=0.005,
                                   niter=niter, local_optimizer=local_optimizer)

                log.append("OU result:")
                Sigma, Lambda, Theta = global_min
                correct_values = array([0.25, 0.12, 0.5])
                exp_o.append([0.25, 0.12, 0.5])
                final_error = global_min - correct_values
                log.append(str("Global min: "+str(global_min)))
                sim_o.append([Sigma, Lambda, Theta])
                log.append(str("f at Global min: "+str(f_at_global_min)))
                nLogLik.append(f_at_global_min)
                # aic calulated with 2*n_params-2*LN(-1*nloglik)
                aic_t = aic(niter,f_at_global_min)
                log.append(str("aic: " +str(aic_t)))
                aic_o.append(aic_t)

                final_errors["%s_%i_%i" % (local_optimizer, niter, n_timepoints)] = final_error
                log.append(str("*" * 80))
                log.append(
                    str("%s error: %.4f,%.4f,%.4f" % (local_optimizer,final_error[0],final_error[1],final_error[2])))
                err_o.append([abs(final_error[0]),abs(final_error[1]),abs(final_error[2])])
                log.append(str("*" * 80))
                log.append("")

    df = pd.DataFrame({
        "model": model,
        "n_timepoints": n_tp,
        "key": key,
        "sim_input": sim_i,
        "sim_output": sim_o,
        "exp_output": exp_o,
        "mag_err": err_o,
        "nLogLik": nLogLik,
        "aic": aic_o},
        columns = ["model", "n_timepoints", "key", "sim_input", "sim_output", "exp_output", "mag_err", "nLogLik", "aic"])
    for opt, err in final_errors.items():
        log.append(str("%s error: %.4f,%.4f,%.4f" % (opt, err[0], err[1], err[2])))

    for line in log:
        if verbose:
            print(line)
        if output:
            f.write(str(line+"\n"))

    if output:
        if verbose:
            print("Output log saved to: "+output+"fit_timeseries_benchmark_log"+str(max_tp)+".txt")
            print("Output saved to: "+output+"fit_timeseries_benchmark"+str(max_tp)+".csv")
        df.to_csv(output+"fit_timeseries_benchmark"+str(max_tp)+".csv", index=False)
        f.close()

    df = df[["model","n_timepoints","mag_err"]]
    return df

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