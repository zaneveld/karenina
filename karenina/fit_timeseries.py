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
from scipy.optimize import basinhopping,brute,differential_evolution
from scipy.stats import norm
from numpy import diff,inf,all,array
from karenina.process import Process
import zipfile
import numpy as np
import pandas as pd
import collections
import os
import math
import sys

def make_option_parser():
    """Return an optparse OptionParser object"""

    parser = OptionParser(usage = "%prog -o ./simulation_results",
        description="This script fits microbiome change over time to Ornstein-Uhlenbeck (OU) models."+
                    "Demo fitting an OU model using default parameters.",
        version= __version__)

    required_options = OptionGroup(parser, "Required options")

    required_options.add_option('-o','--output', type="string",
                                help='the output folder for the simulation results')

    parser.add_option_group(required_options)

    optional_options = OptionGroup(parser, "Optional options")

    optional_options.add_option('--fit_method', default
    ="basinhopping", type="choice", choices=['basinhopping', 'differential_evolution', 'brute'],
                                help="Global optimization_method to use [default:%default]")

    optional_options.add_option('--pcoa_qza', type="string", default=None,
                                help='The input PCoA qza file')
    optional_options.add_option('--metadata', type="string", default=None,
                                help='The input metadata tsv file, if not defined, '
                                     'metadata will be extracted from PCoA qza')

    # Individual identifier will either accept a single column, or optionally combine two columns of categorical features.
    optional_options.add_option('--individual', type="string", default=None,
                                help='The host subject ID column identifier, separate by comma to combine TWO columns')
    optional_options.add_option('--timepoint', type="string", default=None,
                                help='The timepoint column identifier')
    optional_options.add_option('--treatment', type="string", default=None,
                                help='The treatment column identifier')

    optional_options.add_option('-v', '--verbose', action="store_true", dest="verbose", default=False,
                                help='-v, allows for verbose output' +
                                     ' [default: %default]')

    parser.add_option_group(optional_options)

    return parser

def fit_OU_process(data,dts):
    """Return the parameters of an OU process over data
    
    Strategy: this method combines two tools: parameteric
    fitting of normal distributions, and non-parametric 
    global optimization.

    1. generate a relativized version of all data, 
    describing the difference from the last time point. 
    This is needed to fit the Brownian motion/Weiner process
    part of the OU model, which is wholly location-independent.
   

    2.Call OU_fit to get an estimate of the log-likelihood

    If we already know the lambda and theta parameters,
    we should be able to fit the scale of the normal distribution.
    So, in each iteration of the global optimizer, we take input
    lambda and theta, and parametrically fit sigma to the data, 
    returning a negative log-likelihood
    
    To summarize, for fixed lambda, and theta, we get out
    a negative log-likelihood for the best parametric 
    estimate of the normal distribution describing the data.

    3. We can then pass this to basinhopping, getting out a best
    nLogLik for the overall optimization. This can then be fairly
    compared for e.g. BM vs. OU processes using Akaike's Information Criterion

    """
    pass

def get_OU_nlogLik(x,times,Sigma,Lambda,Theta):
    """Return the negative log likelihood for an OU model given data
    x -- an array of x values, which MUST be ordered by time
    times -- the time values at which x was taken. MUST match x in order.
    
    OU model parameters
    Sigma -- estimated Sigma for OU model (extent of change over time)
    Lambda -- estimated Lambda for OU model (tendency to return to average position)
    Theta -- estimated Theta for OU model (average or 'home' position)
    """
    dx = diff(x) 
    dts = diff(times)
    dx_dt = dx/dts #array of changes per unit time
    
    #To match the shape of x and dx_dt, we ignore
    #the last x, which doesn't produce a corresponding
    #dx_dt
    x = x[:-1]    

    #For each rel_change, the extent of change
    #should be dx/dt = W*sigma + lambda(theta - x)       
    #We should have constant vectors for dx/dt and x, 
    # and fixed estimates for lambda, theta, sigma.
    #Rearranging we get:
    W = ((dx_dt) - Lambda*(Theta - x)) #W is inclusive of Sigma so Ws

    #Current version of the Weiner process
    #is fixed in location, with scale depedent only on Sigma and dt.
    #We know these parameters, so the normal distribution for
    #the Weiner process at each data point is fixed given input params!
    dts_list = dts.tolist()
    dx_list = dx.tolist()
    W_list = W.tolist() #Weiner process rolls
    #For each step, the dx/dt will be based
    #on x in the *previous* step, so offset x's by one
    total_nlogLik = None

    #Sum up logLik's for all the differently scaled, 0-centered
    #normal distributions across all timepoints
    #NOTE: for evenly spaced data a more effecient 
    #single vector-based computation could be used.
    for i,dt in enumerate(dts_list):
        if i == 1:
            #First value is not informative
            continue
        loc = 0
        eps = 0.0000001
        scale = Sigma**(2*dt)+eps
        curr_W = W_list[i]
        nlogLik = -1*norm.logpdf(curr_W,loc=loc,scale=scale)
        if not total_nlogLik:
            total_nlogLik = nlogLik
        else:
            total_nlogLik += nlogLik #can validly sum log-likelihoods.
    return total_nlogLik

def make_OU_objective_fn(x,times,verbose=False):
    """Make an objective function for use with basinhopping with data embedded
    
    scipy.optimize.basinhopping needs a single array p, a function, and will
    minimize f(p). So we want to embded or dx data and time data *in* the 
    function, and use the values of p to represent parameter values that
    could produce the data. 
    """
    # x and times for each subject/site
    # x is just the individual x's, x1 and timepoints, x2 and timepoints, x3 and timepoints
    if len(x) <= 1:
        raise ValueError("Single-point observations are not supported.")

    def fn_to_optimize(p):
        fixed_x = x
        fixed_times = times
        if p.shape != (3,):
            raise ValueError("OU optimization must operate on a (3,) array representing Sigma,Lamda,and Theta values")
        #For clarity, binding these to variables
        Sigma,Lambda,Theta = p
        #print([fixed_x,fixed_times,Sigma,Lambda,Theta])
        nlogLik = get_OU_nlogLik(fixed_x,fixed_times,Sigma,Lambda,Theta)
        if verbose:
            print ("\nnlogLik:",nlogLik)
            #print "\t".join(["Sigma","Lambda","Theta"])
            print ("%.2f\t%.2f\t%.2f" % (Sigma,Lambda,Theta))
        return nlogLik
    return fn_to_optimize

#May not actually need the parametric
#normal distrubtion fit

def fit_normal(data):
    """Return the mean and standard deviation of normal data"""
    estimate = norm.fit(data)
    nlogLik = norm.nnlf(estimate,data)
    mu,std = estimate
    return mu,std,nlogLik



def fit_timeseries(fn_to_optimize,x0,xmin=array([-inf,-inf,-inf]),\
  xmax=array([inf,inf,inf]),global_optimizer="basinhopping",\
  local_optimizer="Nelder-Mead",stepsize=0.01,niter=200):
    """Minimize a function returning input & result
    fn_to_optimize -- the function to minimize. Must return a single value or array x.

    x0 -- initial parameter value or array of parameter values
    xmax -- max parameter values (use inf for infinite)
    xmin -- min parameter values (use -inf for infinite)
    global_optimizer -- the global optimization method (see scipy.optimize)
    local_optimizer -- the local optimizer (must be supported by global method)
    """

    if global_optimizer == "basinhopping":
        local_min_bounds = list(zip(xmin.tolist(),xmax.tolist()))
        local_min_kwargs = {"method":local_optimizer,"bounds":local_min_bounds}

        #From the scipy docs, if we want bounds on
        #our basinhopping steps, we need to implement a custom
        #accept test

        class Bounds(object):
          def __init__(self):
            self.xmax = array(xmax)
            self.xmin = array(xmin)
          def __call__(self, **kwargs):
            x = kwargs["x_new"]
            tmax = bool(all(x <= self.xmax))
            tmin = bool(all(x >= self.xmin))
            return tmax and tmin

        bounds_test = Bounds()
        result = basinhopping(fn_to_optimize,x0,minimizer_kwargs=local_min_kwargs,niter = niter,stepsize=0.05,accept_test = bounds_test)
        global_min = result.x
        #Result of evaluating fn_to_optimize at global min
        f_at_global_min = result.fun
    else:
        raise NotImplementedError("Only basinhopping is currently supported")

    return global_min,f_at_global_min


def parse_pcoa(pcoa_qza, individual, timepoint, treatment, metadata):
    """
    Load data from PCoA output in Qiime2 Format
    :param pcoa_qza: Location of PCoA file
    :param individual: Subject column identifier(s) [ex: Subject; Subject,BodySite]
    :param timepoint: Timepoint column identifier
    :param treatment: Treatment column identifier
    :param metadata: optionally defined metadata file location, if not defined, will use metadata from PCoA.qza
    :return: input dataframe, tsv filepath location
    """

    # Checks for column definitions
    if individual is None:
        raise ValueError('Individual Column name must be identified.')
    if " " in individual:
        raise ValueError("Remove all white space from subject identifier.")
    if timepoint is None:
        raise ValueError('Timepoint Column name must be identified.')
    if treatment is None:
        raise ValueError('Treatment Column name must be identified.')

    # elements of PCoA file are defined here, but only pc1, pc2, and pc3 from site are extracted at this time.
    eigs = []
    propEx = []
    species = []
    site = []

    # Open the zip files
    qza = zipfile.ZipFile(pcoa_qza, 'r')
    pcoa = qza.open(qza.namelist()[0].split('/')[0]+'/data/ordination.txt')
    # Define the tsv file location, if not set in opts.metadata
    if metadata is None:
        metadata = qza.open(qza.namelist()[0].split('/')[0]+'/provenance/action/metadata.tsv')
    else:
        pass
    # Parse PCoA Contents
    lines = pcoa.readlines()
    i = 0
    for line in lines:
        line = line.decode("utf-8")
        if line.startswith("Eigvals"):
            eigs = lines[i+1].decode("utf-8").strip('\n').split("\t")
        elif line.startswith("Proportion explained"):
            propEx = lines[i+1].decode("utf-8").strip('\n').split("\t")
        elif line.startswith("Species"):
            species = lines[i + 1].decode("utf-8").strip('\n').split("\t")

        elif line.startswith("Site"):
            # We don't want Site constraints.
            if ("constraints") in line:
                break
            max = int(line.split("\t")[1])+i
            j = i + 1
            while j <= max:
                t_line = lines[j].decode('utf-8')
                site.append(t_line.strip('\n').split("\t"))
                j += 1
        i += 1
    t_site = []
    for item in site:
        t_site.append([item[0],item[1],item[2],item[3]])
    site = t_site

    #We now have variables:
            #First three are stored for now, will be utilized later when standardization of axes can be supported.
        # eigs: every eigenvalue
        # propEx: every proportion explained
        # species: every species

        # site: Every site, with the highest three proportion-explained PCs, in the form:
            # subjectID, x1, x2, x3

    return site, metadata

def parse_metadata(metadata, individual, timepoint, treatment, site):
    """
    Parse relevant contents from metadata file to complete input dataframe
    :param metadata: tsv file location
    :param individual: subject column identifier
    :param timepoint: timepoint column identifier
    :param treatment: treatment column identifier
    :param site: [subjectID, x1, x2, x3]
    :return: input dataframe: [#SampleID,individual,timepoint,treatment,pc1,pc2,pc3]
    """

    # metadata file will be located in provenance/action/metadata.tsv, opens from QZA if not defined in opts
    df = pd.read_csv(metadata, sep="\t")

    # Drop any rows that are informational
    while df.iloc[0][0].startswith("#"):
        df.drop(df.index[:1], inplace=True)

    # Combine Individual columns if multiple subject identifiers defined (Such as individual and site)
    if "," in individual:
        individual_temp = individual.split(",")
        df_ind = df[individual_temp[0]]

        # Iterate over the columns and combine on each pass
        for i in range(len(individual_temp)):
            if i == 0:
                pass
            else:
                df_ind = pd.concat([df_ind,df[individual_temp[i]]],axis=1)

        # Hyphenate the white space between identifiers
        for i in range(len(individual_temp)):
            if i == 0:
                pass
            else:
                df_ind = df_ind[individual_temp[0]].astype(str)+"-"+df_ind[individual_temp[1]]

        # Remove any user-generated white space
        inds = []
        for row in df_ind.iteritems():
            inds.append(row[1].replace(" ","-"))
        df_ind = pd.DataFrame({individual:inds})
        # Reindex to match dataframes to be merged
        df_ind.index += 1
    else:
        df_ind = df[individual]

    df_tp = df[timepoint]
    df_tx = df[treatment]

    # Force timepoint to numerics
    df_tp = df_tp.replace('[^0-9]','', regex=True)

    # Build the return dataframe
    df_ret = pd.concat([df[df.columns[0]], df_ind, df_tp, df_tx], axis=1)
    df_site = pd.DataFrame.from_records(site, columns=["#SampleID","pc1","pc2","pc3"])
    df_ret = pd.merge(df_ret, df_site, on="#SampleID", how='inner')

    if "," in individual:
        df_ret.rename(columns = {0:str(individual)}, inplace=True)

    #Now have a full dataframe with sampleIDs, Subject Identifiers, timepoints, treatment, and initial x coordinates.
        #Preserves only values which have initial positions, drops null coordinates.
        #In example files, L3S242 was in metadata File, but not in the PCOA file, so there was one less in the df_ret.

    return df_ret

def fit_input(input, ind, tp, tx, method):
    """
    Fit and minimize the timeseries for each subject defined
    :param input: dataframe: [#SampleID,individual,timepoint,treatment,pc1,pc2,pc3]
    :param ind: subject column identifier
    :param tp: timepoint column identifier
    :param tx: treatment column identifier
    :param method: "basinhopping" if not defined in opts
    :return: dataframe: [ind,"pc","sigma","lambda","theta","optimizer","nLogLik","n_parameters","aic",tp,tx,"x"]
    """

    subjects = input[ind].unique()
    values = ['pc1','pc2','pc3']
    data = []
    for sj in subjects:
        df = input[input[ind] == sj]
        tp_t = df[tp].tolist()
        tp_t = [float(item) for item in tp_t]
        # Check for duplicate timepoints. If the subject identifier is not properly defined,
        # duplicate entries will be found in the timepoints.
        duplicates = [item for item, count in collections.Counter(tp_t).items() if count > 1]
        if duplicates:
            raise ValueError('Duplicate timepoint entries. Try redefining your subject identifiers.\n'
                             'Duplicates: '+str(duplicates))
        tx_t = df[tx].tolist()
        vals = []

        # Pivot the PCs
        for value in values:
            numeric = df[value].tolist()
            numeric = [float(item) for item in numeric]
            vals.append([numeric, value])
        data.append([sj,tp_t,tx_t,vals])

    # Make an objective function for each row
    fx = []
    for row in data:
        for item in row[3]:
            #[0: function, 1: subject, 2: times, 3: treatment, 4: values, 5:pc1/2/3]
            fx.append([make_OU_objective_fn(item[0], row[1]), row[0], row[1], row[2], item[0], item[1]])
    fit_ts = []

    # Optimize each function
    for row in fx:
        x = np.asarray(row[4])
        if len(x)>1:
            fit_ts.append([fit_timeseries(row[0], [0.1, 0.0, np.mean(x)], global_optimizer=method),
                           row[1], row[2], row[3], x, row[5]])

    # Generate output data from optimized functions.
    sig, lam, the, nlogLik, ind_o, tp_o, tx_o, x_o, pc, op, n_param, aic = ([] for i in range(12))
    for row in fit_ts:
        # global_min = [sigma, lambda, theta]
        sig.append(row[0][0][0])
        lam.append(row[0][0][1])
        the.append(row[0][0][2])
        # f_at_global min = nlogLik
        nlogLik.append(row[0][1])
        ind_o.append(row[1])
        tp_o.append(row[2])
        tx_o.append(row[3])
        x_o.append(row[4])
        n_param.append(len(row[4]))
        pc.append(row[5])
        op.append(method)
        # aic calulated with 2*n_params-2*LN(-1*nloglik)
        aic.append(2*len(row[4])-2*(math.log(-1*row[0][1])))
    output = pd.DataFrame({"sigma":sig,
                           "lambda":lam,
                           "theta":the,
                           "nLogLik":nlogLik,
                           ind:ind_o,
                           tp:tp_o,
                           tx:tx_o,
                           "x": x_o,
                           "n_parameters":n_param,
                           "pc":pc,
                           "optimizer":op,
                           "aic":aic},
                           columns=[ind,"pc","sigma","lambda","theta","optimizer",
                                    "nLogLik","n_parameters","aic",tp,tx,"x"])
    return output

def benchmark(output):
    """
    Verifies that fit_timeseries recovers OU model params
    :param output: location for output log
    """
    f = open(output+"fit_timeseries_benchmark_log.txt","w+")

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
    for n_timepoints in list(range(1, 300)):
        print("Building OU model for %i timepoints" % n_timepoints)
        # print("Building OU model for %i timepoints" % n_timepoints, file=f)

        # run ou_process to get history
        ou = Process(start_coord=0.20, motion="Ornstein-Uhlenbeck", \
                     history=None, params= {"lambda": 0.12, "delta": 0.25, "mu": 0.5})
        for t in range(0, n_timepoints):
            ou.update(dt)
        print(n_timepoints, ou.History)
        # print(n_timepoints, ou.History, file=f)
        xs = array(ou.History)
        ts = np.arange(0, len(ou.History)) * dt
        print(xs, ts, dt)
        # print(xs, ts, dt, file=f)
        fn_to_optimize = make_OU_objective_fn(xs, ts)

        # Estimate correct parameters
        for niter in [5]:
            for local_optimizer in ['L-BFGS-B']:
                print("Running optimizer:", local_optimizer)
                # print("Running optimizer:", local_optimizer, file=f)
                # Using intentionally kinda bad estimates
                start_Sigma = 0.1
                start_Lambda = 0.0
                start_Theta = np.mean(xs)
                print("niter=", niter)
                # print("niter=", niter, file=f)
                print("start_Theta: ", start_Theta)
                # print("start_Theta: ", start_Theta, file=f)
                print("n_timepoints: ", n_timepoints)
                # print("n_timepoints: ", n_timepoints, file=f)
                xmax = array([1.0, 1.0, 1.0])
                xmin = array([0.0, 0.0, -1.0])
                x0 = array([start_Sigma, start_Lambda, start_Theta])

                global_min, f_at_global_min = \
                    fit_timeseries(fn_to_optimize, x0, xmin, xmax, stepsize=0.005, \
                                   niter=niter, local_optimizer=local_optimizer)

                print("OU result:")
                # print("OU result:", file=f)
                Sigma, Lambda, Theta = global_min
                correct_values = array([0.25, 0.12, 0.5])
                final_error = global_min - correct_values
                print("Global min:", global_min)
                # print("Global min:", global_min, file=f)
                print("f at Global min:", f_at_global_min)
                # print("f at Global min:", f_at_global_min, file=f)
                # aic calulated with 2*n_params-2*LN(-1*nloglik)
                print("aic:" , (2 * niter - 2 * (math.log(-1 * f_at_global_min))))
                # print("aic:", (2 * niter - 2 * (math.log(-1 * f_at_global_min))), file=f)

                final_errors["%s_%i_%i" % (local_optimizer, niter, n_timepoints)] = final_error
                print("*" * 80)
                # print("*" * 80, file=f)
                print("%s error: %.4f,%.4f,%.4f" % (local_optimizer,final_error[0],final_error[1],final_error[2]))
                # print("%s error: %.4f,%.4f,%.4f" % (local_optimizer, final_error[0], final_error[1], final_error[2]), file=f)
                print("*" * 80)
                # print("*" * 80, file=f)
                print()
                #print(file=f)
    for opt, err in final_errors.items():
        print("%s error: %.4f,%.4f,%.4f" % (opt, err[0], err[1], err[2]))
        #print("%s error: %.4f,%.4f,%.4f" %(opt,err[0],err[1],err[2]), file=f)

def main():
    parser = make_option_parser()
    opts, args = parser.parse_args()
    verbose = opts.verbose
    if verbose:
        print(opts)

    if opts.output is None:
        raise IOError('Output must be defined.')

    if not os.path.exists(opts.output):
        os.makedirs(opts.output)

    if opts.pcoa_qza is not None:
        # Extract the name of the pcoa qza, and append with '_fit_timeseries'
        out_name = opts.pcoa_qza.split("/")[-1].split(".")[0] + "_fit_timeseries.csv"
        site, metadata = parse_pcoa(opts.pcoa_qza, opts.individual, opts.timepoint, opts.treatment, opts.metadata)
        input = parse_metadata(metadata,opts.individual, opts.timepoint, opts.treatment, site)
        output = fit_input(input, opts.individual, opts.timepoint, opts.treatment, opts.fit_method)
        output.to_csv(opts.output+out_name)
    else:
        benchmark(opts.output)

    # FAIL:
    """
     -o ./data/fit_timeseries/
     --pcoa_qza ./data/fit_timeseries/2485_pcoa.qza
     --individual host_subject_id
     --timepoint sampledate
     --treatment obesitycat
     -v
     
     ValueError: Single-point observations are not supported.
    """

    """
    -o ./data/fit_timeseries/
    --pcoa_qza ./data/fit_timeseries/weighted_unifrac_pcoa_results.qza
    --individual Subject
    --timepoint DaysSinceExperimentStart
    --treatment ReportedAntibioticUsage
    -v
    
    ValueError: Duplicate timepoint entries. Try redefining your subject identifiers.
    Duplicates:[0.0, 84.0, 112.0, 140.0, 168.0]
    """

    # PASS:
    """
    -o ./data/fit_timeseries/
    --pcoa_qza ./data/fit_timeseries/weighted_unifrac_pcoa_results.qza
    --individual Subject,BodySite
    --timepoint DaysSinceExperimentStart
    --treatment ReportedAntibioticUsage
    -v
    """

    """
    -o ./data/fit_timeseries/
    --pcoa_qza ./data/fit_timeseries/unweighted_unifrac_pcoa_results.qza
    --individual Subject,BodySite
    --timepoint DaysSinceExperimentStart
    --treatment ReportedAntibioticUsage
    -v
    """

    """
    -o ./data/fit_timeseries/
    --pcoa_qza ./data/fit_timeseries/jaccard_pcoa_results.qza
    --individual Subject,BodySite
    --timepoint DaysSinceExperimentStart
    --treatment ReportedAntibioticUsage
    -v
    """

    """
    -o ./data/fit_timeseries/
    --pcoa_qza ./data/fit_timeseries/bray_curtis_pcoa_results.qza
    --individual Subject,BodySite
    --timepoint DaysSinceExperimentStart
    --treatment ReportedAntibioticUsage
    -v
    """



if __name__ == "__main__":
    main()
