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
import zipfile
import numpy as np
import pandas as pd
import collections
import os
import math

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

"""
def fit_OU_process(data,dts):
    Return the parameters of an OU process over data
    
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

    pass
"""

def get_OU_nLogLik(x,times,Sigma,Lambda,Theta):
    """
    Return the negative log likelihood for an OU model given data x -- an array of x values,
    which MUST be ordered by time times -- the time values at which x was taken.
    MUST match x in order.

    OU model parameters

    1. Sigma -- estimated Sigma for OU model (extent of change over time)
    2. Lambda -- estimated Lambda for OU model (tendency to return to average position)
    3. Theta -- estimated Theta for OU model (average or 'home' position)

    :param x: array of values ordered by time
    :param times: times associated with x values
    :param Sigma: initial sigma value
    :param Lambda: initial lambda value
    :param Theta: initial theta value
    :return: nLogLik value
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
    total_nLogLik = None

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
        nLogLik = -1*norm.logpdf(curr_W,loc=loc,scale=scale)
        if not total_nLogLik:
            total_nLogLik = nLogLik
        else:
            total_nLogLik += nLogLik #can validly sum log-likelihoods.
    return total_nLogLik

def make_OU_objective_fn(x,times,verbose=False):
    """
    Make an objective function for use with basinhopping with data embedded

    scipy.optimize.basinhopping needs a single array p, a function, and will
    minimize f(p). So we want to embded or dx data and time data *in* the
    function, and use the values of p to represent parameter values that
    could produce the data.

    :param x: individual x values for objective function
    :param times: timepoints associated with passed-in x values
    :param verbose: verbose output, default = False
    :return: fn_to_optimize
    """
    # x and times for each subject/site
    # x is just the individual x's, x1 and timepoints, x2 and timepoints, x3 and timepoints
    if len(x) <= 1:
        raise ValueError("Single-point observations are not supported.")

    def fn_to_optimize(p):
        fixed_x = x
        fixed_times = times
        if p.shape != (3,):
            raise ValueError(
                "OU optimization must operate on a (3,) array representing Sigma,Lamda,and Theta values")
        # For clarity, binding these to variables
        Sigma, Lambda, Theta = p
        # print([fixed_x,fixed_times,Sigma,Lambda,Theta])
        nLogLik = get_OU_nLogLik(fixed_x, fixed_times, Sigma, Lambda, Theta)
        if verbose:
            print("\nnLogLik:", nLogLik)
            # print "\t".join(["Sigma","Lambda","Theta"])
            print("%.2f\t%.2f\t%.2f" % (Sigma, Lambda, Theta))
        return nLogLik
    return fn_to_optimize

#May not actually need the parametric
#normal distrubtion fit

def make_OU_objective_fn_cohort(x, times, verbose=False):
    """
    Sums nLogLik and build objective function based on cohorts. Operates in the
    same manner as make_OU_objective_fn, except that it considers a treatment cohort,
    not just the individuals.

    Make an objective function for use with basinhopping with data embedded

    scipy.optimize.basinhopping needs a single array p, a function, and will
    minimize f(p). So we want to embded or dx data and time data *in* the
    function, and use the values of p to represent parameter values that
    could produce the data.

    Overall strategy: Treat this exactly like the per-individual
    fitting, BUT within the objective function loop over all individuals
    in a given treatment (not just one) to get nLogLiklihoods. The nLogLiklihood
    for the individuals in the treatment is then just the sum of the individual
    nLogLikelihoods for each individuals timeseries.

    data -- a dict of {"individual1": (fixed_x,fixed_times)}

    :param x: x values for cohort objective function
    :param times: times associated with x values from cohort
    :param verbose: verbose output, default = False
    :return: fn_to_optimize
    """

    if len(x) <= 1:
        raise ValueError("Single-point observations are not supported.")

    def fn_to_optimize(p):
        fixed_x = x
        fixed_times = times
        if(len(fixed_x) != len(fixed_times)):
            raise ValueError("Dimensions of X and Times must be equal")
        if p.shape != (3,):
            raise ValueError(
                "OU optimization must operate on a (3,) array representing Sigma,Lamda,and Theta values")
        # For clarity, binding these to variables
        Sigma, Lambda, Theta = p
        #Assuming data is a dict of indiviudal_id: 
        #start with nLogLik = 0
        cohort_nLogLik = 0
        for i in range(len(fixed_x)):
            # print([fixed_x,fixed_times,Sigma,Lambda,Theta])
            nLogLik = get_OU_nLogLik(fixed_x[i], fixed_times[i], Sigma, Lambda, Theta)
            if verbose:
                print("\nnLogLik:", nLogLik)
                # print "\t".join(["Sigma","Lambda","Theta"])
                print("%.2f\t%.2f\t%.2f" % (Sigma, Lambda, Theta))
            cohort_nLogLik += nLogLik #test by passing in individuals and summing
        if verbose:
            print("\nCohort nLogLik:", nLogLik)
            # print "\t".join(["Sigma","Lambda","Theta"])
            print("%.2f\t%.2f\t%.2f" % (Sigma, Lambda, Theta))
        return cohort_nLogLik
    return fn_to_optimize


def fit_normal(data):
    """
    Return the mean and standard deviation of normal data

    :param data: fit data to normal distribution
    :return:  mu, std, nLogLik
    """
    estimate = norm.fit(data)
    nLogLik = norm.nnlf(estimate,data)
    mu,std = estimate
    return mu,std,nLogLik


def fit_timeseries(fn_to_optimize,x0,xmin=array([-inf,-inf,-inf]),
  xmax=array([inf,inf,inf]),global_optimizer="basinhopping",
  local_optimizer="Nelder-Mead",stepsize=0.01,niter=200, verbose=False):
    """
    Minimize a function returning input & result
    fn_to_optimize -- the function to minimize. Must return a single value or array x.

    x0 -- initial parameter value or array of parameter values
    xmax -- max parameter values (use inf for infinite)
    xmin -- min parameter values (use -inf for infinite)
    global_optimizer -- the global optimization method (see scipy.optimize)
    local_optimizer -- the local optimizer (must be supported by global method)

    :param fn_to_optimize: function that is being optimized, generated from fn_to_optimize
    :param x0: initial parameter value or array of parameter values
    :param xmin: min parameter values (-inf for infinite)
    :param xmax: max parameter values (inf for infinite)
    :param global_optimizer: global optimization method
    :param local_optimizer: local optimization method (must be supported by global)
    :param stepsize: size for each step (.01)
    :param niter: number of iterations (200)
    :param verbose: verbose output, default = False
    :return: global_min, f_at_global_min
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
        result = basinhopping(fn_to_optimize, x0, minimizer_kwargs=local_min_kwargs, niter = niter,stepsize=0.05,
                              accept_test = bounds_test)
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
            eigs = lines[i+1].decode("utf-8").strip('\n').strip('\r').split("\t")
        elif line.startswith("Proportion explained"):
            propEx = lines[i+1].decode("utf-8").strip('\n').strip('\r').split("\t")
        elif line.startswith("Species"):
            species = lines[i + 1].decode("utf-8").strip('\n').strip('\r').split("\t")

        elif line.startswith("Site"):
            # We don't want Site constraints.
            if ("constraints") in line:
                break
            max = int(line.split("\t")[1])+i
            j = i + 1
            while j <= max:
                t_line = lines[j].decode('utf-8')
                site.append(t_line.strip('\n').strip('\r').split("\t"))
                j += 1
        i += 1
    t_site = []
    for item in site:
        t_subj_id = item[0]
        t_pc1 = item[1]
        t_pc2 = item[2]
        t_pc3 = item[3]
        t_site.append([t_subj_id, t_pc1, t_pc2, t_pc3])
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
    #if "," in individual:
    #    individual = individual.replace(",","_")
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
                df_ind = df_ind[individual_temp[0]].astype(str)+"_"+df_ind[individual_temp[1]]

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
    
    #Copy individual column into treatment if none assigned
    if treatment is not None:
        df_tx = df[treatment]
    else:
        df_tx = df_ind.copy()
        df_tx.columns=[str(df_tx.columns.values[0])+"_tx"]

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

def aic(n_param, nLogLik):
    """
    calculates AIC with 2 * n_parameters - 2 * LN(-1 * nLogLik)

    :param n_param: number of parameters
    :param nLogLik: negative log likelihood
    :return: aic score
    """
    return 2*n_param-2*(math.log(-1*nLogLik))


def gen_output(fit_ts, ind, tp, tx, method):
    """
    Generate output dataframe for either cohort, or non-cohort data

    :param fit_ts: dataframe [0: [[Sigma, Lambda, Theta], nLogLik], 1: Individuals, 2: Times 3: Treatments, 4: Values, 5: PC Axis
    :param ind: individual identifier
    :param tp: timepoint identifier
    :param tx: treatment identifier
    :param method: optimization method
    :return: Formatted dataframe for csv output
    """
    sig, lam, the, nLogLik, ind_o, tp_o, tx_o, x_o, pc, op, n_param, aic_o = ([] for i in range(12))
    for row in fit_ts:
        # global_min = [sigma, lambda, theta]
        sig.append(row[0][0][0])
        lam.append(row[0][0][1])
        the.append(row[0][0][2])
        # f_at_global min = nLogLik
        nLogLik.append(row[0][1])
        ind_o.append(row[1])
        tp_o.append(row[2])
        tx_o.append(row[3])
        x_o.append(row[4].tolist())
        n_param.append(len(row[4].tolist()))
        pc.append(row[5])
        op.append(method)
        # aic calulated with 2*n_params-2*LN(-1*nLogLik)
        aic_o.append(aic(len(row[4].tolist()), row[0][1]))

    if "," in ind:
        ind = ind.replace(",", "_")

    output = pd.DataFrame({"sigma": sig,
                           "lambda": lam,
                           "theta": the,
                           "nLogLik": nLogLik,
                           ind: ind_o,
                           "n_parameters": n_param,
                           "pc": pc,
                           "optimizer": op,
                           "aic": aic_o},
                          columns=[ind, "pc", "sigma", "lambda", "theta", "nLogLik",
                                   "n_parameters", "aic", "optimizer"])
    return output


def fit_cohorts(input, ind, tp, tx, method, verbose = False):
    """
    Completes the same operation as fit_input, for cohorts.
    Fit and minimize the timeseries for each subject and cohort defined.

    :param input: dataframe: [#SampleID,individual,timepoint,treatment,pc1,pc2,pc3]
    :param ind: subject column identifier
    :param tp: timepoint column identifier
    :param tx: treatment column identifier
    :param method: "basinhopping" if not defined in opts
    :return: dataframe: [ind,"pc","sigma","lambda","theta","optimizer","nLogLik","n_parameters","aic",tp,tx,"x"]
    """

    # Need to pass cohort data to new make_ou_objective_fn_cohort(data, verbose=false)
    # Split the data into each treatment, for each treatment_set, return an objective fn,
    #   assign each row their associated objective fn
    fx_cohorts = []
    #Cycle through all treatments
    for treatment in input[tx].unique():
        df = input[input[tx] == treatment]

        # Create empty dataframes for treatment PC axes, and timepoints
        pc1_values = []
        pc2_values = []
        pc3_values = []
        pc_times = []

        # Append values and times for each individual
        for subject in df[ind].unique():
            pc1_values.append(df[df[ind] == subject]["pc1"].values.tolist())
            pc2_values.append(df[df[ind] == subject]["pc1"].values.tolist())
            pc3_values.append(df[df[ind] == subject]["pc1"].values.tolist())
            pc_times.append(df[df[ind] == subject][tp].values.tolist())

        # Force the values to float
        sub_pc1 = []
        for item in pc1_values:
            sub_pc1.append([float(i) for i in item])
        pc1_values = sub_pc1

        sub_pc2 = []
        for item in pc2_values:
            sub_pc2.append([float(i) for i in item])
        pc2_values = sub_pc2

        sub_pc3 = []
        for item in pc3_values:
            sub_pc3.append([float(i) for i in item])
        pc3_values = sub_pc3

        sub_times = []
        for item in pc_times:
            sub_times.append([float(i) for i in item])
        pc_times = sub_times

        # Append objective functions and metadata to dataframe
        fx_cohorts.append([make_OU_objective_fn_cohort(pc1_values, pc_times), str(treatment + "_pc1"), pc_times,
                           treatment, pc1_values, "pc1"])
        fx_cohorts.append([make_OU_objective_fn_cohort(pc2_values, pc_times), str(treatment + "_pc2"), pc_times,
                           treatment, pc2_values, "pc2"])
        fx_cohorts.append([make_OU_objective_fn_cohort(pc3_values, pc_times), str(treatment + "_pc3"), pc_times,
                           treatment, pc3_values, "pc3"])

    fit_ts_cohorts = []

    # Fit each cohort's timeseries
    for row in fx_cohorts:
        x = np.asarray(row[4])
        if len(x) > 1:
            row_ou_fn = row[0]
            row_subj_id = row[1]
            row_times = row[2]
            row_tx = row[3]
            row_pc_val = x
            row_pc_axis = row[5]

            fit_ts_cohorts.append([fit_timeseries(row_ou_fn, [0.1, 0.0, np.mean(x)], global_optimizer=method),
                                   row_subj_id, row_times, row_tx, row_pc_val, row_pc_axis])

    return gen_output(fit_ts_cohorts, ind, tp, tx, method)

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
    if tx is not None:
        cohorts = True
    else:
        cohorts = False
        tx = str(ind+"_tx")
    subjects = input[ind].unique()
    values = ['pc1', 'pc2', 'pc3']
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
        data.append([sj, tp_t, tx_t, vals])

    fx = []
    for row in data:
        axis = row[3]
        for item in axis:
            #[0: function, 1: subject, 2: times, 3: treatment, 4: values, 5:pc1/2/3]
            axis_data = item[0]
            row_times = row[1]
            row_subj_id = row[0]
            row_tx = row[2]
            axis_label = item[1]
            fx.append([make_OU_objective_fn(axis_data, row_times), row_subj_id, row_times, row_tx, axis_data, axis_label])

    if cohorts:
        cohort_output = fit_cohorts(input, ind, tp, tx, method)

    fit_ts = []

    # Optimize each function
    for row in fx:
        x = np.asarray(row[4])
        if len(x)>1:
            row_ou_fn = row[0]
            row_subj_id = row[1]
            row_times = row[2]
            row_tx = row[3]
            row_val = x
            row_axis_label = row[5]
            fit_ts.append([fit_timeseries(row_ou_fn, [0.1, 0.0, np.mean(x)], global_optimizer=method),
                           row_subj_id, row_times, row_tx, row_val, row_axis_label])

    output = gen_output(fit_ts, ind, tp, tx, method)

    if cohorts:
        return output, cohort_output
    else:
        return output


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
        cohort_out_name = opts.pcoa_qza.split("/")[-1].split(".")[0] + "_fit_timeseries_cohorts.csv"
        site, metadata = parse_pcoa(opts.pcoa_qza, opts.individual, opts.timepoint, opts.treatment, opts.metadata)
        input = parse_metadata(metadata,opts.individual, opts.timepoint, opts.treatment, site)
        if opts.treatment is not None:
            output, cohort_output = fit_input(input, opts.individual, opts.timepoint, opts.treatment, opts.fit_method)
            cohort_output.to_csv(opts.output+cohort_out_name, index=False)
        else:
            output = fit_input(input, opts.individual, opts.timepoint, opts.treatment, opts.fit_method)
        output.to_csv(opts.output+out_name, index=False)
    else:
        parser.print_help()
        exit(0)

if __name__ == "__main__":
    main()
