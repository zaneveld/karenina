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

import qiime2
import q2templates
from q2_types.ordination import PCoAResults

TEMPLATES = pkg_resources.resource_filename('q2_karenina')
"""Changes: 
	removed options definitions
	changed each non-q2 visualization / method function prefix to "_"
	adjusted variable definitions
"""
def fit_timeseries(pcoa : PCoAResults, metadata : qiime2.Metadata,
					method : str, individual: str, timepoint: str, treatment: str):
    
	# Extract the name of the pcoa qza, and append with '_fit_timeseries'
	out_name = pcoa.split("/")[-1].split(".")[0] + "_fit_timeseries.csv"
	cohort_out_name = pcoa.split("/")[-1].split(".")[0] + "_fit_timeseries_cohorts.csv"
	site, metadata = _parse_pcoa(pcoa, individual, timepoint, treatment, metadata)
	input = _parse_metadata(metadata,individual, timepoint, treatment, site)
	if treatment is not None:
		output, cohort_output = _fit_input(input, individual, timepoint, treatment, method)
		cohort_output.to_csv(opts.output+cohort_out_name, index=False)
	else:
		output = _fit_input(input, individual, timepoint, treatment, fit_method)
	output.to_csv(opts.output+out_name, index=False)
	#q2templates.render()????


def _get_OU_nLogLik(x,times,Sigma,Lambda,Theta):
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

def _make_OU_objective_fn(x,times,verbose=False):
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
            raise ValueError(
                "OU optimization must operate on a (3,) array representing Sigma,Lamda,and Theta values")
        # For clarity, binding these to variables
        Sigma, Lambda, Theta = p
        # print([fixed_x,fixed_times,Sigma,Lambda,Theta])
        nLogLik = _get_OU_nLogLik(fixed_x, fixed_times, Sigma, Lambda, Theta)
        if verbose:
            print("\nnLogLik:", nLogLik)
            # print "\t".join(["Sigma","Lambda","Theta"])
            print("%.2f\t%.2f\t%.2f" % (Sigma, Lambda, Theta))
        return nLogLik
    return fn_to_optimize

#May not actually need the parametric
#normal distrubtion fit

def _make_OU_objective_fn_cohort(x, times, verbose=False):
    """Make an objective function for use with basinhopping with data embedded

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
            nLogLik = _get_OU_nLogLik(fixed_x[i], fixed_times[i], Sigma, Lambda, Theta)
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


def _fit_normal(data):
    """Return the mean and standard deviation of normal data"""
    estimate = norm.fit(data)
    nLogLik = norm.nnlf(estimate,data)
    mu,std = estimate
    return mu,std,nLogLik


def _fit_timeseries(fn_to_optimize,x0,xmin=array([-inf,-inf,-inf]),
  xmax=array([inf,inf,inf]),global_optimizer="basinhopping",
  local_optimizer="Nelder-Mead",stepsize=0.01,niter=200, verbose=False):
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
        result = basinhopping(fn_to_optimize, x0, minimizer_kwargs=local_min_kwargs, niter = niter,stepsize=0.05,
                              accept_test = bounds_test)
        global_min = result.x
        #Result of evaluating fn_to_optimize at global min
        f_at_global_min = result.fun
    else:
        raise NotImplementedError("Only basinhopping is currently supported")

    return global_min,f_at_global_min


def _parse_pcoa(pcoa_qza, individual, timepoint, treatment, metadata):
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

def _parse_metadata(metadata, individual, timepoint, treatment, site):
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

def _aic(n_param, nLogLik):
    """
    calculates AIC with 2 * n_parameters - 2 * LN(-1 * nLogLik)
    :param n_param: number of parameters
    :param nLogLik: negative log likelihood
    :return: aic score
    """
    return 2*n_param-2*(math.log(-1*nLogLik))


def _gen_output(fit_ts, ind, tp, tx, method):
    """
    Generate output dataframe for either cohort, or non-cohort data
    :param fit_ts: dataframe [0: [[Sigma, Lambda, Theta], nLogLik], 1: Individuals, 2: Times
        3: Treatments, 4: Values, 5: PC Axis
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
        aic_o.append(_aic(len(row[4].tolist()), row[0][1]))

    if "," in ind:
        ind = ind.replace(",", "_")

    output = pd.DataFrame({"sigma": sig,
                           "lambda": lam,
                           "theta": the,
                           "nLogLik": nLogLik,
                           ind: ind_o,
                           tp: tp_o,
                           tx: tx_o,
                           "x": x_o,
                           "n_parameters": n_param,
                           "pc": pc,
                           "optimizer": op,
                           "aic": aic_o},
                          columns=[ind, "pc", "sigma", "lambda", "theta", "nLogLik",
                                   "n_parameters", "aic", "optimizer", tp, tx, "x"])
    return output


def _fit_cohorts(input, ind, tp, tx, method, verbose = False):
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
        fx_cohorts.append([_make_OU_objective_fn_cohort(pc1_values, pc_times), str(treatment + "_pc1"), pc_times,
                           treatment, pc1_values, "pc1"])
        fx_cohorts.append([_make_OU_objective_fn_cohort(pc2_values, pc_times), str(treatment + "_pc2"), pc_times,
                           treatment, pc2_values, "pc2"])
        fx_cohorts.append([_make_OU_objective_fn_cohort(pc3_values, pc_times), str(treatment + "_pc3"), pc_times,
                           treatment, pc3_values, "pc3"])

    fit_ts_cohorts = []

    # Fit each cohort's timeseries
    for row in fx_cohorts:
        x = np.asarray(row[4])
        if len(x) > 1:
            fit_ts_cohorts.append([_fit_timeseries(row[0], [0.1, 0.0, np.mean(x)], global_optimizer=method),
                                   row[1], row[2], row[3], x, row[5]])

    return _gen_output(fit_ts_cohorts, ind, tp, tx, method)

def _fit_input(input, ind, tp, tx, method):
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
        for item in row[3]:
            #[0: function, 1: subject, 2: times, 3: treatment, 4: values, 5:pc1/2/3]
            fx.append([_make_OU_objective_fn(item[0], row[1]), row[0], row[1], row[2], item[0], item[1]])

    if cohorts:
        cohort_output = _fit_cohorts(input, ind, tp, tx, method)

    fit_ts = []

    # Optimize each function
    for row in fx:
        x = np.asarray(row[4])
        if len(x)>1:
            fit_ts.append([_fit_timeseries(row[0], [0.1, 0.0, np.mean(x)], global_optimizer=method),
                           row[1], row[2], row[3], x, row[5]])

    output = _gen_output(fit_ts, ind, tp, tx, method)

    if cohorts:
        return output, cohort_output
    else:
        return output
