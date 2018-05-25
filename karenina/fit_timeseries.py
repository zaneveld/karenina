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


def make_option_parser():
    """Return an optparse OptionParser object"""

    parser = OptionParser(usage = [
                               ("","Demo fitting an OU model using default parameters.", "%prog -o ./simulation_results")
                                ],
    description = "This script fits microbiome change over time to Ornstein-Uhlenbeck (OU) models.",
    version = __version__)

    required_options = OptionGroup(parser, "Required options")

    required_options.add_option('-o','--output', type="new_filepath",
    help='the output folder for the simulation results')

    parser.add_option_group(required_options)

    optional_options = OptionGroup(parser, "Optional options")

    optional_options.add_option('--fit_method', default
    ="basinhopping", type="choice", choices=['basinhopping', 'differential_evolution', 'brute'], help="Global optimization_method to use [default:%default]")

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
    
    def fn_to_optimize(p):
        fixed_x = x
        fixed_times = times
        if p.shape != (3,):
            raise ValueError("OU optimization must operate on a (3,) array representing Sigma,Lamda,and Theta values")
        #For clarity, binding these to variables
        Sigma,Lambda,Theta = p
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

 
def main():
    parser = make_option_parser()
    opts, args = parser.parse_args()
    if opts.verbose:
        print(opts)

    #TODO: Make if not exists output directory
    #if exists(opts.output):
    #    print("Output saved to: " +str(opts.output))


if __name__ == "__main__":
    main()
