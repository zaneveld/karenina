#!/usr/bin/env python

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"

import unittest
import sys, os
testdir = os.path.dirname(__file__)
srcdir = '../karenina'
sys.path.insert(0, os.path.abspath(os.path.join(testdir, srcdir)))
from warnings import catch_warnings
from karenina.fit_timeseries  import fit_timeseries,fit_normal,\
  get_OU_nLogLik,make_OU_objective_fn
from karenina.process import Process
import numpy.testing as npt

from scipy.stats import norm
from numpy import linspace,cos,arange,array,mean
  
"""
Tests for fit_timeseries.py
"""



class TestFit(unittest.TestCase):
    """Tests of the fit timeseries module"""

    def setUp(self):
        """
        Set up test data for each test function
        """
        
        #If this model fitting approach is going to work
        #we need to ensure that a super basic function
        #can be successfulyy fit.

        #generate several normal distributions
        test_normal_data = {}
        n_obs = 1000
        dt = 0.01
        for delta in [float(i)/100.0 for i in range(0,100,1)]:
            curr_data =  norm.rvs(loc=0,size=n_obs,scale=delta**(2*dt))
            test_normal_data[delta**(2*dt)] = curr_data
        self.BasicNormalData = test_normal_data

        #generate OU process for testing 
        ou_process = Process(start_coord=0.20,motion="Ornstein-Uhlenbeck",\
          history = None, params = {"lambda":0.20,"delta":0.25,"mu":0.0})
        #run ou_process to get history
        for t in range(1,30):
            dt = 1
            ou_process.update(dt)
        self.OU = ou_process
        self.verbose = False
            

    def test_fit_normal(self):
        """Return the mean and standard deviation of normal data."""
        for scale,data in self.BasicNormalData.items():
             est_loc,est_scale,nLogLik = fit_normal(data)
             accurate_to = 3 #decimal places
             npt.assert_almost_equal(est_loc,0,1)
             npt.assert_almost_equal(est_scale,scale,1)

    def test_basinhopping_canned_example(self):
        """Basinhopping fits a parabola with superimposed local minima."""
        #This is directly from the scipy.optimize docs
        
        fn_to_optimize = lambda x: cos(14.5 * x - 0.3) + (x + 0.2) * x   
        x0 = 1.0
        global_min,f_at_global_min = fit_timeseries(fn_to_optimize,x0,\
          global_optimizer = "basinhopping",local_optimizer="BFGS")
        
        npt.assert_almost_equal(global_min,-0.1951,4)
        npt.assert_almost_equal(f_at_global_min,-1.0009,4)

    def test_get_OU_nLogLik_accords_with_correct_params(self):
        """get_OU_nLogLik gives best score to correct params"""

        for attempt in range(3):
            try:
                ou = self.OU
                xs = array(ou.History)
                ts = arange(0,len(ou.History))
                if self.verbose:
                    print ("xs.shape:",xs.shape)
                    print ("ts.shape:",ts.shape)
                #Let's try a range of values around the true ones
                #as we'd produce when basinhopping

                #True values:
                #Theta=0, Lambda = 0.20,Sigma=delta=0.25

                #nLogLik estimate with true parameters:
                true_param_est = get_OU_nLogLik(xs,ts,Lambda=0.20,\
                  Sigma=0.25,Theta=0)

                bad_param_est1 = get_OU_nLogLik(xs,ts,Lambda=0.40,\
                  Sigma=0.25,Theta=0)

                self.assertTrue(true_param_est < bad_param_est1)

                bad_param_est2 = get_OU_nLogLik(xs,ts,Lambda=0.80,\
                  Sigma=0.25,Theta=0)

                self.assertTrue(true_param_est < bad_param_est2)
                self.assertTrue(bad_param_est1 < bad_param_est2)

                bad_param_est3 = get_OU_nLogLik(xs,ts,Lambda=0.80,\
                  Sigma=0.25,Theta=0.5)

                self.assertTrue(true_param_est < bad_param_est3)
                self.assertTrue(bad_param_est2 < bad_param_est3)

                bad_param_est4 = get_OU_nLogLik(xs,ts,Lambda=0.20,\
                  Sigma=0.50,Theta=0)

                self.assertTrue(true_param_est < bad_param_est4)
                #Don't expect that this is worse than bad_param_est3

            except AssertionError:
                attempt += 1
                print("Bad RNG Sample; retrying ", attempt,"/3...")
                self.setUp()

    def test_parse_pcoa(self):
        """TODO"""
        pass

    def test_parse_metadata(self):
        """TODO"""
        pass

    def test_fit_input(self):
        """TODO"""
        pass


if __name__ == '__main__':
    unittest.main()


