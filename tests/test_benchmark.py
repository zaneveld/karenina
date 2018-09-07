#!/usr/bin/env python

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Jesse Zaneveld","Samuel L. Peoples"]
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"

import unittest
import sys, os

from os.path import join
from os import listdir

testdir = os.path.dirname(__file__)
srcdir = '../karenina'
sys.path.insert(0, os.path.abspath(os.path.join(testdir, srcdir)))

from warnings import catch_warnings
from karenina.benchmark import benchmark_simulated_datasets,\
  graph_absolute_errors

"""
Tests for fit_timeseries_benchmark.py
"""

class TestFit(unittest.TestCase):
    """Tests of the fit timeseries benchmarking module"""

    def setUp(self):
        self.OutputDir = "./test_benchmark_results/"
        if not os.path.exists(self.OutputDir):
            os.makedirs(self.OutputDir)
        

    def test_benchmark_simulated_datasets_outputs_expected_files(self):
        """benchmark_simulated_datasets outputs expected files"""
        self.df = benchmark_simulated_datasets(3, niter=[1],\
          output_dir=self.OutputDir, verbose=False)

        expected_files = ['fit_timeseries_benchmark_log_3.txt',\
         "fit_timeseries_benchmark_3_results.csv"]
        observed_files = listdir(self.OutputDir)
        for exp_filename in expected_files:
            self.assertIn(exp_filename,observed_files)

    
    def test_benchmark_simulated_datasets_has_no_null_values(self):
        """benchmark_simulated datasets returns a dataframe with no null values"""
        benchmark_result_df =\
          benchmark_simulated_datasets(3, output_dir=self.OutputDir, verbose=False)          
        
        assert benchmark_result_df.isnull().any().any() == False

    #matplotlib.use('Agg') in test_experiment breaks this because backend is changed.

    def test_graph_absolute_errors_outputs_graphs(self):
        """graph_absolute_errors outputs image files"""
        
        benchmark_result_df =\
          benchmark_simulated_datasets(3, niter=1,\
          output_dir=self.OutputDir, verbose=False)  

        graph_absolute_errors(benchmark_result_df, self.OutputDir)
        
        expected_files = ["benchmark_theta_error.png",\
         "benchmark_sigma_error.png","benchmark_lambda_error.png"]

        observed_files = os.listdir(self.OutputDir)
        for exp_filename in expected_files:
            self.assertIn(exp_filename,observed_files)


    def tearDown(self):
        """Removes output files"""

        expected_filenames = ['fit_timeseries_benchmark_log_3.txt',\
         "fit_timeseries_benchmark_3_results.csv","benchmark_theta_error.png",\
         "benchmark_sigma_error.png","benchmark_lambda_error.png"]

        expected_filepaths = [join(self.OutputDir,f) for f in expected_filenames]
        
        for observed_file in listdir(self.OutputDir):
            if observed_file not in expected_filenames:
                raise ValueError("Found an unexpected file - {} not in {}".format(\
                  str(observed_file),str(expected_filenames)))                

        for exp_file in expected_filepaths:
            if(os.path.exists(exp_file)):
                os.remove(exp_file)
       
        os.rmdir(self.OutputDir)


if __name__ == '__main__':
    unittest.main()


