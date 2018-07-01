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
from karenina.fit_timeseries_benchmark import benchmark,vis

"""
Tests for fit_timeseries_benchmark.py
"""


class TestFit(unittest.TestCase):
    """Tests of the fit timeseries benchmarking module"""

    def setUp(self):
        if not os.path.exists("./test_benchmark/"):
            os.makedirs("./test_benchmark/")
        self.df = benchmark(3, output="./test_benchmark/", verbose=False)

    def test_benchmark(self):
        """Tests that the benchmarking output is complete and files are saved"""
        assert os.path.exists("./test_benchmark/fit_timeseries_benchmark3_log.txt")
        assert os.path.exists("./test_benchmark/fit_timeseries_benchmark3.csv")
        assert self.df.isnull().any().any() == False
        """TODO"""
        pass

    def test_vis(self):
        """Ensures that the output visualizations are generated"""
        vis(self.df, "./test_benchmark/")
        assert os.path.exists("./test_benchmark/benchmark_sigma_err.png")
        assert os.path.exists("./test_benchmark/benchmark_lambda_err.png")
        assert os.path.exists("./test_benchmark/benchmark_theta_err.png")

        """TODO"""
        pass

    def tearDown(self):
        """Removes output files"""
        if(os.path.exists("./test_benchmark/fit_timeseries_benchmark3_log.txt")):
            os.remove("./test_benchmark/fit_timeseries_benchmark3_log.txt")
        if (os.path.exists("./test_benchmark/fit_timeseries_benchmark3.csv")):
            os.remove("./test_benchmark/fit_timeseries_benchmark3.csv")
        if (os.path.exists("./test_benchmark/benchmark_sigma_err.png")):
            os.remove("./test_benchmark/benchmark_sigma_err.png")
        if (os.path.exists("./test_benchmark/benchmark_lambda_err.png")):
            os.remove("./test_benchmark/benchmark_lambda_err.png")
        if (os.path.exists("./test_benchmark/benchmark_theta_err.png")):
            os.remove("./test_benchmark/benchmark_theta_err.png")

        os.rmdir("./test_benchmark/")


if __name__ == '__main__':
    unittest.main()


