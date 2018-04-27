#!/usr/bin/env python

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2011-2013, The PICRUSt Project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"

import sys, os
testdir = os.path.dirname(__file__)
srcdir = '../karenina'
sys.path.insert(0, os.path.abspath(os.path.join(testdir, srcdir)))
import unittest
from warnings import catch_warnings
from karenina.experiment import Experiment
import numpy.testing as npt

"""
Tests for spatial_ornstein_uhlenbeck.py
"""


class TestExperiment(unittest.TestCase):
    # TODO: Tests

    def setUp(self):
        pass

    def test_check_variable_specified_per_treatment(self):
        pass

    def test_check_n_timepoints_is_int(self):
        pass

    def test_simulate_timesteps(self):
        pass

    def test_simulate_timestep(self):
        pass

    def test_writeToMovieFile(self):
        pass

if __name__ == '__main__':
    unittest.main()


