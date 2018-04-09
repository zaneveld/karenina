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
from warnings import catch_warnings
import numpy.testing as npt
import sys, os
testdir = os.path.dirname(__file__)
srcdir = '../karenina'
sys.path.insert(0, os.path.abspath(os.path.join(testdir, srcdir)))
from visualization import get_timeseries_data, save_simulation_figure, save_simulation_movie, update_3d_plot

"""
Tests for spatial_ornstein_uhlenbeck.py
"""


class TestExperiment(unittest.TestCase):
    # TODO: Tests

    def setUp(self):
        pass

    def test_get_timeseries_data(self):
        pass

    def test_save_simulation_figure(self):
        pass

    def test_save_simulation_movie(self):
        pass

    def test_update_3d_plot(self):
        pass

if __name__ == '__main__':
    unittest.main()


