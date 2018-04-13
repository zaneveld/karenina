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
from individual import Individual
import numpy.testing as npt

"""
Tests for spatial_ornstein_uhlenbeck.py
"""


class TestIndividual(unittest.TestCase):
    # TODO: Tests

    def setUp(self):
        pass

    def test_apply_perturbation(self):
        pass

    def test_remove_perturbation(self):
        pass

    def test_remove_perturbation_from_axis(self):
        pass

    def test_apply_perturbation_to_axis(self):
        pass

    def test_check_identity(self):
        pass

    def test_simulate_movement(self):
        pass

    def test_get_data(self):
        pass

if __name__ == '__main__':
    unittest.main()


