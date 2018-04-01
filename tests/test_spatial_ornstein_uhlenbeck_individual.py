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
from karenina.spatial_ornstein_uhlenbeck import Process, Individual, Experiment, Perturbation
import numpy.testing as npt

"""
Tests for spatial_ornstein_uhlenbeck.py
"""


class TestIndividual(unittest.TestCase):
    # TODO: Tests

    def setUp(self):
        pass

    def test_applyPerturbation(self):
        pass

    def test_removePerturbation(self):
        pass

    def test_removePerturbationFromAxis(self):
        pass

    def test_applyPerturbationToAxis(self):
        pass

    def test_check_identity(self):
        pass

    def test_simulate_movement(self):
        pass

    def test_get_data(self):
        pass

if __name__ == '__main__':
    unittest.main()


