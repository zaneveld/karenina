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
import sys, os
from os.path import join,realpath,dirname
testdir = os.path.dirname(__file__)
srcdir = '../karenina'
sys.path.insert(0, os.path.abspath(os.path.join(testdir, srcdir)))
from karenina.spatial_ornstein_uhlenbeck import check_perturbation_timepoint, write_options_to_log, parse_perturbation_file
import numpy.testing as npt
from pkg_resources import resource_filename

"""
Tests for spatial_ornstein_uhlenbeck.py
"""

class TestPrimary(unittest.TestCase):
    # TODO: Tests

    def setUp(self):
        self.verbose = False
        self.perturbation_file_path = os.path.abspath(resource_filename('karenina.data', 'set_x_lambda_small.tsv'))
        self.perturbation_timepoint = 5
        self.perturbation_duration = 100

    def test_check_perturbation_timepoint(self):
        """TODO"""
        pass

    def test_write_options_to_log(self):
        """TODO"""
        pass

    def test_parse_perturbation_file(self):
        """parse_perturbation_file functions with valid input"""

        obs_perturbations_list = parse_perturbation_file(self.perturbation_file_path, self.perturbation_timepoint,self.perturbation_duration)

        exp_perturbations_list = [{'start': 5, 'end': 105,\
        'params': {'lambda': 0.008, 'mu': -0.08},\
        'update_mode': 'replace', 'axes': ['x']}]

        self.assertEqual(obs_perturbations_list, exp_perturbations_list)


if __name__ == '__main__':
    unittest.main()
