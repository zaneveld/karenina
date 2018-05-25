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
testdir = os.path.dirname(__file__)
srcdir = '../karenina'
sys.path.insert(0, os.path.abspath(os.path.join(testdir, srcdir)))
from karenina.perturbation import Perturbation
import numpy.testing as npt

"""
Tests for spatial_ornstein_uhlenbeck.py
"""

class TestPerturbation(unittest.TestCase):
    """Tests of the Perturbation class"""

    def setUp(self):
        # Set up some common perturbations to test

        perturbation_start = 50
        perturbation_end = 75

        set_low_lambda = Perturbation(perturbation_start, perturbation_end, params={"lambda": 0.005},
                                      update_mode="replace")
        set_zero_lambda = Perturbation(perturbation_start, perturbation_end, params={"lambda": 0.0},
                                       update_mode="replace")
        double_lambda = Perturbation(perturbation_start, perturbation_end, params={"lambda": 2.0},
                                     update_mode="multiply")
        double_lambda_and_mu = Perturbation(perturbation_start, perturbation_end, params={"lambda": 2.0, "mu": 2.0},
                                            update_mode="multiply")
        halve_lambda = Perturbation(perturbation_start, perturbation_end, params={"lambda": 0.5},
                                    update_mode="multiply")
        add_to_lambda = Perturbation(perturbation_start, perturbation_end, params={"lambda": 0.1}, update_mode="add")
        subtract_from_lambda = Perturbation(perturbation_start, perturbation_end, params={"lambda": 0.1},
                                            update_mode="add")
        set_high_lambda_low_delta = Perturbation(perturbation_start, perturbation_end,
                                                 params={"lambda": 0.5, "delta": 0.1}, update_mode="replace")
        self.TestPerturbations = {"set_low_lambda": set_low_lambda, "set_zero_lambda": set_zero_lambda,
                                  "double_lambda": double_lambda, \
                                  "halve_lambda": halve_lambda, "add_to_lambda": add_to_lambda,
                                  "subtract_from_lambda": subtract_from_lambda, \
                                  "set_high_lambda_low_delta": set_high_lambda_low_delta,
                                  "double_lambda_and_mu": double_lambda_and_mu}

    def test_is_active(self):
        """Perturbation is_active tests whether it is active based on the current time"""
        # perturbations start at 50 and end at 75
        exact_start = 50
        exact_end = 75
        in_the_middle = 55
        start_of_experiment = 0

        curr_perturbation = self.TestPerturbations["set_zero_lambda"]
        self.assertTrue(curr_perturbation.is_active(in_the_middle))
        self.assertFalse(curr_perturbation.is_active(start_of_experiment))

        # let's test the borderline cases, which should be active
        self.assertTrue(curr_perturbation.is_active(exact_start))
        self.assertTrue(curr_perturbation.is_active(exact_end))

        # Test that if we change Start and End we change
        # is_active results

        curr_perturbation.Start = 90
        curr_perturbation.End = 100

        self.assertFalse(curr_perturbation.is_active(55))
        self.assertTrue(curr_perturbation.is_active(95))

    def test_update_param_updates_single_param(self):
        """Perturbation.update_params updates a single parameter correctly"""

        base_params = {"mu": 0.1, "lambda": 0.25, "delta": 0.18}

        curr_perturbation = self.TestPerturbations["set_low_lambda"]
        obs = curr_perturbation.update_params(base_params)
        exp = {"mu": 0.1, "lambda": 0.005, "delta": 0.18}
        self.assertEqual(obs, exp)

        curr_perturbation = self.TestPerturbations["set_zero_lambda"]
        obs = curr_perturbation.update_params(base_params)
        exp = {"mu": 0.1, "lambda": 0.0, "delta": 0.18}
        self.assertEqual(obs, exp)

        curr_perturbation = self.TestPerturbations["double_lambda"]
        obs = curr_perturbation.update_params(base_params)
        exp = {"mu": 0.1, "lambda": 0.50, "delta": 0.18}
        self.assertEqual(obs, exp)

        curr_perturbation = self.TestPerturbations["halve_lambda"]
        obs = curr_perturbation.update_params(base_params)
        exp = {"mu": 0.1, "lambda": 0.125, "delta": 0.18}
        self.assertEqual(obs, exp)

    def test_update_param_updates_multiple_params(self):
        """Perturbation.update_params updates a single parameter correctly"""

        base_params = {"mu": 0.1, "lambda": 0.25, "delta": 0.18}

        curr_perturbation = self.TestPerturbations["set_high_lambda_low_delta"]
        obs = curr_perturbation.update_params(base_params)
        exp = {"mu": 0.1, "lambda": 0.5, "delta": 0.1}
        self.assertEqual(obs, exp)

        curr_perturbation = self.TestPerturbations["double_lambda_and_mu"]
        obs = curr_perturbation.update_params(base_params)
        exp = {"mu": 0.2, "lambda": 0.5, "delta": 0.18}
        self.assertEqual(obs, exp)


if __name__ == '__main__':
    unittest.main()
