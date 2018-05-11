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
from karenina.individual import Individual
from karenina.perturbation import Perturbation
import numpy.testing as npt
from copy import copy

"""
Tests for spatial_ornstein_uhlenbeck.py
"""


class TestExperiment(unittest.TestCase):
    # TODO: Tests

    def setUp(self):
        self.TreatmentNames = ['control','destabilizing_treatment']
        self.Treatments = [{'treatment_name': 'control'}, {'treatment_name': 'destabilizing_treatment'}]
        self.BaseParams = {'lambda': 0.2, 'delta': 0.25, 'interindividual_variation': 0.01}
        self.NIndividuals = [35,35]
        self.n_timepoints = 10
        self.treatment_params = [[],[]]
        self.interindividual_variation = 0.01

        self.exp = Experiment(self.TreatmentNames, self.NIndividuals, self.n_timepoints, self.BaseParams,
                         self.treatment_params, self.interindividual_variation)

    def test_check_variable_specified_per_treatment(self):
        Experiment.check_variable_specified_per_treatment(self.exp,
                                                          self.NIndividuals)
        self.TreatmentNames.append('Error')
        with self.assertRaises(ValueError):
            self.exp = Experiment(self.TreatmentNames, self.NIndividuals, self.n_timepoints, self.BaseParams,
                                  self.treatment_params, self.interindividual_variation)
            Experiment.check_variable_specified_per_treatment(self.exp,self.NIndividuals)

    def test_check_n_timepoints_is_int(self):
        Experiment.check_n_timepoints_is_int(self.exp,
                                             self.n_timepoints)
        self.n_timepoints = [10]
        with self.assertRaises(ValueError):
            self.exp = Experiment(self.TreatmentNames, self.NIndividuals, self.n_timepoints, self.BaseParams,
                                  self.treatment_params, self.interindividual_variation)
            Experiment.check_n_timepoints_is_int(self.exp, self.n_timepoints)

    def test_simulate_timesteps(self):
        pass

    def test_simulate_timestep(self):
        #Experiment.simulate_timestep(self, 0)
        pass

    def test_writeToMovieFile(self):
        # Travis-CI Uses Xwindows backend, this prevents that issue.
        import matplotlib
        matplotlib.use('Agg')

        self.output_folder = "./"
        Experiment.write_to_movie_file(self.exp,
                                       self.output_folder)
        assert os.path.exists("./simulation_video.mp4")
        os.remove("./simulation_video.mp4")

if __name__ == '__main__':
    unittest.main()


