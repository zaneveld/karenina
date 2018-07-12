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


class TestExperiment(unittest.TestCase):
    """Tests for the Experiment class"""
    def setUp(self):
        """
        Creates default local variables for use in the tests
        """
        self.TreatmentNames = ['control','destabilizing_treatment']
        self.Treatments = [{'treatment_name': 'control'}, {'treatment_name': 'destabilizing_treatment'}]
        self.BaseParams = {'lambda': 0.2, 'delta': 0.25, 'interindividual_variation': 0.01}
        self.NIndividuals = [35,35]
        self.n_timepoints = 10
        self.treatment_params = [[],[]]
        self.interindividual_variation = 0.01
        self.verbose = False

        self.exp = Experiment(self.TreatmentNames, self.NIndividuals, self.n_timepoints, self.BaseParams,
                         self.treatment_params, self.interindividual_variation, self.verbose)


    def test_check_variable_specified_per_treatment(self):
        """Tests that the NIndividuals is equal to the number of Treatment Names."""
        self.exp.check_variable_specified_per_treatment(self.NIndividuals, self.verbose)
        self.TreatmentNames.append('Error')
        with self.assertRaises(ValueError):
            self.exp = Experiment(self.TreatmentNames, self.NIndividuals, self.n_timepoints, self.BaseParams,
                                  self.treatment_params, self.interindividual_variation, self.verbose)
            self.exp.check_variable_specified_per_treatment(self.NIndividuals, self.verbose)


    def test_check_n_timepoints_is_int(self):
        """Tests that the n_timepoints is of the int datatype."""
        self.exp.check_n_timepoints_is_int(self.n_timepoints)
        self.n_timepoints = [10]
        with self.assertRaises(ValueError):
            self.exp = Experiment(self.TreatmentNames, self.NIndividuals, self.n_timepoints, self.BaseParams,
                                  self.treatment_params, self.interindividual_variation, self.verbose)
            self.exp.check_n_timepoints_is_int(self.n_timepoints)


    def test_simulate_timesteps(self):
        #simulate_timestep only implements simulate_movement
        #   expected value could be tested there.

        """Tests that the timesteps are successfully completed, populating Data with expected number of entries."""
        assert len(self.exp.Data) == 1
        self.exp.simulate_timesteps(0,self.n_timepoints, self.verbose)
        assert len(self.exp.Data) == 701

    # Travis CI does not like using matplotlib. alternate solution required to test this.

    def test_writeToMovieFile(self):
        """TODO"""
        """Tests that the output movie file is successfully written, then removes the file."""
        pass
    """
        # Travis-CI Uses Xwindows backend, this prevents that issue.
        import os
        import matplotlib as mpl
        if os.environ.get('DISPLAY', '') == '':
            mpl.use('Agg')

        self.output_folder = "./"
        Experiment.write_to_movie_file(self.exp,
                                       self.output_folder, self.verbose)
        assert os.path.exists("./simulation_video.mp4")
        os.remove("./simulation_video.mp4")
    """

if __name__ == '__main__':
    unittest.main()


