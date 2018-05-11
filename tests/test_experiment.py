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

        # Prepare variables from init

        #for i,n in enumerate(n_individuals):
        for i,n in enumerate(self.NIndividuals):
            self.Treatments[i]["n_individuals"] = n

        self.NTimepoints = self.n_timepoints

        colors = ['fuchsia','cyan','darkorange','blue','yellow']
        #Set up the experimental subjects
        for treatment_idx,treatment in enumerate(self.Treatments):

            #Set a color for each treatment
            individuals = []
            params = copy(self.BaseParams)

            if treatment_idx < len(colors):
                params['color'] = colors[treatment_idx]
            else:
                params['color'] = 'lightgray'

            print(treatment)
            for i in range(treatment["n_individuals"]):

                curr_subject_id = "%s_%i" %(treatment["treatment_name"],i)
                curr_subject = Individual(subject_id = curr_subject_id,
                  params = params,\
                  metadata={"treatment":treatment["treatment_name"]},\
                  interindividual_variation=self.interindividual_variation)
                individuals.append(curr_subject)
            treatment["individuals"] = individuals


        #Set up the treatment parameters
        for treatment_idx,treatment in enumerate(self.Treatments):
            treatment["perturbations"] = []
            treatment["active_perturbations"] = []
            raw_perturbation_info = self.treatment_params[treatment_idx]
            #We should have a dict with start, end, and parms for each perturbation
            print ("raw_perturbation_info:",raw_perturbation_info)
            for p in raw_perturbation_info:
                print ("params:",p)
                curr_perturbation = Perturbation(p["start"],p["end"],p["params"],p["update_mode"],p["axes"])
                treatment["perturbations"].append(curr_perturbation)

        #Set up a place to hold data on the experiment outcome

        coords = ["x","y","z"]
        headers = "\t".join(["SampleID"]+coords)+"\n"
        self.Data = [headers]

    def test_check_variable_specified_per_treatment(self):
        Experiment.check_variable_specified_per_treatment(self,self.NIndividuals)
        self.TreatmentNames.append('Error')
        self.assertRaises(ValueError, Experiment.check_variable_specified_per_treatment,self,self.NIndividuals)

    def test_check_n_timepoints_is_int(self):
        Experiment.check_n_timepoints_is_int(self,self.n_timepoints)
        self.n_timepoints = [10]
        self.assertRaises(ValueError, Experiment.check_n_timepoints_is_int, self, self.n_timepoints)

    def test_simulate_timesteps(self):
        pass

    def test_simulate_timestep(self):
        #Experiment.simulate_timestep(self, 0)
        pass

    def test_writeToMovieFile(self):
        self.output_folder = "./"
        Experiment.write_to_movie_file(self, self.output_folder)
        assert os.path.exists("./simulation_video.mp4")
        os.remove("./simulation_video.mp4")

if __name__ == '__main__':
    unittest.main()


