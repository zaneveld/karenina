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
from process import Process
from perturbation import Perturbation
import numpy.testing as npt

"""
Tests for spatial_ornstein_uhlenbeck.py
"""

class TestProcess(unittest.TestCase):
    """Tests of the Process class"""

    def setUp(self):
        self.TestProcesses = {}
        # Note that each Process is 1-D, so only one coordinate is provided

        # typical parameters
        start_coord = 0.0
        attractor_pos = 0.0
        history = None
        process_type = "Ornstein-Uhlenbeck"

        # specific parameters

        # First let's generate a stable process that has lambda = 1
        # This should mean that the process ALWAYS revers to its attractor
        # at every timestep.  So delta shouldn't matter.
        params = {"lambda": 1.0, "delta": 0.0, "mu": attractor_pos}
        stable_process = Process(start_coord, \
                                 motion=process_type, params=params)
        self.TestProcesses["stable_process"] = stable_process

    def test_stable_process_update(self):
        """A stable OU process (lambda =1) is invariable in position"""

        stable_process = self.TestProcesses["stable_process"]
        start_coord = stable_process.Coord
        a_long_time = 1000
        for i in range(a_long_time):
            stable_process.update(1.0)
        end_coord = stable_process.Coord
        npt.assert_almost_equal(start_coord, end_coord)

    def test_perturb_stable_process(self):
        """Process can be perturbed to become more variable"""
        stable_process = self.TestProcesses["stable_process"]
        start_coord = stable_process.Coord

        a_long_time = 1000

        for i in range(a_long_time):
            stable_process.update(1.0)
        unperturbed_coord = stable_process.Coord

        # Start and end shouldn't matter as these are just
        # stored in 'dumb' variables accessed by external objects
        perturbation = Perturbation(start=0, end=0, params={"delta": 1.0, "lambda": 0.0}, update_mode="replace")
        stable_process.Perturbations.append(perturbation)
        for i in range(a_long_time):
            stable_process.update(1.0)
        end_coord = stable_process.Coord

        # We will test that Perturbations actually update
        # process parameters during simulation elsewhere

        # Here I just want to confirm that the unstable process
        # which becomes pure Brownian motion
        # diverges more over time than the stable process
        stable_change = abs(unperturbed_coord - start_coord)

        unstable_change = abs(end_coord - unperturbed_coord)

        self.assertTrue(stable_change < unstable_change)

    def test_perturb_alter_mean(self):
        """Process can be perturbed to take on a new mean location"""
        stable_process = self.TestProcesses["stable_process"]
        start_coord = stable_process.Coord

        a_long_time = 1000

        for i in range(a_long_time):
            stable_process.update(1.0)
        unperturbed_coord = stable_process.Coord

        new_mu = -0.5
        perturbation = Perturbation(start=0, end=0, params={"mu": new_mu}, update_mode="replace")
        stable_process.Perturbations.append(perturbation)
        for i in range(a_long_time):
            stable_process.update(1.0)
        end_coord = stable_process.Coord

        # We will test that Perturbations actually update
        # process parameters during simulation elsewhere

        # Here I just want to confirm that after running for a while
        # altering mu moves the process closer to its new mean
        stable_diff = abs(unperturbed_coord - new_mu)
        unstable_diff = abs(end_coord - new_mu)

        self.assertTrue(unstable_diff < stable_diff)


if __name__ == '__main__':
    unittest.main()