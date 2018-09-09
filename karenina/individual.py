#/usr/bin/env python

from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2016, The Karenina Project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "0.0.1-dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"

from karenina.process import Process
from random import random,randint
from copy import copy

class Individual(object):
    """
    Generates an individual for OU simulation
    """
    def __init__(self,subject_id,coords=["x","y","z"],metadata={},params={},interindividual_variation=0.01, verbose=False):
        self.SubjectId = subject_id
        self.Metadata = metadata
        self.MovementProcesses = {}
        self.BaseParams = params
        self.verbose = verbose
        for c in coords:
            if verbose:
                print("STARTING PROCESS for Axis:",c)

            #simulate minor interindividual variation
            if c in params.keys():
                start_coord = params[c]
                start_mu = params[c]
            else:
                start_coord =  (random()-0.50) *\
                     2.0 * ( interindividual_variation)
                start_mu = start_coord
            if verbose:
                print("START COORD %s: %f" %(c,start_coord))

            #For OU processes, assume the process reverts
            # to its starting coordinate
            curr_params = copy(self.BaseParams)
            curr_params["mu"] = start_coord
            if verbose:
                print ("start_coord:",start_coord)
                print ("curr_params['mu']",curr_params['mu'])
            self.MovementProcesses[c] = Process(start_coord = start_coord,params=curr_params,\
              motion = "Ornstein-Uhlenbeck")

    def apply_perturbation(self,perturbation):
        """
        Apply a perturbation to the appropriate axes

        :param perturbation: perturbation to apply
        """
        for axis in perturbation.Axes:
            self.apply_perturbation_to_axis(axis,perturbation)

    def remove_perturbation(self,perturbation):
        """
        Remove a perturbation from the appropriate axes

        :param perturbation: perturbation to remove
        """
        for axis in perturbation.Axes:
            self.remove_perturbation_from_axis(axis,perturbation)

    def remove_perturbation_from_axis(self,axis,perturbation):
        """
        Remove a perturbation from one or more Process objects

        :param axis: axis to remove perturbation from
        :param perturbation: perturbation to remove from axis
        """
        self.MovementProcesses[axis].Perturbations.remove(perturbation)

    def apply_perturbation_to_axis(self,axis,perturbation):
        """
        Apply a perturbation to a Processes objects

        :param axis: Axis to apply perturbation to
        :param perturbation: perturbation to apply
        """
        self.MovementProcesses[axis].Perturbations.append(perturbation)

    def check_identity(self, verbose=False):
        """
        Check identity of movement process

        :param verbose: verbose output, default = False
        :return: True if processes are equivalent, False if not
        """
        for coord_name,movement_process in self.MovementProcesses.iteritems():
            for coord_name2,movement_process2 in self.MovementProcesses.iteritems():
                if coord_name == coord_name2:
                    continue
                if movement_process.History == movement_process2.History:
                    if verbose:
                        print("History1:",movement_process.History)
                        print("History2:",movement_process2.History)
                    return True
        return False

    def simulate_movement(self,n_timepoints,params=None):
        """
        Simulate movement over timepoints

        :param n_timepoints: number of timepoints to simulate
        :param params: parameters to change from baseparams
        """
        if not params:
            params = self.BaseParams

        for c in self.MovementProcesses.keys():
            for t in range(n_timepoints):
                mu = self.MovementProcesses[c].StartCoord
                self.MovementProcesses[c].update(dt=1.0)

    def get_data(self,n_timepoints):
        """
        get data from movement processes

        :param n_timepoints: number of timepoints to gather data for
        :return: data for timepoints
        """
        result = []
        coords = self.MovementProcesses.keys()
        for t in range(n_timepoints):
            timepoint = [self.MovementProcesses[c].History[t] for c in coords]
            curr_sampleid = "S%s_t%i"%(self.SubjectId,t)
            curr_data = [curr_sampleid]
            curr_data.extend(map(float,timepoint))
            result.append(curr_data)
        return result
