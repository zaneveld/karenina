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

from process import Process
from random import random,randint
from copy import copy

class Individual(object):
    def __init__(self,subject_id,coords=["x","y","z"],metadata={},params={},interindividual_variation=0.01):
        self.SubjectId = subject_id
        self.Metadata = metadata
        self.MovementProcesses = {}
        self.BaseParams = params
        for c in coords:
            #print "STARTING PROCESS for Axis:",c

            #simulate minor interindividual variation
            if c in params.keys():
                start_coord = params[c]
                start_mu = params[c]
            else:
                start_coord =  (random()-0.50) *\
                     2.0 * ( interindividual_variation)
                start_mu = start_coord
            #print "START COORD %s: %f" %(c,start_coord)

            #For OU processes, assume the process reverts
            # to its starting coordinate
            curr_params = copy(self.BaseParams)
            curr_params["mu"] = start_coord
            print ("start_coord:",start_coord)
            print ("curr_params['mu']",curr_params['mu'])
            self.MovementProcesses[c] = Process(start_coord = start_coord,params=curr_params,\
              motion = "Ornstein-Uhlenbeck")

    def apply_perturbation(self,perturbation):
        """Apply a perturbation to the appropriate axes"""
        for axis in perturbation.Axes:
            self.apply_perturbation_to_axis(axis,perturbation)

    def remove_perturbation(self,perturbation):
        for axis in perturbation.Axes:
            self.remove_perturbation_from_axis(axis,perturbation)

    def remove_perturbation_from_axis(self,axis,perturbation):
        """Remove a perturbation from one or more Process objects"""
        self.MovementProcesses[axis].Perturbations.remove(perturbation)

    def apply_perturbation_to_axis(self,axis,perturbation):
        """Apply a perturbation to a Processes objects"""

        self.MovementProcesses[axis].Perturbations.append(perturbation)

    def check_identity(self):
        for coord_name,movement_process in self.MovementProcesses.iteritems():
            for coord_name2,movement_process2 in self.MovementProcesses.iteritems():
                if coord_name == coord_name2:
                    continue
                if movement_process.History == movement_process2.History:
                    #print "History1:",movement_process.History
                    #print "History2:",movement_process2.History
                    return True
        return False

    def simulate_movement(self,n_timepoints,params=None):
        if not params:
            params = self.BaseParams

        for c in self.MovementProcesses.keys():
            for t in range(n_timepoints):
                mu = self.MovementProcesses[c].StartCoord
                self.MovementProcesses[c].update(dt=1.0)

    def get_data(self,n_timepoints):
        result = []
        coords = self.MovementProcesses.keys()
        for t in range(n_timepoints):
            timepoint = [self.MovementProcesses[c].History[t] for c in coords]
            curr_sampleid = "S%s_t%i"%(self.SubjectId,t)
            curr_data = [curr_sampleid]
            curr_data.extend(map(float,timepoint))
            result.append(curr_data)
        return result