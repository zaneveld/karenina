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

from scipy.stats import norm
from copy import copy

class Process(object):
    """
    Represents a 1d process in a Euclidean space
    """

    def __init__(self,start_coord, motion = "Ornstein-Uhlenbeck",\
        history = None,params={"L":0.20,"delta":0.25}):
        """
        start_coords - float starting coordinate for the particle
        """
        if history is None:
            history = []
        self.StartCoord = start_coord
        self.Coord = start_coord
        self.History = history
        self.History.append(start_coord)
        self.Params = params
        self.ProcessType = motion
        self.Perturbations = []

    def update(self,dt):
        """
        Update the process

        :param dt: time elapsed since last update
        """
        curr_params = copy(self.Params)
        for p in self.Perturbations:
            curr_params = p.update_params(curr_params)
        
        if self.ProcessType == "Brownian":
            self.bm_update(dt,delta=curr_params["delta"])
        elif self.ProcessType == "Ornstein-Uhlenbeck":
            self.ou_update(dt,mu=curr_params["mu"],\
            delta = curr_params["delta"],\
            L=curr_params["lambda"])
        else:
            raise ValueError("only 'Brownian' or 'Ornstein-Uhlenbeck' processes are supported, not %s"%self.ProcessType)

    def bm_change(self,dt,delta):
        """
        Change the Brownian motion process

        :param dt: time elapsed since last update
        :param delta: delta value to adjust
        :return: change
        """
        change =  norm.rvs(loc=0,size=1,scale=delta**2*dt)
        return change

    def bm_update(self,dt,delta):
        """
        Update the Brownian motion process

        :param dt: time elapsed since last update
        :param delta: delta value to adjust
        :return: change
        """
        curr_coord = self.Coord
        self.History.append(curr_coord)
        change = self.bm_change(dt,delta)
        self.Coord = curr_coord + change


    def ou_change(self,dt,mu,L,delta):
        """
        Change the Ornstein Uhlenbeck motion process

        The Ornstein Uhlenbeck process is modelled as:

        ds = lambda * (mu - s) * dt + dW

        ds -- change in our process from the last timepoint
        L -- lambda, the speed of reversion to mean
        (NOTE: lambda is a reserved keyword in Python so I use L)
        mu -- mean position
        s -- current position
        dt -- how much time has elapsed since last update
        dW -- the Weiner Process (basic Brownian motion)

        This says we update as usual for Brownian motion,
        but add in a term that reverts us to some
        mean position (mu) over time (dt) at some speed (lambda)

        :param dt: time elapsed since last update
        :param delta: delta value to adjust
        :return: change in process since last timepoint
        """

        dW = self.bm_change(dt=dt,delta=delta)
        ds = L * (mu - self.Coord) * dt + dW
        return ds

    def ou_update(self,dt,mu,\
        L,delta,min_bound=-1.0,max_bound=1.0):
        """
        Update the Brownian motion process

        :param dt: time elapsed since last update
        :param mu: mean position
        :param L: speed of reversion to the mean
        :param delta: variance from the mean
        :param min_bound: minimum bound
        :param max_bound: maximum bound
        """
        curr_coord = self.Coord
        self.History.append(self.Coord)
        change = self.ou_change(dt=dt,mu=mu,L=L,delta=delta)
        self.Coord = curr_coord + change
        if min_bound is not None:
            self.Coord = max(self.Coord,min_bound)
        if max_bound is not None:
            self.Coord = min(self.Coord,max_bound)
