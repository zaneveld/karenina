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
from spatial_ornstein_uhlenbeck  import Process,Individual,Experiment

  
"""
Tests for spatial_ornstein_uhlenbeck.py
"""



class TestProcess(unittest.TestCase):
    """Tests of the Process class"""

    def setUp(self):
        pass 


