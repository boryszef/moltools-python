import unittest
import numpy
import mdarray as mt


class TestTrajectoryXTC(unittest.TestCase):

    def setUp(self):
        if not mt.__config__['gromacs']:
            self.skipTest('Skipping tests since gromacs support was not compiled')
        
