import mdarray
import unittest
import os

class TestFindBonds(unittest.TestCase):

    def setUp(self):

        self.testDir = os.path.dirname(os.path.realpath(__file__))

    def test_findBonds(self):

        traj = mdarray.Trajectory(self.testDir+'/3mols.xyz')
        frame = traj.read()
        bonds = mdarray.findBonds(traj.symbols, frame['coordinates'])
        expected = [ (0, 1), (0, 2), (0, 3), (0, 4), (5, 6), (5, 7) ]
        self.assertEqual(bonds, expected)

        bonds = mdarray.findBonds(traj.symbols,
            frame['coordinates'], format="atoms")

        expected = [ [1,2,3,4], [0], [0], [0], [0], [6,7], [5], [5], [] ]
        self.assertEqual(bonds, expected)
