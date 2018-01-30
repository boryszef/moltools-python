import mdarray
import unittest
import os

class TestTopology(unittest.TestCase):

    def setUp(self):

        self.testDir = os.path.dirname(os.path.realpath(__file__))
        traj = mdarray.Trajectory(self.testDir+'/3mols.xyz')
        frame = traj.read()
        self.natoms = traj.nAtoms
        self.symbols = traj.symbols
        self.crd = frame['coordinates']

    def test_findBonds(self):

        bonds = mdarray.findBonds(self.symbols, self.crd)
        expected = [ (0, 1), (0, 2), (0, 3), (0, 4), (5, 6), (5, 7) ]
        self.assertEqual(bonds, expected)

        bonds = mdarray.findBonds(self.symbols,
            self.crd, format="atoms")

        expected = [ [1,2,3,4], [0], [0], [0], [0], [6,7], [5], [5], [] ]
        self.assertEqual(bonds, expected)

    def test_findMolecules(self):

        bonds = mdarray.findBonds(self.symbols, self.crd)
        mols = mdarray.findMolecules(self.natoms, bonds)

        self.assertEqual(len(mols), 3)
        s = [ set(x) for x in mols ]
        self.assertTrue(set((0,1,2,3,4)) in s)
        self.assertTrue(set((5,6,7)) in s)
        self.assertTrue(set((8,)) in s)
