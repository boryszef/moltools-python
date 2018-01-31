import mdarray
import unittest
import os

class TestTopology(unittest.TestCase):

    def setUp(self):

        systems = [ '3mols.xyz', 'pyr.xyz', ]#'nh3h2so4.xyz' ]

        self.expected = {
            '3mols.xyz' : [
                [ (0, 1), (0, 2), (0, 3), (0, 4), (5, 6), (5, 7) ],
                [ [1,2,3,4], [0], [0], [0], [0], [6,7], [5], [5], [] ] ],
            'pyr.xyz' : [
                [ (0, 1), (0, 4), (0, 9), (1, 2), (1, 8), (2, 3), (2, 7),
                  (3, 4), (3, 6), (4, 5) ],
                [ [1,4,9], [0,2,8], [1,3,7], [2,4,6], [0,3,5], [4], [3],
                  [2], [1], [0] ] ]
        }

        self.testDir = os.path.dirname(os.path.realpath(__file__))

        self.systems = {}
        for s in systems:
            traj = mdarray.Trajectory(self.testDir+'/'+s)
            frame = traj.read()
            sys = {}
            sys['natoms'] = traj.nAtoms
            sys['symbols'] = traj.symbols
            sys['crd'] = frame['coordinates']
            self.systems[s] = sys

    def test_findBonds(self):

        for k, v in self.systems.items():
            bonds = mdarray.findBonds(v['symbols'], v['crd'])
            print(k, bonds)
            self.assertEqual(bonds, self.expected[k][0])

            bonds = mdarray.findBonds(v['symbols'],
                v['crd'], format="atoms")
            print(k, bonds)
            self.assertEqual(bonds, self.expected[k][1])

    def test_findMolecules(self):

        #bonds = mdarray.findBonds(self.symbols, self.crd)
        #mols = mdarray.findMolecules(self.natoms, bonds)

        #self.assertEqual(len(mols), 3)
        #s = [ set(x) for x in mols ]
        #self.assertTrue(set((0,1,2,3,4)) in s)
        #self.assertTrue(set((5,6,7)) in s)
        #self.assertTrue(set((8,)) in s)
        pass
