import unittest
import numpy
import os
import mdarray as mt


def diff(a, b):
    ma = numpy.array(a)
    mb = numpy.array(b)
    return numpy.max(numpy.abs(a-b))


class TestTrajectoryMolden(unittest.TestCase):

    def setUp(self):

        self.testDir = os.path.dirname(os.path.realpath(__file__))


    def test_styleGeometries(self):

        fp = self.testDir + "/3wat.molden"
        traj = mt.Trajectory(fp)
        self.assertEqual(traj.nAtoms, 9)
        symbols = ['O', 'H', 'H' ] * 3
        self.assertEqual(traj.symbols, symbols)
        num = [ 8, 1, 1 ] * 3
        self.assertEqual(list(traj.aNumbers), num)

        frame = traj.read()
        self.assertTrue('coordinates' in frame)

        first = [ -0.000582, -1.638651, 0.057954]
        last = [ -0.498708, 1.159692, -0.064495]
        crd = frame['coordinates']
        self.assertTrue(diff(crd[0,:], first) < 1e-5)
        self.assertTrue(diff(crd[-1,:], last) < 1e-5)

        # Rewind to last frame
        nextFrame = traj.read()
        while nextFrame:
            frame = nextFrame
            nextFrame = traj.read()

        self.assertEqual(traj.lastFrame, 31)
        first = [ -0.063459, -1.525383, 0.106962]
        last = [ -0.509157, 1.190487, -0.171253]
        crd = frame['coordinates']
        self.assertTrue(diff(crd[0,:], first) < 1e-5)
        self.assertTrue(diff(crd[-1,:], last) < 1e-5)


    def test_styleAtoms(self):

        fp = self.testDir + "/bf3nh3.molden"
        traj = mt.Trajectory(fp)
        self.assertEqual(traj.nAtoms, 8)
        symbols = ['B', 'N', 'F', 'F', 'F', 'H', 'H', 'H' ]
        self.assertEqual(traj.symbols, symbols)
        num = [ 5, 7, 9, 9, 9, 1, 1, 1 ]
        self.assertEqual(list(traj.aNumbers), num)

        frame = traj.read()
        self.assertTrue('coordinates' in frame)

        first = [ 0.000000, 0.000000, -0.133009]
        last = [ -0.951454, 0.000000, 1.903657]
        crd = frame['coordinates']
        self.assertTrue(diff(crd[0,:], first) < 1e-5)
        self.assertTrue(diff(crd[-1,:], last) < 1e-5)

        # Try to read next frame
        nextFrame = traj.read()
        self.assertEqual(nextFrame, None)


    def test_styleFreq(self):

        fp = self.testDir + "/freq.molden"
        traj = mt.Trajectory(fp)
        self.assertEqual(traj.nAtoms, 8)
        symbols = ['B', 'N', 'F', 'F', 'F', 'H', 'H', 'H' ]
        self.assertEqual(traj.symbols, symbols)
        num = [ 5, 7, 9, 9, 9, 1, 1, 1 ]
        self.assertEqual(list(traj.aNumbers), num)

        frame = traj.read()
        self.assertTrue('coordinates' in frame)

        first = [ 0.000000, 0.000000, -0.133009]
        last = [ -0.951454, 0.000000, 1.903657]
        crd = frame['coordinates']
        self.assertTrue(diff(crd[0,:], first) < 1e-5)
        self.assertTrue(diff(crd[-1,:], last) < 1e-5)

        # Try to read next frame
        nextFrame = traj.read()
        self.assertEqual(nextFrame, None)
