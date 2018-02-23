import unittest
import numpy
import os
import mdarray as mt


class TestTrajectoryXTC(unittest.TestCase):

    def setUp(self):

        if not mt.__config__['gromacs']:
            self.skipTest('Skipping tests since gromacs support was not compiled')
        self.testDir = os.path.dirname(os.path.realpath(__file__))


    def test_readXTC(self):

        fp = self.testDir + "/traj.xtc"
        traj = mt.Trajectory(fp)
        self.assertEqual(traj.nAtoms, 10)

        frame = traj.read()
        self.assertTrue('coordinates' in frame)
        self.assertTrue('step' in frame)
        self.assertTrue('time' in frame)
        self.assertTrue('box' in frame)

        first = [ 9.9900, 13.8300, 15.0400]
        last = [ 10.3500, 14.5100, 12.1400]
        crd = frame['coordinates']
        self.assertAlmostEqual(crd[0,0], first[0], places=3)
        self.assertAlmostEqual(crd[0,1], first[1], places=3)
        self.assertAlmostEqual(crd[0,2], first[2], places=3)
        self.assertAlmostEqual(crd[-1,0], last[0], places=3)
        self.assertAlmostEqual(crd[-1,1], last[1], places=3)
        self.assertAlmostEqual(crd[-1,2], last[2], places=3)
        self.assertEqual(frame['step'], 0)
        self.assertAlmostEqual(frame['time'], 0.0000, places=3)
        box = [35.091, 34.575, 37.084]
        self.assertAlmostEqual(frame['box'][0,0], box[0], places=3)
        self.assertAlmostEqual(frame['box'][1,1], box[1], places=3)
        self.assertAlmostEqual(frame['box'][2,2], box[2], places=3)

        # Rewind to last frame
        nextFrame = traj.read()
        while nextFrame:
            frame = nextFrame
            nextFrame = traj.read()

        self.assertEqual(traj.lastFrame, 25)
        first = [ 9.9200, 14.7100, 13.3400]
        last = [ 11.1800, 13.6000, 16.0000]
        crd = frame['coordinates']
        self.assertAlmostEqual(crd[0,0], first[0], places=3)
        self.assertAlmostEqual(crd[0,1], first[1], places=3)
        self.assertAlmostEqual(crd[0,2], first[2], places=3)
        self.assertAlmostEqual(crd[-1,0], last[0], places=3)
        self.assertAlmostEqual(crd[-1,1], last[1], places=3)
        self.assertAlmostEqual(crd[-1,2], last[2], places=3)
        self.assertEqual(frame['step'], 125000)
        self.assertAlmostEqual(frame['time'], 250.0, places=3)
        self.assertAlmostEqual(frame['box'][0,0], box[0], places=3)
        self.assertAlmostEqual(frame['box'][1,1], box[1], places=3)
        self.assertAlmostEqual(frame['box'][2,2], box[2], places=3)
