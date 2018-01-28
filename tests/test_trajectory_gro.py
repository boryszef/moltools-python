import unittest
import tempfile
import os
import numpy
import mdarray as mt

DATA = """GROningen Mixture of Alchemy and Childrens' Stories
   23
    1ALA      N    1   0.000   0.000   0.000
    1ALA     CA    2   0.000   0.000   0.146
    1ALA      C    3   0.140   0.000   0.202
    1ALA      O    4   0.176   0.075   0.291
    1ALA     CB    5  -0.078  -0.121   0.199
    1ALA     1H    6   0.042  -0.076  -0.052
    1ALA     2H    7  -0.097   0.002  -0.031
    1ALA     3H    8   0.050   0.083  -0.031
    1ALA     HA    9  -0.050   0.092   0.179
    1ALA    1HB   10  -0.180  -0.121   0.159
    1ALA    2HB   11  -0.083  -0.117   0.309
    1ALA    3HB   12  -0.029  -0.215   0.169
    2ALA      N   13   0.220  -0.090   0.142
    2ALA     CA   14   0.359  -0.106   0.182
    2ALA      C   15   0.436   0.023   0.167
    2ALA      O   16   0.509   0.067   0.254
    2ALA     CB   17   0.426  -0.217   0.100
    2ALA    OXT   18   0.427   0.090   0.065
    2ALA      H   19   0.179  -0.148   0.068
    2ALA     HA   20   0.360  -0.134   0.288
    2ALA    1HB   21   0.370  -0.311   0.111
    2ALA    2HB   22   0.529  -0.233   0.134
    2ALA    3HB   23   0.428  -0.190  -0.007
   1.00000   1.50000   0.50000
"""


DATAV = """GROningen Mixture of Alchemy and Childrens' Stories
   23
    1ALA      N    1   0.000   0.000   0.000  0.0000  0.0000  0.0000
    1ALA     CA    2   0.000   0.000   0.146  0.0060  0.0006  0.1466
    1ALA      C    3   0.140   0.000   0.202  0.1420  0.0002  0.2022
    1ALA      O    4   0.176   0.075   0.291  0.1716  0.0751  0.2911
    1ALA     CB    5  -0.078  -0.121   0.199 -0.0798 -0.1219  0.1999
    1ALA     1H    6   0.042  -0.076  -0.052  0.0422 -0.0762 -0.0522
    1ALA     2H    7  -0.097   0.002  -0.031 -0.0917  0.0021 -0.0311
    1ALA     3H    8   0.050   0.083  -0.031  0.0510  0.0831 -0.0311
    1ALA     HA    9  -0.050   0.092   0.179 -0.0590  0.0929  0.1799
    1ALA    1HB   10  -0.180  -0.121   0.159 -0.1890 -0.1219  0.1599
    1ALA    2HB   11  -0.083  -0.117   0.309 -0.0893 -0.1179  0.3099
    1ALA    3HB   12  -0.029  -0.215   0.169 -0.0299 -0.2159  0.1699
    2ALA      N   13   0.220  -0.090   0.142  0.2220 -0.0902  0.1422
    2ALA     CA   14   0.359  -0.106   0.182  0.3529 -0.1062  0.1822
    2ALA      C   15   0.436   0.023   0.167  0.4376  0.0237  0.1677
    2ALA      O   16   0.509   0.067   0.254  0.5049  0.0674  0.2544
    2ALA     CB   17   0.426  -0.217   0.100  0.4206 -0.2170  0.1000
    2ALA    OXT   18   0.427   0.090   0.065  0.4257  0.0905  0.0655
    2ALA      H   19   0.179  -0.148   0.068  0.1789 -0.1488  0.0688
    2ALA     HA   20   0.360  -0.134   0.288  0.3680 -0.1348  0.2888
    2ALA    1HB   21   0.370  -0.311   0.111  0.3710 -0.3111  0.1111
    2ALA    2HB   22   0.529  -0.233   0.134  0.5249 -0.2334  0.1344
    2ALA    3HB   23   0.428  -0.190  -0.007  0.4278 -0.1907 -0.0077
   1.00000   1.50000   0.50000
"""


class TestTrajectoryGRO(unittest.TestCase):

    def setUp(self):
        
        # Create tempdir for XYZ files
        self.tmpDir = tempfile.mkdtemp()

        lines = DATAV.splitlines()
        self.nAtoms = int(lines[1])
        self.comment = lines[0].strip()
        resids = []
        self.resnames = []
        self.symbols = []
        crd = []
        vel = []
        for i in range(2, 2+self.nAtoms):
            l = lines[i]
            resids.append(int(l[:5]))
            self.resnames.append(l[5:10].strip())
            self.symbols.append(l[10:15].strip())
            crd.append((float(l[20:28]), float(l[28:36]), float(l[36:44])))
            vel.append((float(l[44:52]), float(l[52:60]), float(l[60:68])))
        self.resids = numpy.array(resids)
        self.crd = numpy.array(crd)
        self.crd *= 10
        self.vel = numpy.array(vel)
        l = lines[-1].split()
        box = list(map(float, l))
        self.box = numpy.zeros((3,3))
        self.box[0,0] = box[0]
        self.box[1,1] = box[1]
        self.box[2,2] = box[2]
        self.box *= 10
        self.fileName = "test.gro"


    def tearDown(self):

        # Remove files
        for f in os.listdir(self.tmpDir):
            if not f.startswith("."): os.remove(self.tmpDir+"/"+f)

        # Remove directory
        os.rmdir(self.tmpDir)


    def test_write(self):

        full = "%s/%s" % (self.tmpDir, self.fileName)
        traj = mt.Trajectory(full, "w", self.symbols, self.resids, self.resnames)
        traj.write(self.crd, box=self.box, comment=self.comment)
        del traj
        with open(full) as f: saved = f.read()
        self.assertEqual(saved, DATA)
        os.remove(full)


    def test_writeVelocities(self):

        full = "%s/read.gro" % self.tmpDir
        traj = mt.Trajectory(full, "w", self.symbols, self.resids, self.resnames)
        traj.write(self.crd, self.vel, self.box, self.comment)
        del traj
        with open(full) as f: saved = f.read()
        self.assertEqual(saved, DATAV)
        os.remove(full)


    def test_read(self):

        full = "%s/read.gro" % self.tmpDir
        with open(full, 'w') as f:
            f.write(DATA)
        traj = mt.Trajectory(full)
        self.assertEqual(traj.nAtoms, self.nAtoms)
        self.assertTrue(numpy.all(traj.resids == self.resids))
        self.assertEqual(traj.resNames, self.resnames)
        self.assertEqual(traj.symbols, self.symbols)
        frame = traj.read()
        self.assertEqual(frame['comment'], self.comment)
        diff = self.crd - frame['coordinates']
        maxDiff = numpy.max(numpy.abs(diff))
        self.assertTrue(maxDiff <= 0.01)
        diff = self.box - frame['box']
        maxDiff = numpy.max(numpy.abs(diff))
        self.assertTrue(maxDiff <= 0.00001)
        os.remove(full)


    def test_readVelocities(self):

        full = "%s/read.gro" % self.tmpDir
        with open(full, 'w') as f:
            f.write(DATAV)
        traj = mt.Trajectory(full)
        frame = traj.read()
        diff = self.vel - frame['velocities']
        maxDiff = numpy.max(numpy.abs(diff))
        self.assertTrue(maxDiff <= 0.0001)
        os.remove(full)
