import unittest
import tempfile
import random
import os
import stat
import numpy
import mdarray as mt

# Return random string of spaces and tabs
def rndb():
    blanks = [' ', '\t']
    l = [ random.choice(blanks) for i in range(random.randint(0, 5)) ]
    random.shuffle(l)
    return "".join(l)

atomicMasses = {
    'H':(1,1.008), 'C':(6,12.011), 'N':(7,14.007), 'O':(8,15.999),
    'P':(15,30.973762), 'S':(16,32.06), 'Cl':(17,35.45), 'Na':(11,22.98976928),
    'Cu':(29,63.546), 'Fe':(26,55.845) }
characters = "".join([chr(x) for x in range(ord('A'), ord('z')+1)])

class TestTrajectoryXYZ(unittest.TestCase):

    def setUp(self):
        
        # Create tempdir for XYZ files
        self.tmpDir = tempfile.mkdtemp()

        # Keep for comparison
        self.data = []

        self.nFiles = 10
        for i in range(self.nFiles):

            structure = {}

            nAtoms = random.randint(1, 20)
            structure['nAtoms'] = nAtoms

            nFrames = random.randint(1, 10)
            structure['nFrames'] = nFrames

            xyz = open("%s/%d.xyz" % (self.tmpDir, i), "w")
            asym = list(atomicMasses.keys())
            symbols = [random.choice(asym) for x in range(nAtoms)]
            structure['symbols'] = symbols[:]

            structure['coordinates'] = []

            structure['comments'] = []

            for f in range(nFrames):
                xyz.write("%s%d%s\n" % (rndb(), nAtoms, rndb()))
                comLen = random.randint(0, 100)
                comment = "".join([random.choice(characters) for x in range(comLen)])
                structure['comments'].append(comment[:])
                xyz.write("%s\n" % comment)
                crd = numpy.random.uniform(-100, 100, (nAtoms, 3))
                structure['coordinates'].append(crd.copy())
                for a in range(nAtoms):
                    xyz.write("%s%s%s %.6f%s %.6f%s %.6f%s\n" % (rndb(),
                        symbols[a], rndb(), crd[a,0], rndb(), crd[a,1], rndb(),
                        crd[a,2], rndb()))
            xyz.close()
            self.data.append(structure)

        # Use various units
        nAtoms = 10
        units = { "angs" : 1, "bohr" : 0.529177209, "nm" : 10 }
        self.units_crd = numpy.random.uniform(-10, 10, (nAtoms, 3))
        for u,f in units.items():
            xyz = open("%s/%s.xyz" % (self.tmpDir, u), "w")
            xyz.write("%d\n\n" % nAtoms)
            for a in range(nAtoms):
                xyz.write("He %.7f %.7f %.7f\n" % tuple(self.units_crd[a,:]/f))
            xyz.close()


    def tearDown(self):

        # Remove files
        #for f in os.listdir(self.tmpDir):
        #    if not f.startswith("."): os.remove(self.tmpDir+"/"+f)

        # Remove directory
        #os.rmdir(self.tmpDir)
        pass


    def test_initResult(self):

        for i in range(self.nFiles):
            absolute = "%s/%d.xyz" % (self.tmpDir, i)
            traj = mt.Trajectory(absolute)
            self.assertEqual(traj.fileName, absolute)
            self.assertEqual(traj.nAtoms, self.data[i]['nAtoms'])
            symbols = self.data[i]['symbols']
            for a in range(self.data[i]['nAtoms']):
                self.assertEqual(traj.symbols[a], symbols[a])
                m = atomicMasses[symbols[a]]
                self.assertEqual(traj.masses[a], m[1])
                self.assertEqual(traj.aNumbers[a], m[0])

            frameNo = 0
            frame = traj.read()
            while frame:
                diff = frame['coordinates'] - self.data[i]['coordinates'][frameNo]
                maxDiff = numpy.max(numpy.abs(diff))
                self.assertTrue(maxDiff <= 1e-6)
                comment = frame['comment']
                self.assertEqual(comment, self.data[i]['comments'][frameNo])
                frameNo += 1
                frame = traj.read()
            self.assertEqual(frameNo, len(self.data[i]['coordinates']))

    def test_Units(self):

        for u in ["angs", "bohr", "nm"]:
            absolute = "%s/%s.xyz" % (self.tmpDir, u)
            traj = mt.Trajectory(absolute, units=u)
            frame = traj.read()
            diff = frame['coordinates'] - self.units_crd
            maxDiff = numpy.max(numpy.abs(diff))
            self.assertTrue(maxDiff <= 1e-6)


    def test_initExceptions(self):
        self.assertRaises(TypeError, mt.Trajectory)
        self.assertRaises(FileNotFoundError, mt.Trajectory, 'nonexistent.xyz')
        denied = self.tmpDir+"/denied.xyz"
        open(denied, "w").close()
        os.chmod(denied, 0)
        self.assertRaises(PermissionError, mt.Trajectory, denied)
        os.chmod(denied, stat.S_IRWXU)
        xyz = self.tmpDir+"/angs.xyz"
        self.assertRaises(ValueError, mt.Trajectory, xyz, 'x')
        self.assertRaises(ValueError, mt.Trajectory, xyz, format='x')
        self.assertRaises(ValueError, mt.Trajectory, xyz, 'XYZ', 'x')
        self.assertRaises(ValueError, mt.Trajectory, xyz, mode='x')
        self.assertRaises(ValueError, mt.Trajectory, xyz, 'XYZ', 'r', 'x')
        self.assertRaises(ValueError, mt.Trajectory, xyz, units='x')
        self.assertRaises(ValueError, mt.Trajectory, xyz, mode='r', symbols=[])
        self.assertRaises(FileExistsError, mt.Trajectory, xyz, mode='w')
        self.assertRaises(ValueError, mt.Trajectory, xyz, mode='a')
        empty = self.tmpDir+"/empty"
        open(empty, "w").close()
        self.assertRaises(OSError, mt.Trajectory, empty)
        nonempty = self.tmpDir+"/nonempty"
        f = open(nonempty, "w")
        f.write("x\n")
        f.close()
        self.assertRaises(RuntimeError, mt.Trajectory, nonempty)
        

    def test_WriteXYZ(self):

        s=['C']
        x=numpy.array([[1,2,3]], dtype=numpy.float16)
        xyz = self.tmpDir+"/out.xyz"
        t=mt.Trajectory(xyz, 'XYZ', 'w', 'angs', s)
        t.write(x,'blah')


