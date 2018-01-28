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

        # Write ten, very messy xyz files
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

        # Write file with extra data
        nAtoms = 10
        crd = numpy.random.uniform(-10, 10, (nAtoms, 4))
        self.extra_data = crd[:,3]
        xyz = open("%s/extra.xyz" % self.tmpDir, "w")
        xyz.write("%d\n\n" % nAtoms)
        for a in range(nAtoms):
            xyz.write("Ar %.3f %.3f %.3f %.6f\n" % tuple(crd[a,:]))
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
                self.assertAlmostEqual(traj.masses[a], m[1])
                self.assertEqual(traj.aNumbers[a], m[0])

    def test_readXYZ(self):

        for i in range(self.nFiles):
            absolute = "%s/%d.xyz" % (self.tmpDir, i)
            traj = mt.Trajectory(absolute)
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
        # Wrong format
        self.assertRaises(ValueError, mt.Trajectory, xyz, format='fmt')
        # Wrong mode
        self.assertRaises(ValueError, mt.Trajectory, xyz, 'mode')
        self.assertRaises(ValueError, mt.Trajectory, xyz, mode='mode')
        # Wrong units
        self.assertRaises(ValueError, mt.Trajectory, xyz, units='units')
        # Symbols in 'read' mode
        self.assertRaises(ValueError, mt.Trajectory, xyz, 'r', [])
        self.assertRaises(ValueError, mt.Trajectory, xyz, 'r', symbols=[])
        # 'append' mode but no symbols
        self.assertRaises(ValueError, mt.Trajectory, xyz, 'a')
        self.assertRaises(ValueError, mt.Trajectory, xyz, mode='a')

        empty = self.tmpDir+"/empty"
        open(empty, "w").close()
        # Empty file and no format given
        self.assertRaises(OSError, mt.Trajectory, empty)

        nonempty = self.tmpDir+"/nonempty"
        f = open(nonempty, "w")
        f.write("x\n")
        f.close()
        # Guessing format fails
        self.assertRaises(RuntimeError, mt.Trajectory, nonempty)
        self.assertRaises(RuntimeError, mt.Trajectory, nonempty, 'a', [])
        # Writing to an existing file
        self.assertRaises(FileExistsError, mt.Trajectory, nonempty, 'w', [], 'XYZ')
        # Wrong format
        self.assertRaises(ValueError, mt.Trajectory, nonempty, 'a', [], 'fmt')
        # Wrong units
        self.assertRaises(ValueError, mt.Trajectory, nonempty, 'a', [], 'XYZ', 'units')
        

    def test_WriteXYZ(self):

        sym = ['C']
        comment = 'blah'
        crd = [1.1111111111,2,3]

        for dtype in [ numpy.float16, numpy.float32,
                       numpy.float64, numpy.float128 ]:

            x = numpy.array([crd], dtype=dtype)
            xyz = self.tmpDir+"/out.xyz"
            traj = mt.Trajectory(xyz, 'w', sym, 'XYZ', 'angs')
            traj.write(x, comment)
            del traj

            traj = mt.Trajectory(xyz)
            self.assertEqual(traj.nAtoms, 1)
            self.assertEqual(traj.symbols, sym)
            frame = traj.read()
            self.assertEqual(frame['comment'], comment)
            diff = crd - frame['coordinates'][0,:]
            maxDiff = numpy.max(numpy.abs(diff))
            # 1e-8 is the precision of write_frame_to_xyz due to the format % 12.8f
            epsilon = max(1e-8, numpy.finfo(dtype).resolution)
            self.assertTrue(maxDiff <= epsilon)

            # Try (and fail) to read more frames
            frame = traj.read()
            self.assertIsNone(frame)

            os.remove(xyz)

    def test_writeExceptions(self):

        sym = ['C']
        comment = 'blah'
        crd = [1,2,3]
        xyz = self.tmpDir+"/out.xyz"
        # No symbols given
        self.assertRaises(ValueError, mt.Trajectory, xyz, 'w', format='XYZ')
        self.assertRaises(ValueError, mt.Trajectory, xyz, mode='w', format='XYZ')
        traj = mt.Trajectory(xyz, 'w', sym, 'XYZ')
        # Wrong shape of the matrix
        x = numpy.array([], dtype=numpy.float)
        self.assertRaises(RuntimeError, traj.write, x, comment)
        x = numpy.array([[[]]], dtype=numpy.float)
        self.assertRaises(RuntimeError, traj.write, x, comment)
        # Wrong dimensions of the matrix
        x = numpy.array([crd, crd], dtype=numpy.float)
        self.assertRaises(RuntimeError, traj.write, x, comment)
        # Try to read in write mode
        self.assertRaises(RuntimeError, traj.read)
        del traj

    def test_writeMultiple(self):

        nFrames = 10
        sym = ['C', 'O']
        comment = 'blah'
        crd = [[1,2,3], [4,5,6]]
        x = numpy.array(crd, dtype=numpy.float)
        xyz = self.tmpDir+"/out.xyz"
        traj = mt.Trajectory(xyz, 'w', sym, 'XYZ', 'angs')
        for i in range(nFrames):
            traj.write(x, str(i))
        del traj

        traj = mt.Trajectory(xyz)
        count = 0
        frame = traj.read()
        while frame:
            self.assertEqual(frame['comment'], str(count))
            count += 1
            frame = traj.read()
        self.assertEqual(traj.lastFrame, nFrames-1)

    def test_readExtraData(self):
        traj = mt.Trajectory(self.tmpDir+"/extra.xyz")
        frame = traj.read()
        self.assertTrue('extra' in frame)
        diff = frame['extra'] - self.extra_data
        maxDiff = numpy.max(numpy.abs(diff))
        epsilon = max(1e-6, numpy.finfo(numpy.float).resolution)
        self.assertTrue(maxDiff <= epsilon)
