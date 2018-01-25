import unittest
import tempfile
import random
import os
import numpy
import moltools as mt

# Return random string of spaces and tabs
def rndb():
    blanks = [' ', '\t']
    l = [ random.choice(blanks) for i in range(random.randint(0, 5)) ]
    random.shuffle(l)
    return "".join(l)

atomicSymbols = [ 'H', 'C', 'N', 'O', 'P', 'S', 'Cl', 'Na', 'Cu', 'Fe' ]
characters = "".join([chr(x) for x in range(ord('A'), ord('z')+1)])

class TestTrajectoryTopologyXYZ(unittest.TestCase):

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
            symbols = [random.choice(atomicSymbols) for x in range(nAtoms)]
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


    def tearDown(self):

        # Remove files
        for f in os.listdir(self.tmpDir):
            if not f.startswith("."): os.remove(self.tmpDir+"/"+f)

        # Remove directory
        os.rmdir(self.tmpDir)


    def test_topoComponents(self):
        for i in range(self.nFiles):
            absolute = "%s/%d.xyz" % (self.tmpDir, i)
            traj = mt.Trajectory(absolute)
            print(traj.nOfAtoms, self.data[i]['nAtoms'])
            self.assertEqual(traj.fileName, absolute)
            self.assertEqual(traj.nOfAtoms, self.data[i]['nAtoms'])


if __name__ == '__main__':
    unittest.main()
