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

atomicMasses = {
    'H':1.008, 'C':12.011, 'N':14.007, 'O':15.999,
    'P':30.973762, 'S':32.06, 'Cl':35.45, 'Na':22.98976928,
    'Cu':63.546, 'Fe':55.845 }
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


    def tearDown(self):

        # Remove files
        #for f in os.listdir(self.tmpDir):
        #    if not f.startswith("."): os.remove(self.tmpDir+"/"+f)

        # Remove directory
        #os.rmdir(self.tmpDir)
        pass


    def test_topoComponents(self):
        for i in range(self.nFiles):
            absolute = "%s/%d.xyz" % (self.tmpDir, i)
            traj = mt.Trajectory(absolute)
            self.assertEqual(traj.fileName, absolute)
            self.assertEqual(traj.nOfAtoms, self.data[i]['nAtoms'])
            symbols = self.data[i]['symbols']
            for a in range(self.data[i]['nAtoms']):
                self.assertEqual(traj.symbols[a], symbols[a])
                m = atomicMasses[symbols[a]]
                self.assertEqual(traj.atomicMasses[a], m)

if __name__ == '__main__':
    unittest.main()
