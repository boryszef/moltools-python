import unittest
import mdarray as mt


class TestAtomicLists(unittest.TestCase):

    def test_AtomicNumbers(self):

        anum = mt.AtomicNumbers
        self.assertEqual(anum['H'], 1)
        self.assertEqual(anum['C'], 6)
        self.assertEqual(anum['Og'], 118)

    def test_AtomicSymbols(self):

        asym = mt.AtomicSymbols
        self.assertEqual(asym[0], '')
        self.assertEqual(asym[1], 'H')
        self.assertEqual(asym[79], 'Au')
        self.assertEqual(asym[118], 'Og')

    def test_AtomicNames(self):

        anam = mt.AtomicNames
        self.assertEqual(anam[0], '')
        self.assertEqual(anam[1], 'Hydrogen')
        self.assertEqual(anam[48], 'Cadmium')
        self.assertEqual(anam[118], 'Oganesson')

    def test_AtomicMasses(self):

        ams = mt.AtomicMasses
        self.assertAlmostEqual(ams[1], 1.008)
        self.assertAlmostEqual(ams[19], 39.0983)
        self.assertAlmostEqual(ams[80], 200.592)
        self.assertEqual(ams[0], None)
        self.assertEqual(ams[118], None)

    def test_CovalentRadii(self):

        acr = mt.CovalentRadii
        self.assertAlmostEqual(acr['H'], 0.31)
        self.assertAlmostEqual(acr['K'], 2.03)
