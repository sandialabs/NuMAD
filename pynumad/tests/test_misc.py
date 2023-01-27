import unittest
import numpy as np

from pynumad.utils.misc_utils import LARCetaL, LARCetaT, _parse_data

class TestMisc(unittest.TestCase):

    def test_larcetat(self):
        alp0 = 53.0
        etat = LARCetaT(alp0)
        self.assertAlmostEqual(etat, 0.2867453857588078)

    def test_larcetal(self):
        SL = 10000000000.0
        YC = -10000000000.0
        alp0 = 53.0
        etal = LARCetaL(SL, YC, alp0)
        self.assertAlmostEqual(etal,-0.7610479585895458)

    def test_parse(self):
        cases = [50.0, '1e10', '6.4023E8']
        self.assertEqual(_parse_data(cases[0]),50.0)
        self.assertEqual(_parse_data(cases[1]),10000000000.0)
        truths = np.isclose(_parse_data(cases), np.array([50.0, 10000000000.0, 640230000.0]))
        self.assertTrue(truths.all())