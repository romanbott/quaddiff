import unittest

import os
import sys
import time
import cmath as cm
import numpy as np
sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import quaddiff as qd  # pylint: disable=wrong-import-position
import quaddiff.core.constants as constants
from quaddiff.plot.matplotlibplotter import MatplotlibPlotter



class MatplotlibPlotterTests(unittest.TestCase):
    def setUp(self):
        self.qd = qd.QuadraticDifferential()
        self.plot = qd.MatplotlibPlotter(self.qd)



    def test_complex2XY(self):
        X, Y = self.plot._complex2XY([1+1j, 2j])
        print X, Y
        self.assertEqual(X, (1.0,0.0))
        self.assertEqual(Y, (1.0,2.0))
        
        


if __name__ == '__main__':
    unittest.main(verbosity=2)
