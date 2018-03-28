import unittest
import os
import sys
import time
import cmath as cm
import numpy as np
sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from quad_tests import QuadraticDifferentialTests
from monodromy_tests import MonodromyTests
from base_plotter_tests import BasePlotterTests
from matplotlib_plotter_tests import MatplotlibPlotterTests



if __name__ == '__main__':
    unittest.main(verbosity=2)
