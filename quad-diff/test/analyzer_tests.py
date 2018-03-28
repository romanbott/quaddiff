import unittest

import os
import sys
import time
import cmath as cm
import logging
import numpy as np
sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import quaddiff as qd  # pylint: disable=wrong-import-position,import-error
import quaddiff.core.constants as constants  # pylint: disable=wrong-import-position,import-error


class AnalyzerTests(unittest.TestCase):
    def setUp(self):
        self.qd = qd.QuadraticDifferential()
        self.analyzer = qd.Analyzer(self.qd)

    def test_critical_trajectories_zero(self):
        num_zeros = np.random.randint(1, 10)
        num_poles = np.random.randint(1, 10)

        random_zeros = np.random.uniform(-5, 5, size=[num_zeros, 2])
        random_poles = np.random.uniform(-5, 5, size=[num_poles, 2])

        for zero in random_zeros:
            self.qd.add_zero(complex(*zero))

        for pole in random_poles:
            self.qd.add_dblpole(complex(*pole))

        critical_trajectories = self.analyzer.critical_trajectories_zeros(1)

        plotter = qd.MatplotlibPlotter(self.qd)
        plotter.plot(critical_trajectories)

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    unittest.main()
