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

    @unittest.skip("")
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

    @unittest.skip("")
    def test_critical_trajectories_zero_large_epsilon(self):

        self.qd.add_zero(-1.0)
        self.qd.add_zero(1.0)
        self.analyzer.epsilon = 0.1
        critical_trajectories = self.analyzer.critical_trajectories_zeros(-1)

        plotter = qd.MatplotlibPlotter(self.qd)
        plotter.plot(critical_trajectories)

    def test_saddle_trajectories_zero_large_epsilon(self):

        self.qd.add_zero(0.0)
        self.qd.add_zero(1.0+1j)
        self.qd.add_zero(1j)
        #self.qd.add_zero(2.5+2.5j)
        self.analyzer.epsilon = 0.001
        self.analyzer.max_step = 0.01
        self.analyzer.factor = 2
        curva = [t*(1j+1) for t in np.linspace(0.00000001, .9999999, 50000)]
        #curva = [1j+t for t in np.linspace(0.00000001, .9999999, 50000)]
        phase = self.qd.integrate(curva)
        phase /= abs(phase)
        phase = phase**-2
        print(phase)
        self.qd.phase = phase
        phase2 = self.qd.integrate(curva)
        phase2/= abs(phase2)
        phase2 = phase2**-2
        print((phase2))


        critical_trajectories = self.analyzer.critical_trajectories_zeros(phase)

        plotter = qd.MatplotlibPlotter(self.qd)
        plotter.plot(critical_trajectories)



if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    unittest.main()
