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
    """Analyzer test class"""
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


    def test_saddle_trajectories(self):

        zero1 = 0
        zero2 = 1+1j
        zero3 = 5j
        self.qd.dblpoles = [0.5+2j, -3j-2, -4j+5]
        pole = -0.5-0.5j
        self.qd.add_smplpole(pole)
        self.qd.add_dblpole(pole+3j)
        self.qd.zeros = [zero1, zero2, zero3, zero2+zero3, -zero3]
        #self.qd.zeros = [zero1]
        self.analyzer.epsilon = 0.001
        self.analyzer.max_step = 0.01
        self.analyzer.factor = 2
        #integration_curve = qd.Trajectory([zero1, pole])
        #phase = self.qd.integrate(integration_curve.refine(max_distance=.000001))
        #print(phase)
        #phase /= abs(phase)
        #phase = phase**-2
        #print(phase)
        #self.qd.phase = phase
        #phase2 = self.qd.integrate(integration_curve.refine(max_distance=.00003))
        #phase2 /= abs(phase2)
        #phase2 = phase2**-2
        #print((phase2))
        self.qd.phase = -1
        saddles = self.analyzer.saddle_trajectory(zero2+zero3, zero3, [zero2+zero3, zero3])


        critical_trajectories = self.analyzer.critical_trajectories_zeros(-1)

        plotter = qd.MatplotlibPlotter(self.qd, trajectories=critical_trajectories)
        plotter.saddle_trajectories.update(saddles)
        plotter.plot(plotter.get_trajectories(simplify=False))


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    unittest.main()
