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


class BasePlotterTests(unittest.TestCase):
    def setUp(self):
        self.qd = qd.QuadraticDifferential()
        self.plot = qd.BasePlotter(self.qd)

    def test_plot_point_empty_quad(self):
        self.plot.add_plotpoint(1)
        self.plot.calculate_trajectories()

        number_of_phases = len(self.plot.phases)
        number_of_trajectories = len(self.plot.trajectories.keys())
        self.assertEqual(number_of_phases, number_of_trajectories)

    def test_plot_point_1pole_quad(self):
        self.qd.add_smplpole(0)
        self.plot.add_plotpoint(1)
        self.plot.calculate_trajectories()

        number_of_phases = len(self.plot.phases)
        number_of_trajectories = len(self.plot.trajectories.keys())
        self.assertEqual(number_of_phases, number_of_trajectories)

    def test_change_in_quad_reflected_in_plotter(self):
        self.qd.add_smplpole(0)
        self.qd.add_dblpole(1j)

        qd_repr = str(self.qd)
        qd_plotter_repr = str(self.plot.qd)

        self.assertEqual(qd_repr, qd_plotter_repr)




if __name__ == '__main__':
    unittest.main(verbosity=2)
