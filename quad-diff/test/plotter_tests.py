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


class BasePlotterTests(unittest.TestCase):
    def setUp(self):
        self.qd = qd.QuadraticDifferential()
        self.plot = qd.BasePlotter(self.qd)

    @unittest.skip("")
    def test_plot_point_empty_quad(self):
        self.plot.add_plotpoint(1)
        self.plot.calculate_trajectories()

        number_of_phases = len(self.plot.phases)
        number_of_trajectories = len(self.plot.trajectories.keys())
        self.assertEqual(number_of_phases, number_of_trajectories)

    @unittest.skip("")
    def test_plot_point_1pole_quad(self):
        self.qd.add_smplpole(0)
        self.plot.add_plotpoint(1)
        self.plot.calculate_trajectories()

        number_of_phases = len(self.plot.phases)
        number_of_trajectories = len(self.plot.trajectories.keys())
        self.assertEqual(number_of_phases, number_of_trajectories)

    @unittest.skip("")
    def test_change_in_quad_reflected_in_plotter(self):
        self.qd.add_smplpole(0)
        self.qd.add_dblpole(1j)

        qd_repr = str(self.qd)
        qd_plotter_repr = str(self.plot.qd)

        self.assertEqual(qd_repr, qd_plotter_repr)

    @unittest.skip("")
    def test_plot_4points_4smplpoles(self):
        self.qd.add_smplpole(0)
        self.qd.add_smplpole(1 + 1j)
        self.qd.add_smplpole(2 - 1j)
        self.qd.add_smplpole(3j)

        self.plot.add_plotpoint(1)
        self.plot.add_plotpoint(-1)
        self.plot.add_plotpoint(1j)
        self.plot.add_plotpoint(-1j)
        self.plot.calculate_trajectories()

        number_of_phases = len(self.plot.phases)
        number_of_trajectories = len(self.plot.trajectories.keys())
        self.assertEqual(number_of_phases * 4, number_of_trajectories)

    @unittest.skip("")
    def test_save_and_load_trajectories(self):
        self.qd.add_smplpole(0)
        self.qd.add_smplpole(1 + 1j)
        self.qd.add_smplpole(2 - 1j)
        self.qd.add_smplpole(3j)

        self.plot.add_plotpoint(1)
        self.plot.add_plotpoint(1j)
        self.plot.add_plotpoint(-1j)
        self.plot.calculate_trajectories()

        prev_trajectories = self.plot.trajectories

        name = 'test_quad'
        self.plot.save_trajectories('/tmp', name=name)
        loaded_trajectories = self.plot.from_file('/tmp', name=name)

        self.assertEqual(prev_trajectories, loaded_trajectories)


class MatplotlibPlotterTests(unittest.TestCase):
    def setUp(self):
        self.qd = qd.QuadraticDifferential()
        self.plot = qd.MatplotlibPlotter(self.qd)

    def test_plot(self):
        self.qd.add_zero(2+.01j)
        self.qd.add_smplpole(1j-0.1)
        self.qd.add_dblpole(-1j+0.1)

        self.plot.add_plotpoint(0)
        self.plot.add_plotpoint(1)
        self.plot.add_plotpoint(-1)

        self.plot.calculate_trajectories()

        self.plot.get_phase_plot(0)



if __name__ == '__main__':
    unittest.main(verbosity=2)
