import unittest

import os
import sys
import time
import cmath as cm
import numpy as np
import logging
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

    # @unittest.skip("")
    def test_phase_plot(self):
        random = np.random.uniform(-5, 5, size=[8, 2])
        self.qd.add_zero(complex(*random[0]))
        self.qd.add_smplpole(complex(*random[1]))
        self.qd.add_dblpole(complex(*random[2]))
        self.qd.add_zero(complex(*random[3]))
        self.qd.add_smplpole(complex(*random[4]))
        self.qd.add_dblpole(complex(*random[5]))
        self.qd.add_smplpole(complex(*random[6]))
        self.qd.add_dblpole(complex(*random[7]))

        self.plot.phases = [1j]
        self.plot.add_point_mesh()
        self.plot.format = 'svg'

        self.plot.get_phase_plot(1j, save='phase_plot', show=False)

    # @unittest.skip("")
    def test_animate(self):
        random = np.random.uniform(-5, 5, size=[8, 2])
        self.qd.add_zero(complex(*random[0]))
        self.qd.add_smplpole(complex(*random[1]))
        self.qd.add_dblpole(complex(*random[2]))
        self.qd.add_zero(complex(*random[3]))
        self.qd.add_smplpole(complex(*random[4]))
        self.qd.add_dblpole(complex(*random[5]))
        self.qd.add_smplpole(complex(*random[6]))
        self.qd.add_dblpole(complex(*random[7]))

        self.plot.add_point_mesh(N=5)
        self.plot.add_phase_mesh(N=8)

        self.plot.animate(save='animate_plot', show=False)

    # @unittest.skip("")
    def test_point_plot(self):
        random = np.random.uniform(-5, 5, size=[8, 2])
        self.qd.add_zero(complex(*random[0]))
        self.qd.add_smplpole(complex(*random[1]))
        self.qd.add_dblpole(complex(*random[2]))
        self.qd.add_zero(complex(*random[3]))
        self.qd.add_smplpole(complex(*random[4]))
        self.qd.add_dblpole(complex(*random[5]))
        self.qd.add_smplpole(complex(*random[6]))
        self.qd.add_dblpole(complex(*random[7]))

        point = complex(*np.random.uniform(-5, 5, size=[2]))
        self.plot.add_plotpoint(point)
        self.plot.add_phase_mesh(N=10)

        self.plot.get_point_plot(point, show=False, save='point_plot')

    def tests_point_vicinity_plot(self):
        random = np.random.uniform(-5, 5, size=[8, 2])
        self.qd.add_zero(complex(*random[0]))
        self.qd.add_smplpole(complex(*random[1]))
        self.qd.add_dblpole(complex(*random[2]))
        self.qd.add_zero(complex(*random[3]))
        self.qd.add_smplpole(complex(*random[4]))
        self.qd.add_dblpole(complex(*random[5]))
        self.qd.add_smplpole(complex(*random[6]))
        pole = complex(*random[7])
        self.qd.add_dblpole(pole)

        self.plot.add_vecinity_points(pole, N=100)
        self.plot.get_phase_plot(1, save='vicinity_plot', show=False)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    unittest.main(verbosity=2)
