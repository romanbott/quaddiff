"""Quadratic Differential unit tester."""
import unittest  # pylint: disable=wrong-import-position

import os
import sys
import time
import cmath as cm
import numpy as np
sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import quaddiff as qd  # pylint: disable=wrong-import-position
import quaddiff.core.constants as constants


class QuadraticDifferentialTests(unittest.TestCase):
    def setUp(self):
        self.qd = qd.QuadraticDifferential()  # pylint: disable=invalid-name

    def test_non_empty_quaddiff(self):
        self.qd.add_zero(1)
        self.assertEqual(-1, self.qd(0))

    def test_quaddiff_zeros(self):
        self.qd.zeros += [1, 2, 3, 4, 5, 6, 7]
        for zero in self.qd.zeros:
            self.assertEqual(0, self.qd(zero))

    def test_poles(self):
        self.qd.add_smplpole(1)
        self.assertEqual(self.qd(1), qd.INF)

        self.qd.add_dblpole(2)
        self.assertEqual(self.qd(2), qd.INF)

    def test_proximity(self):
        self.qd.add_smplpole(1)
        self.qd.add_dblpole(-1)
        self.assertFalse(self.qd.close_2pole(2))
        self.assertTrue(self.qd.close_2pole(1.001))
        self.assertTrue(self.qd.close_2pole(-1.001))
        self.qd.sensitivity = 1
        self.assertTrue(self.qd.close_2pole(1))

    def test_phase_change(self):
        self.qd.add_zero(1)
        first = self.qd(1 + 1j)

        self.qd.phase = 1j
        second = self.qd(1 + 1j)

        self.assertNotEqual(first, second)
        self.assertEqual(1j*first, second)

    def test_save_and_load(self):
        self.qd.add_zero(1)
        self.qd.add_zero(0)
        self.qd.add_zero(-1j)

        self.qd.add_smplpole(1+1j)
        self.qd.add_smplpole(-2 + 1j)
        self.qd.add_smplpole(3j)
        self.qd.add_smplpole(2 - 1j)

        self.qd.add_dblpole(2)
        self.qd.add_dblpole(3 + 0.5j)

        self.qd.phase = cm.exp(.45 * cm.pi * 2j)

        qd_repr = str(self.qd)

        name = 'test_quad'
        self.qd.save('/tmp', name=name)

        reconstructed = qd.QuadraticDifferential.from_file('/tmp', name=name)
        reconstructed_repr = str(reconstructed)

        self.assertEqual(qd_repr, reconstructed_repr)
        os.remove(os.path.join('/tmp', name + '.json'))




class MonodromyTests(unittest.TestCase):
    def setUp(self):
        self.mono = qd.Monodromy(1+0*1j)

    def test_trivial_dist(self):
        self.assertEqual(self.mono(1), 1)

    def test_monodromy_change(self):
        point = 1 + 1j
        dist1 = self.mono(point)

        trajectory = [cm.exp(t*cm.pi*2j) for t in np.linspace(0, 1, 10)]
        for point in trajectory:
            self.mono.update(point)

        dist2 = self.mono(point)

        self.assertNotEqual(dist1, dist2)


class TrajectoryTests(unittest.TestCase):
    def setUp(self):
        self.qd = qd.QuadraticDifferential()

        self.point = 1 + 0j
        self.trajectory = qd.TrajectorySolver(self.qd)

    def test_trivial_trajectory_calculation_1(self):
        trajectory = self.trajectory.calculate(self.point)

        for point in trajectory:
            self.assertEqual(point[1], 0)

    def test_trivial_trajectory_calculation_2(self):
        trajectory = self.trajectory.calculate(self.point, phase=-1)

        for point in trajectory:
            self.assertEqual(point[0], 1)

    def test_trivial_trajectory_calculation_3(self):
        point = 0 + 0*1j
        trajectory = self.trajectory.calculate(point, phase=1j)

        for point in trajectory:
            self.assertEqual(round(point[0], 4), -round(point[1], 4))

    def test_boundary_condition(self):
        trajectory = self.trajectory.calculate(self.point)

        last = complex(*trajectory[-1])
        first = complex(*trajectory[0])

        first_close_2_boundary = constants.LIM - 0.5 <= abs(first) and \
            abs(first) <= constants.LIM + 0.5
        last_close_2_boundary = constants.LIM - 0.5 <= abs(last) and \
            abs(last) <= constants.LIM + 0.5

        self.assertTrue(first_close_2_boundary)
        self.assertTrue(last_close_2_boundary)

    def test_zdz_trajectory(self):
        self.qd.phase = 1
        self.qd.add_zero(0)

        trajectory = self.trajectory.calculate(1)

        last = complex(*trajectory[-1])
        first = complex(*trajectory[0])
        first_is_in_x_axis = round(first.imag, 4) == 0.0
        first_close_2_boundary = constants.LIM - 0.5 <= abs(first) and \
            abs(first) <= constants.LIM + 1.0
        self.assertTrue(first_close_2_boundary)
        self.assertTrue(first_is_in_x_axis)

    def test_z2dz_trajectory(self):
        self.qd.phase = 1
        self.qd.add_zero(0)
        self.qd.add_zero(0)

        trajectory = self.trajectory.calculate(1+0.001*1j)

        last = complex(*trajectory[-1])
        first = complex(*trajectory[0])
        first_is_in_y_axis = round(first.real, 2) == 0.0
        first_close_2_boundary = constants.LIM - 0.5 <= abs(first) and \
            abs(first) <= constants.LIM + 1.0
        last_is_in_x_axis = round(last.imag, 2) == 0.0
        last_close_2_boundary = constants.LIM - 0.5 <= abs(last) and \
            abs(last) <= constants.LIM + 0.5
        self.assertTrue(first_close_2_boundary)
        self.assertTrue(first_is_in_y_axis)
        self.assertTrue(last_close_2_boundary)
        self.assertTrue(last_is_in_x_axis)

    def test_qdiff_with_dblpole_trajectory(self):
        self.qd.phase = cm.rect(1, cm.pi/3) 
        self.qd.zeros = []
        self.qd.add_dblpole(0)

        points = [cm.rect(1, coords[0]) for coords in np.random.random((10, 2))]
        for point in points:
            trajectory = self.trajectory.calculate(point)
            last = complex(*trajectory[-1])
            first = complex(*trajectory[0])
            first_close_2_zero = abs(first) <= 0.1
            last_close_2_zero = abs(first) <= 0.1
            converges_2_zero = first_close_2_zero | last_close_2_zero
            self.assertTrue(converges_2_zero)


class BasePlotterTests(unittest.TestCase):
    def setUp(self):
        self.qd = qd.QuadraticDifferential()
        self.plot = qd.BasePlotter(self.qd)

    @unittest.skip("")
    def test_plot_point_empty_quad(self):
        self.plot.add_plotpoint(1)
        self.plot.compute_trajectories()

        number_of_phases = len(self.plot.phases)
        number_of_trajectories = len(self.plot.trajectories.keys())
        self.assertEqual(number_of_phases, number_of_trajectories)

    @unittest.skip("")
    def test_plot_point_1pole_quad(self):
        self.qd.add_smplpole(0)
        self.plot.add_plotpoint(1)
        self.plot.compute_trajectories()

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
        self.plot.compute_trajectories()

        number_of_phases = len(self.plot.phases)
        number_of_trajectories = len(self.plot.trajectories.keys())
        self.assertEqual(number_of_phases * 4, number_of_trajectories)

    def test_save_and_load_trajectories(self):
        self.qd.add_smplpole(0)
        self.qd.add_smplpole(1 + 1j)
        self.qd.add_smplpole(2 - 1j)
        self.qd.add_smplpole(3j)

        self.plot.add_plotpoint(1)
        self.plot.add_plotpoint(1j)
        self.plot.add_plotpoint(-1j)
        self.plot.compute_trajectories()

        prev_trajectories = self.plot.trajectories

        name = 'test_quad'
        self.plot.save_trajectories('/tmp', name=name)
        loaded_trajectories = self.plot.from_file('/tmp', name=name)

        self.assertEqual(prev_trajectories, loaded_trajectories)

if __name__ == '__main__':
    unittest.main(verbosity=2)
