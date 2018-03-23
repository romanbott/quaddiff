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

    # def test_empty_quaddiff(self):
    #     self.assertRaises(ValueError, self.qd, 0)

    def test_non_empty_quaddiff(self):
        self.qd.zeros.append(1)
        self.assertEqual(-1, self.qd(0))

    def test_quaddiff_zeros(self):
        self.qd.zeros += [1, 2, 3, 4, 5, 6, 7]
        for zero in self.qd.zeros:
            self.assertEqual(0, self.qd(zero))

    def test_poles(self):
        self.qd.smplpoles.append(1)
        self.assertEqual(self.qd(1), qd.INF)

        self.qd.dblpoles.append(2)
        self.assertEqual(self.qd(2), qd.INF)

    def test_proximity(self):
        self.qd.smplpoles.append(1)
        self.qd.dblpoles.append(-1)
        self.assertFalse(self.qd.close_2pole(2))
        self.assertTrue(self.qd.close_2pole(1.001))
        self.assertTrue(self.qd.close_2pole(-1.001))
        self.qd.sensitivity = 1
        self.assertTrue(self.qd.close_2pole(1))

    def test_phase_change(self):
        self.qd.zeros.append(1)
        first = self.qd(1 + 1j)

        self.qd.phase = 1j
        second = self.qd(1 + 1j)

        self.assertNotEqual(first, second)
        self.assertEqual(1j*first, second)


class MonodromyTests(unittest.TestCase):
    def setUp(self):
        self.mono = qd.Monodromy(1+0*1j)

    def test_trivial_dist(self):
        self.assertEqual(self.mono.dist(1), 1)

    def test_monodromy_change(self):
        point = 1 + 1j
        dist1 = self.mono.dist(point)

        trajectory = [cm.exp(t*cm.pi*2j) for t in np.linspace(0, 1, 10)]
        for point in trajectory:
            self.mono.update(point)

        dist2 = self.mono.dist(point)

        self.assertNotEqual(dist1, dist2)

    # def test_retrace_trajectory(self):
    #     point = 1 + 1j
    #     dist1 = self.mono.dist(point)

    #     trajectory = [cm.exp(t*cm.pi*2j) for t in np.linspace(0, 2, 1000)]
    #     for point in trajectory:
    #         self.mono.update(point)

    #     dist2 = self.mono.dist(point)

    #     self.assertEqual(dist1, dist2)


class TrajectoryTests(unittest.TestCase):
    def setUp(self):
        self.qd = qd.QuadraticDifferential()

        initial_point = 1 + 0j
        self.trajectory = qd.Trajectory(self.qd, initial_point)

    def test_trivial_trajectory_calculation_1(self):
        trajectory = self.trajectory.calculate()

        for point in trajectory:
            self.assertEqual(point[1], 0)

    def test_trivial_trajectory_calculation_2(self):
        trajectory = self.trajectory.calculate(phase=-1.01)

        for point in trajectory:
            self.assertEqual(point[0], 1)

    def test_trivial_trajectory_calculation_3(self):
        self.trajectory.point = 0 + 0*1j
        trajectory = self.trajectory.calculate(phase=1j)

        for point in trajectory:
            self.assertEqual(round(point[0], 4), -round(point[1], 4))

    def test_boundary_condition(self):
        trajectory = self.trajectory.calculate()

        last = complex(*trajectory[-1])
        first = complex(*trajectory[0])

        first_close_2_boundary = constants.LIM - 0.5 <= abs(first) and \
            abs(first) <= constants.LIM + 0.5
        last_close_2_boundary = constants.LIM - 0.5 <= abs(last) and \
            abs(last) <= constants.LIM + 0.5

        self.assertTrue(first_close_2_boundary)
        self.assertTrue(last_close_2_boundary)

    def test_trivial_trajectory_phase_changed_from_quad(self):
        self.qd.phase = 1j
        self.trajectory.point = 0 + 0*1j

        trajectory = self.trajectory.calculate()

        for point in trajectory:
            self.assertEqual(round(point[0], 4), -round(point[1], 4))

    def test_zdz_trajectory(self):
        self.qd.phase = 1
        self.qd.add_zero(0)
        self.trajectory.point = 1

        trajectory = self.trajectory.calculate()

        last = complex(*trajectory[-1])
        first = complex(*trajectory[0])
        first_is_in_x_axis = round(first.imag,4) == 0.0
        first_close_2_boundary = constants.LIM - 0.5 <= abs(first) and \
            abs(first) <= constants.LIM + 1.0
        self.assertTrue(first_close_2_boundary)
        self.assertTrue(first_is_in_x_axis)

    def test_z2dz_trajectory(self):
        self.qd.phase = 1
        self.qd.add_zero(0)
        self.qd.add_zero(0)
        self.trajectory.point = 1+0.001*1j

        trajectory = self.trajectory.calculate()

        last = complex(*trajectory[-1])
        first = complex(*trajectory[0])
        first_is_in_y_axis = round(first.real,2) == 0.0
        first_close_2_boundary = constants.LIM - 0.5 <= abs(first) and \
            abs(first) <= constants.LIM + 1.0
        last_is_in_x_axis = round(last.imag,2) == 0.0
        last_close_2_boundary = constants.LIM - 0.5 <= abs(last) and \
            abs(last) <= constants.LIM + 0.5
        self.assertTrue(first_close_2_boundary)
        self.assertTrue(first_is_in_y_axis)
        self.assertTrue(last_close_2_boundary)
        self.assertTrue(last_is_in_x_axis)

    def test_qdiff_with_dblpole_trajectory(self):
        self.qd.phase = cm.rect(1,cm.pi/3) 
        self.qd.zeros = []
        self.qd.add_dblpole(0)

        points = [cm.rect(1, coords[0]) for coords in np.random.random((10,2))]
        for point in points:
            self.trajectory.point = point
            trajectory = self.trajectory.calculate()
            last = complex(*trajectory[-1])
            first = complex(*trajectory[0])
            first_close_2_zero = abs(first) <= 0.1
            last_close_2_zero = abs(first) <= 0.1
            converges_2_zero = first_close_2_zero | last_close_2_zero
            print(first, last)
            self.assertTrue(converges_2_zero)



if __name__ == '__main__':
    unittest.main(verbosity=2)
