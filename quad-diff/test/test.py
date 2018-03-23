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
            self.assertEqual(point.imag, 0)

    def test_trivial_trajectory_calculation_2(self):
        """Trajectories of constant qudartical differntials"""
        trajectory = self.trajectory.calculate(phase=1j)

        for point in trajectory:
            self.assertEqual(abs(point), 1)



if __name__ == '__main__':
    unittest.main()
