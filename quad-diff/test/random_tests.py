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


class RandomTests(unittest.TestCase):
    def test_trajectory_simplify(self):
        trajectory_list = [complex(2*x, 3*x) for x in np.linspace(0, 4, 1000)]
        trajectory = qd.Trajectory(trajectory_list)
        trajectory.min_distance = 10

        simplified = trajectory.simplify()
        self.assertEqual(len(simplified), 2)

    def test_trajectory_getitem(self):
        trajectory_list = [complex(2*x, 3*x) for x in np.linspace(0, 4, 1000)]
        trajectory = qd.Trajectory(trajectory_list)

        for i in range(1000):
            original_point = trajectory_list[i]
            new_point = trajectory_list[i]
            self.assertEqual(original_point, new_point)

    def test_trajectory_iteration(self):
        N = 100
        trajectory_list = [complex(2*x, 3*x) for x in np.linspace(0, 4, N)]
        trajectory = qd.Trajectory(trajectory_list)

        for x, y in zip(trajectory_list, trajectory):
            self.assertEqual(x, y)

    def test_intersection(self):
        N = 100
        trajectory_list = [complex(2*x, 3*x) for x in np.linspace(1, -1, N)]
        trajectory = qd.Trajectory(trajectory_list)

        M = 50
        trajectory2_list = [complex(-2*x, 3*x) for x in np.linspace(-1, 1, M)]
        trajectory2 = qd.Trajectory(trajectory2_list)

        trajectory3_list = [complex(-5 + x, 0) for x in np.linspace(0, 3, M)]
        trajectory3 = qd.Trajectory(trajectory3_list)

        trajectory4_list = [complex(-5 + 2*x, 3*x) for x in np.linspace(0, 3, M)]
        trajectory4 = qd.Trajectory(trajectory4_list)

        self.assertTrue(trajectory.intersects(trajectory2))
        self.assertFalse(trajectory.intersects(trajectory3))
        self.assertFalse(trajectory.intersects(trajectory4))

    def test_integration(self):
        N = 100
        quad = qd.QuadraticDifferential()
        trajectory = qd.Trajectory(
            [complex(2*x, 3*x) for x in np.linspace(0, 1, N)])

        integral = quad.integrate(trajectory)
        self.assertEqual(integral, 2 + 3j)

    def test_intrinsic_length(self):
        N = 100
        quad = qd.QuadraticDifferential()
        trajectory = qd.Trajectory(
            [complex(2*x, 3*x) for x in np.linspace(0, 1, N)])

        lenght = quad.intrinsic_length(trajectory)
        self.assertEqual(round(abs(lenght - abs(2 + 3j)), 3), 0)

        length2 = quad.intrinsic_length_2(trajectory)
        self.assertEqual(round(abs(length2 - abs(2 + 3j)), 3), 0)


if __name__ == '__main__':
    unittest.main()
