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

class TrajectoryTests(unittest.TestCase):
    def setUp(self):
        self.qd = qd.QuadraticDifferential()

        self.point = 1 + 0j
        self.trajectory = qd.TrajectorySolver(self.qd)

    def test_trivial_trajectory_calculation_1(self):
        trajectory = self.trajectory.calculate(self.point)

        for point in trajectory:
            self.assertEqual(point.imag, 0)

    # @unittest.skip("")
    def test_trivial_trajectory_calculation_2(self):
        trajectory = self.trajectory.calculate(self.point, phase=-1)

        for point in trajectory:
            self.assertEqual(point.real, 1)

    # @unittest.skip("")
    def test_trivial_trajectory_calculation_3(self):
        point = 0 + 0*1j
        trajectory = self.trajectory.calculate(point, phase=1j)

        for point in trajectory:
            self.assertEqual(round(point.real, 4), -round(point.imag, 4))

    def test_boundary_condition(self):
        trajectory = self.trajectory.calculate(self.point)

        last = trajectory[-1]
        first = trajectory[0]

        first_close_2_boundary = constants.LIM - 1.5 <= abs(first) and \
            abs(first) <= constants.LIM + 1.5
        last_close_2_boundary = constants.LIM - 1.5 <= abs(last) and \
            abs(last) <= constants.LIM + 1.5

        self.assertTrue(first_close_2_boundary)
        self.assertTrue(last_close_2_boundary)

    # @unittest.skip("")
    def test_zdz_trajectory(self):
        self.qd.phase = 1
        self.qd.add_zero(0)

        trajectory = self.trajectory.calculate(1)

        last = trajectory[-1]
        first = trajectory[0]
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

        last = trajectory[-1]
        first = trajectory[0]
        first_is_in_y_axis = round(first.real, 2) == 0.0
        first_close_2_boundary = constants.LIM - 0.5 <= abs(first) and \
            abs(first) <= constants.LIM + 1.0
        last_is_in_x_axis = round(last.imag, 2) == 0.0
        last_close_2_boundary = constants.LIM - 1.5 <= abs(last) and \
            abs(last) <= constants.LIM + 1.5
        self.assertTrue(first_close_2_boundary)
        self.assertTrue(first_is_in_y_axis)
        self.assertTrue(last_close_2_boundary)
        self.assertTrue(last_is_in_x_axis)

    #@unittest.skip("")
    def test_qdiff_with_dblpole_trajectory(self):
        self.qd.phase = cm.rect(1, cm.pi) 
        self.qd.zeros = []
        self.qd.add_dblpole(0)

        points = [cm.rect(1, phase) for phase in np.random.random((5, 1))]
        for point in points:
            trajectory = self.trajectory.calculate(point, phase = cm.rect(1, 0.99*cm.pi) )
            last = trajectory[-1]
            first = trajectory[0]
            first_close_2_zero = abs(first) <= 0.15
            last_close_2_zero = abs(last) <= 0.15
            first_close_2_infty = abs(first) >= 29.15
            last_close_2_infty = abs(last) >= 29.15
            converges_2_zero = first_close_2_zero | last_close_2_zero
            converges_2_infty = first_close_2_infty | last_close_2_infty
            self.assertTrue(converges_2_zero)
            self.assertTrue(converges_2_infty)


    def test_qdiff_with_dblpole_closed(self):
        self.qd.phase = cm.rect(1, cm.pi) 
        self.qd.zeros = []
        self.qd.add_dblpole(0)

        points = [cm.rect(norm, phase) for phase, norm in np.random.random((50, 2))]
        for point in points:
            trajectory = self.trajectory.calculate(point, phase = cm.rect(1, cm.pi) )
            last = trajectory[-1]
            first = trajectory[0]
            first_close_2_start = abs(first-point) <= 0.15
            last_close_2_start = abs(last-point) <= 0.15
            is_closed_loop = first_close_2_start | last_close_2_start
            self.assertTrue(is_closed_loop)


    def test_problematic_trajectory(self):
        self.qd.phase = cm.rect(1, cm.pi/3) 
        self.qd.zeros = []
        self.qd.add_dblpole(0)

        points = [0.6634486760787831+0.748221794797044j]
        for point in points:
            trajectory = self.trajectory.calculate(point)
            last = trajectory[-1]
            first = trajectory[0]
            first_close_2_zero = abs(first) <= 0.1
            last_close_2_zero = abs(last) <= 0.1
            converges_2_zero = first_close_2_zero | last_close_2_zero
            self.assertTrue(converges_2_zero)

    def test_qdiff_with_dblpole_trajectory(self):
        self.qd.phase = cm.rect(1, cm.pi) 
        self.qd.zeros = []
        self.qd.add_dblpole(0)

        points = [cm.rect(1, phase) for phase in np.random.random((5, 1))]
        for point in points:
            trajectory = self.trajectory.calculate(point, phase = cm.rect(1, 0.99*cm.pi) )
            last = trajectory[-1]
            first = trajectory[0]
            first_close_2_zero = abs(first) <= 0.15
            last_close_2_zero = abs(last) <= 0.15
            first_close_2_infty = abs(first) >= 29.15
            last_close_2_infty = abs(last) >= 29.15
            converges_2_zero = first_close_2_zero | last_close_2_zero
            converges_2_infty = first_close_2_infty | last_close_2_infty
            self.assertTrue(converges_2_zero)
            self.assertTrue(converges_2_infty)

    def test_qdiff_with_dblpole_closed(self):
        self.qd.phase = cm.rect(1, cm.pi) 
        self.qd.zeros = []
        self.qd.add_dblpole(0)

        points = [cm.rect(norm, phase) for phase, norm in np.random.random((50, 2))]
        for point in points:
            trajectory = self.trajectory.calculate(point, phase = cm.rect(1, cm.pi) )
            last = trajectory[-1]
            first = trajectory[0]
            first_close_2_start = abs(first-point) <= 0.15
            last_close_2_start = abs(last-point) <= 0.15
            is_closed_loop = first_close_2_start | last_close_2_start
            self.assertTrue(is_closed_loop)


if __name__ == '__main__':
    unittest.main(verbosity=2)
