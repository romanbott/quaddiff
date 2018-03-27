"""Trajectory Module."""

import sys
import numpy as np
import cmath as cm
from scipy.integrate import solve_ivp

from ..utils import simplify_trajectory
from .monodromy import Monodromy
from .constants import *


class Trajectory(object):
    distance_2line = DISTANCE_2LINE
    min_distance = MIN_DISTANCE
    distance_2limit = DISTANCE_2LIMIT

    def __init__(self, trajectory, basepoint):
        self.trajectory = trajectory
        self.basepoint = basepoint

    def simplify(self):
        simplified = simplify_trajectory(
            self.trajectory,
            distance_2line=self.distance_2line,
            min_distance=self.min_distance)
        return simplified

    def __getitem__(self, key):
        return self.trajectory[key]

    def __setitem__(self, key, value):
        self.trajectory[key] = value

    def __iter__(self):
        return iter(self.trajectory)

    def __reversed__(self):
        return reversed(self.trajectory)

    def __contains__(self, item):
        return item in self.trajectory

    def __len__(self):
        return len(self.trajectory)

    def intersects(self, other, simplify=False):
        if simplify:
            line = self.simplify()
            other = other.simplify()
        else:
            line = self.trajectory
            other = other.trajectory

        line_array = np.array(line)
        other_array = np.array(other)

        # Create meshgrid to compare all segments simultaneously
        # Mesgrids are shifted by one to obtain segments ends.
        Z1, W1 = np.meshgrid(line_array[:-1], other_array[:-1])
        Z2, W2 = np.meshgrid(line_array[1:], other_array[1:])

        # Transform to take one segment to [0, 1] (x -> (x - z1)/ (z2 - z1))
        T1 = ((W2 - Z1) / (Z2 - Z1))
        T2 = ((W1 - Z1) / (Z2 - Z1))

        # First vector crosses the line defined by the second?
        condition1 = T1.imag * T2.imag <= 0

        # Crosses in between?
        zero = np.divide(  # pylint: ignore=no-member
            T2.imag * T1.real - T1.imag * T2.real,
            T2.imag - T1.imag,
            where=condition1)
        condition2 = (zero >= 0) & (zero <= 1)

        intersections = condition1 & condition2

        return intersections.any()

    def converges(self, point):
        start = self[0]
        end = self[-1]

        start_close_to_point = abs(point - start) <= self.distance_2limit
        end_close_to_point = abs(point - end) <= self.distance_2limit
        return start_close_to_point | end_close_to_point


class TrajectorySolver(object):
    """clase que representa una trayectoria y los metodos para calcularla."""

    parameters = {
        'max_time': MAX_TIME,
        'velocity_scale': VELOCITY_SCALE,
        'lim': LIM,
        'max_step': MAX_STEP,
        'close_2pole': CLOSE_2POLE,
        'center': 0j}

    def __init__(self, quad):
        self.qd = quad

    def calculate(self, point, phase=None):
        """Calculate trajectory."""
        if phase is None:
            phase = self.qd.phase

        positive_trajectory = calculate_ray(
            point, self.qd, parameters=self.parameters, phase=phase)
        negative_trajectory = calculate_ray(
            point, self.qd, sign=-1, parameters=self.parameters, phase=phase)

        trajectory = list(reversed(negative_trajectory)) + positive_trajectory[1:]
        return Trajectory(trajectory, point)

    def _calculate(self, arg):
        point, phase = arg
        return self.calculate(point, phase=phase)


def calculate_ray(
        starting_point,
        quad,
        sign=1,
        parameters={},
        phase=None):
    """Calculate a ray solution to the Quadratic Differential."""
    max_time = parameters.get('max_time', MAX_TIME)
    max_step = parameters.get('max_step', MAX_STEP)
    velocity_scale = parameters.get('velocity_scale', VELOCITY_SCALE)
    lim = parameters.get('lim', LIM)
    center = parameters.get('center', 0j)
    pole_sensitivity = parameters.get('close_2pole', CLOSE_2POLE)

    # Monodromy
    sqrt_monodromy = Monodromy(quad(starting_point))

    # Iteration counter
    iteration = 0

    # Reference points
    initial = starting_point
    last = starting_point

    # Trajectory points
    trajectory = []

    # Check for new phase
    if phase is None:
        phase = quad.phase

    def vector_field(t, y):  # pylint: disable=invalid-name
        real, img = y
        comp = complex(real, img)
        value = sqrt_monodromy(quad(comp, phase=phase).conjugate())
        value *= velocity_scale
        if abs(comp) > 1:
            value *= abs(comp)
        value *= sign
        return value.real, value.imag

    # Termination events
    def far_away(t, y):  # pylint: disable=invalid-name
        comp = complex(*y)
        if abs(comp - center) > lim:
            return 0
        else:
            return 1
    far_away.terminal = True

    def close_2pole(t, y):  # pylint: disable=invalid-name
        comp = complex(*y)
        if quad.close_2pole(comp, pole_sensitivity):
            return 0
        else:
            return 1
    close_2pole.terminal = True

    def close_2start(t, y): # pylint: disable=invalid-name
        comp = complex(*y)
        if (abs(comp - starting_point) <= CLOSE_2START) and t > 100:
            return 0
        else:
            return 1
    close_2start.terminal = True

    # Calculate solution with solve_ivp
    solution = solve_ivp(
        vector_field,
        (0, max_time),
        np.array([initial.real, initial.imag]),
        events=[far_away, close_2pole, close_2start],
        max_step=max_step)

    solution_complex_list = [complex(*point) for point in solution['y'].T]
    return solution_complex_list
