"""Trajectory Module."""

import sys
import numpy as np
import cmath as cm
import logging
from scipy.integrate import solve_ivp

from ..utils import simplify_trajectory
from .monodromy import Monodromy
from .constants import *

class Trajectory(object):

    def __init__(self, trajectory, basepoint=None):
        self.trajectory = trajectory
        self.basepoint = basepoint

    def simplify(self, distance_2line=DISTANCE_2LINE, min_distance=MIN_DISTANCE):
        simplified = simplify_trajectory(
            self.trajectory,
            distance_2line=distance_2line,
            min_distance=min_distance)
        return Trajectory(simplified, basepoint=self.basepoint)

    def refine(self, max_distance=MAX_DISTANCE):
        refined_trajectory = []
        for point, next_point in zip(self[:-1], self[1:]):
            dist = abs(point - next_point)
            if dist > max_distance:
                times = np.arange(0, dist, max_distance)
                direction = (next_point - point)
                direction /= abs(direction)
                refined_subinterval = point + times * direction
                refined_trajectory += refined_subinterval.tolist()
            else:
                refined_trajectory += [point]
        print(len(refined_trajectory))
        return refined_trajectory

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
        else:
            line = self.trajectory

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

    def converges(self, point, distance_2limit=DISTANCE_2LIMIT):
        start = self[0]
        end = self[-1]

        start_close_to_point = abs(point - start) <= distance_2limit
        end_close_to_point = abs(point - end) <= distance_2limit
        return start_close_to_point | end_close_to_point


class TrajectorySolver(object):
    """clase que representa una trayectoria y los metodos para calcularla."""

    max_time = MAX_TIME
    velocity_scale = VELOCITY_SCALE
    lim = LIM
    max_step = MAX_STEP
    close_2pole = CLOSE_2POLE
    close_2start = CLOSE_2START
    close_2zero = CLOSE_2ZERO
    center = 0j

    def __init__(self, quad):
        self.qd = quad  # pylint: disable=invalid-name

    def get_parameters(self):
        parameters = {
            'max_time': self.max_time,
            'max_step': self.max_step,
            'velocity_scale': self.velocity_scale,
            'lim': self.lim,
            'close_2pole': self.close_2pole,
            'close_2start': self.close_2start,
            'close_2zero': self.close_2zero,
            'center': self.center}
        return parameters

    def calculate(self, point, phase=None):
        """Calculate trajectory."""
        if phase is None:
            phase = self.qd.phase

        parameters = self.get_parameters()

        positive_trajectory = calculate_ray(
            point, self.qd, parameters=parameters, phase=phase)
        negative_trajectory = calculate_ray(
            point, self.qd, sign=-1, parameters=parameters, phase=phase)

        trajectory = list(reversed(negative_trajectory)) + \
            positive_trajectory[1:]
        return Trajectory(trajectory, point)

    def _calculate(self, arg):
        point, phase = arg
        trajectory = self.calculate(point, phase=phase)
        return (arg, trajectory)


def calculate_ray(
        starting_point,
        quad,
        sign=1,
        parameters=None,
        phase=None):

    """Calculate a ray solution to the Quadratic Differential."""
    # Check for new phase
    if phase is None:
        phase = quad.phase

    # Check for parameters
    if parameters is None:
        parameters = {}

    # Monodromy
    sqrt_monodromy = Monodromy(quad(starting_point))

    velocity_scale = parameters.get('velocity_scale', VELOCITY_SCALE)
    def vector_field(t, y):  # pylint: disable=invalid-name
        comp = complex(*y)
        value = sqrt_monodromy(quad(comp, phase=phase, normalize=True).conjugate())
        value *= velocity_scale
        if abs(comp) > parameters.get('lim', LIM) / 3.0:
            value *= abs(comp)
        value *= sign
        return value.real, value.imag

    # Termination events
    far = parameters.get('lim', LIM)
    center = parameters.get('center', 0j)
    def far_away(t, y):  # pylint: disable=invalid-name
        comp = complex(*y)
        if abs(comp - center) > far:
            return 0
        else:
            return 1
    far_away.terminal = True

    close = parameters.get('close_2pole', CLOSE_2POLE)
    def close_2pole(t, y):  # pylint: disable=invalid-name
        comp = complex(*y)
        return quad.distance_2poles(comp) - close
    close_2pole.terminal = True

    close = parameters.get('close_2start', CLOSE_2START)
    def close_2start(t, y):  # pylint: disable=invalid-name
        comp = complex(*y)
        if t >= 100:
            return abs(comp - starting_point) - close
        else:
            return 1
    close_2start.terminal = True

    close = parameters.get('close_2zero', CLOSE_2ZERO)
    def close_2zero(t, y):
        comp = complex(*y)
        return quad.distance_2zeros(comp) - close
    close_2zero.terminal = True

    # Calculate solution with solve_ivp
    max_time = parameters.get('max_time', MAX_TIME)
    max_step = parameters.get('max_step', MAX_STEP)
    solution = solve_ivp(
        vector_field,
        (0, max_time),
        np.array([starting_point.real, starting_point.imag]),
        events=[far_away, close_2pole, close_2start, close_2zero],
        max_step=max_step)

    return [complex(*point) for point in solution['y'].T]
