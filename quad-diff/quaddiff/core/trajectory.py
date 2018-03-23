"""Trajectory Module."""

import sys
import numpy as np
from cmath import *
from scipy.integrate import odeint

from .monodromy import Monodromy

NUM_POINTS = 10000
MAX_TIME = 100
VELOCITY_SCALE = 0.02
MAX_INT = 500


class Trajectory(object):
    """clase que representa una trayectoria y los metodos para calcularla."""

    parameters = {
        'max_time': MAX_TIME,
        'velocity_scale': VELOCITY_SCALE,
        'max_int': MAX_INT,
        'num_points': NUM_POINTS}

    def __init__(self, quad, point, phase=None):
        self.qd = quad
        self.point = point

        if phase:
            self._phase = phase
        else:
            self._phase = quad.phase

    def calculate(self, phase=None):
        """Calculate trayectory."""
        if phase is None:
            phase = self._phase

        positive_trajectory = calculate_ray(
            self.point, self.qd, parameters=self.parameters, phase=phase)
        negative_trajectory = calculate_ray(
            self.point, self.qd, sign=-1, parameters=self.parameters, phase=phase)

        trajectory = np.concatenate(
            [negative_trajectory[::-1],
            [[self.point.real, self.point.imag]],
            positive_trajectory])

        return trajectory


def calculate_ray(
        starting_point,
        quad,
        sign=1,
        parameters={},
        phase=None):
    """Calculate a ray solution to the Quadratic Differential."""
    max_time = parameters.get('max_time', MAX_TIME)
    num_points = parameters.get('num_points', NUM_POINTS)
    max_int = parameters.get('max_int', MAX_INT)
    velocity_scale = parameters.get('velocity_scale', VELOCITY_SCALE)

    # Monodromy
    monodromy = Monodromy(starting_point)

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

    def vector_field(y, t, monodromy):
        real, img = y
        comp = complex(real, img)
        monodromy.update(quad(comp, phase=phase).conjugate())
        value = monodromy.dist(quad(comp, phase=phase).conjugate())
        value *= velocity_scale
        if abs(comp) > 1:
            value *= abs(comp)
        value *= sign
        return value.real, value.imag

    # Calculate solution with odeint
    time_axis = np.linspace(0, max_time, num_points)
    solution = odeint(
        vector_field,
        [initial.real, initial.imag],
        time_axis,
        mxstep=max_int,
        args=(monodromy,))

    return solution
