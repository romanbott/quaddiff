"""Trajectory Module."""

import sys
import numpy as np
from cmath import *
from scipy.integrate import solve_ivp

from .monodromy import Monodromy
from .constants import *


class TrajectorySolver(object):
    """clase que representa una trayectoria y los metodos para calcularla."""

    parameters = {
        'max_time': MAX_TIME,
        'velocity_scale': VELOCITY_SCALE,
        'num_points': NUM_POINTS,
        'lim': LIM,
        'max_step': MAX_STEP}

    def __init__(self, quad):
        self.qd = quad

    def calculate(self, point, phase=None):
        """Calculate trayectory."""
        if phase is None:
            phase = self.qd.phase

        positive_trajectory = calculate_ray(
            point, self.qd, parameters=self.parameters, phase=phase)
        negative_trajectory = calculate_ray(
            point, self.qd, sign=-1, parameters=self.parameters, phase=phase)

        trajectory = list(reversed(negative_trajectory)) + positive_trajectory[1:]

        return trajectory

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
    num_points = parameters.get('num_points', NUM_POINTS)
    max_step = parameters.get('max_step', MAX_STEP)
    velocity_scale = parameters.get('velocity_scale', VELOCITY_SCALE)
    lim = parameters.get('lim', LIM)

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

    def vector_field(t, y):
        real, img = y
        comp = complex(real, img)
        value = sqrt_monodromy(quad(comp, phase=phase).conjugate())
        value *= velocity_scale
        if abs(comp) > 1:
            value *= abs(comp)
        value *= sign
        return value.real, value.imag

    # Termination events
    def far_away(t, y):
        comp = complex(*y)
        if abs(comp) > lim:
            return 0
        else:
            return 1
    far_away.terminal = True

    def close_2pole(t, y):
        comp = complex(*y)
        if quad.close_2pole(comp):
            return 0
        else:
            return 1
    close_2pole.terminal = True

    def close_2start(t, y):
        comp = complex(*y)
        if (abs(comp - starting_point) <= 1e-2) and t > 100:
            return 0
        else:
            return 1
    close_2start.terminal = True

    # Calculate solution with solve_ivp
    # time_axis = np.linspace(0, max_time, num_points)
    solution = solve_ivp(
        vector_field,
        (0, max_time),
        np.array([initial.real, initial.imag]),
        events=[far_away, close_2pole, close_2start],
        max_step=max_step)

    solution_complex_list = [complex(*point) for point in solution['y'].T]
    return solution_complex_list
