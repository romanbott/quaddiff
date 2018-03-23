"""Trajectory Module."""

import sys
import numpy as np
from cmath import *
from scipy.integrate import solve_ivp

from .monodromy import Monodromy
from .constants import *


class Trajectory(object):
    """clase que representa una trayectoria y los metodos para calcularla."""

    parameters = {
        'max_time': MAX_TIME,
        'velocity_scale': VELOCITY_SCALE,
        'max_int': MAX_INT,
        'num_points': NUM_POINTS,
        'lim': LIM,
        'max_step': MAX_STEP}

    def __init__(self, quad, point, phase=None):
        self.qd = quad
        self.point = point

    def calculate(self, phase=None):
        """Calculate trayectory."""
        if phase is None:
            phase = self.qd.phase

        positive_trajectory = calculate_ray(
            self.point, self.qd, parameters=self.parameters, phase=phase)
        negative_trajectory = calculate_ray(
            self.point, self.qd, sign=-1, parameters=self.parameters, phase=phase)

        trajectory = np.concatenate(
            [negative_trajectory.T[::-1],
            [[self.point.real, self.point.imag]],
            positive_trajectory.T])

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
    max_step = parameters.get('max_step', MAX_STEP)
    velocity_scale = parameters.get('velocity_scale', VELOCITY_SCALE)
    lim = parameters.get('lim', LIM)

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

    def vector_field(t, y):
        real, img = y
        comp = complex(real, img)
        monodromy.update(quad(comp, phase=phase).conjugate())
        value = monodromy.dist(quad(comp, phase=phase).conjugate())
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

    # Calculate solution with solve_ivp
    # time_axis = np.linspace(0, max_time, num_points)
    solution = solve_ivp(
        vector_field,
        (0, max_time),
        np.array([initial.real, initial.imag]),
        events=[far_away, close_2pole],
        max_step=max_step)


    return solution['y']
