"""Trajectory Module."""

import sys
import numpy as np
from cmath import *
from scipy.integrate import solve_ivp, odeint

from .obj import Monodromy
from .constants import *

#hay que pasar las constantes constants.py
#MAX_ITERS = 1000000
#MIN_ITERS = 5000
#VELOCITY_SCALE = 0.002
#MAX_INT = 5000
#POINT_DENSITY = 0.01
#LIM = 30
#time  = np.linspace(0, 0.5, 1000)
#solver_rtol = .00001

class TrajectorySolver(object):
    """clase que representa una trayectoria y los metodos para calcularla."""

    parameters = {
        'max_iters': MAX_ITERS,
        'min_iters': MIN_ITERS,
        'velocity_scale': VELOCITY_SCALE,
        'max_int': MAX_INT,
        'point_density': POINT_DENSITY,
        'lim': LIM,
        'time_axis': TIME_AXIS,
        'solver_rtol': SOLVER_RTOL}

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

        trajectory = negative_trajectory[::-1] \
            + [point] + positive_trajectory

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
    lim = parameters.get('lim', LIM)
    time_axis = np.linspace(*parameters.get('time_axis', TIME_AXIS) )
    solver_rtol = parameters.get('solver_rtol', SOLVER_RTOL)
    max_iters = parameters.get('max_iters', MAX_ITERS)
    min_iters = parameters.get('min_iters', MIN_ITERS)
    max_int = parameters.get('max_int', MAX_INT)
    point_density = parameters.get('point_density', POINT_DENSITY)
    velocity_scale = parameters.get('velocity_scale', VELOCITY_SCALE)
    speed_up_lim = lim/10

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
        value = monodromy.dist(quad(comp, phase=phase).conjugate())
        value *= velocity_scale
        if abs(comp) > speed_up_lim:
            value *= abs(comp)
        value *= sign
        return value.real, value.imag

    # Solve for next point while far from poles, far from infinty
    # and less than maxreps interations.
    while ((not quad.close_2pole(initial)) and
           (abs(initial) < lim) and
           (iteration < max_iters)):

        # Calculate solution with odeint
        solution = odeint(
            vector_field,
            [initial.real, initial.imag],
            time_axis,
            mxstep=max_int,
            rtol = solver_rtol,
            args=(monodromy,))

        # Select final point of solution trayectory
        current_point = complex(solution[-1, 0], solution[-1, 1])
        # Update monodromy
        monodromy.update(quad(current_point, phase=phase).conjugate())

        # Add point to trayectory list if far enough from last
        # point added to trayectory list
        if abs(current_point - last) > point_density:
            trajectory.append(current_point)
            last = current_point

        initial = current_point

        # End iteration if current point is close to starting point and
        # min_iterations have passed
        if (iteration > min_iters and
                abs(current_point - starting_point) < 0.01):
            break
        iteration += 1

    return trajectory
