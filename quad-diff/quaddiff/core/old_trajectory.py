"""Trajectory Module."""

import sys
import numpy as np
from cmath import *
from scipy.integrate import odeint

#from .monodromy import Monodromy
from .obj import Monodromy

MAX_ITERS = 1000000
MIN_ITERS = 5000
VELOCITY_SCALE = 0.002
MAX_INT = 5000
POINT_DENSITY = 0.01
LIM = 30
#maxreps = 500000
#maxint = 500
#t= np.linspace(0,.1,10) #intervalo temporal
#densidadPuntos=0.01
#normav = 0.001
#normav = 0.0002
time  = np.linspace(0, 0.5, 1000)

class Trajectory(object):
    """clase que representa una trayectoria y los metodos para calcularla."""

    parameters = {
        'max_iters': MAX_ITERS,
        'min_iters': MIN_ITERS,
        'velocity_scale': VELOCITY_SCALE,
        'max_int': MAX_INT,
        'point_density': POINT_DENSITY,
        'lim': LIM}

    def __init__(self, quad,  phase=None):
        self.qd = quad

        if phase:
            self._phase = phase
        else:
            self._phase = quad.phase

    def calculate(self, point, phase=None):
        """Calculate trayectory."""

        self.point = point
        if phase is None:
            phase = self._phase

        positive_trajectory = calculate_ray(
            self.point, self.qd, parameters=self.parameters, phase=phase)
        negative_trajectory = calculate_ray(
            self.point, self.qd, sign=-1, parameters=self.parameters, phase=phase)

        trajectory = negative_trajectory[::-1] \
            + [self.point] + positive_trajectory

        return trajectory


def calculate_ray(
        starting_point,
        quad,
        sign=1,
        parameters={},
        phase=None):
    """Calculate a ray solution to the Quadratic Differential."""
    lim = parameters.get('lim', LIM)
    max_iters = parameters.get('max_iters', MAX_ITERS)
    min_iters = parameters.get('min_iters', MIN_ITERS)
    max_int = parameters.get('max_int', MAX_INT)
    point_density = parameters.get('point_density', POINT_DENSITY)
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
        #monodromy.update(quad(comp, phase = phase))
        value = monodromy.dist(quad(comp, phase=phase).conjugate())
        value *= velocity_scale
        if abs(comp) > 3:
            value *= abs(comp)
        value *= sign
        return value.real, value.imag

    # Solve for next point while far from poles, far from infinty
    # and less than maxreps interations.
    while ((not quad.close_2pole(initial)) and
           (abs(initial) < lim) and
           (iteration < max_iters)):

        # Calculate solution with odeint
        time_axis = time
        solution = odeint(
            vector_field,
            [initial.real, initial.imag],
            time_axis,
            mxstep=max_int,
            rtol = .00001,
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
