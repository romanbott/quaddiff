"""Module with Quadratic Differential class definition"""
# -*- coding: utf-8 -*-

import os
import json
import numpy as np
from cmath import *
from constants import *

from .monodromy import  Monodromy
from ..utils import INF


class QuadraticDifferential(object):
    """Esta clase codifica una diferencial cuadratica en el plano"""

    def __init__(self, quad=None, phase=None):
        if quad is not None:
            self.zeros = quad.get('zeros')
            self.dblpoles = quad.get('double_poles')
            self.smplpoles = quad.get('simple_poles')
            self.phase = quad.get('phase', complex(1))
        else:
            self.zeros = []
            self.dblpoles = []
            self.smplpoles = []
            self.phase = complex(1)

        if phase is not None:
            self.phase = phase

    def __repr__(self):
        """Quadratic Differential representation string."""
        msg = "Quadratic Differential:\n"
        msg += "\t Zeros: {}\n".format(self.zeros)
        msg += "\t Simple Poles: {}\n".format(self.smplpoles)
        msg += "\t Double Poles: {}\n".format(self.dblpoles)
        return msg

    @property
    def size(self):
        return len(self.zeros) + len(self.dblpoles) + len(self.smplpoles)

    def __call__(self, z, ignore_zero=False, phase=None):
        # Check if any zeros or poles have been defined
        if self.size == 0:
            msg = "Quadratic Differential is empty:\n {}".format(self)

        # Remove argument from list of zeros if ignore_zero
        zeros = self.zeros
        if z in self.zeros:
            if ignore_zero:
                zeros = [zero for zero in self.zeros if zero != z]
            else:
                return 0.0

        # Select Quad Diff phase if no phase has been given
        phase = self.phase if phase is None else phase

        # Multiply all normalized zero's monomials
        zero_contrib = reduce(
            lambda x, y: x * (z - y) / abs(z - y),
            zeros, 1)

        try:
            dblpoles_contrib = reduce(
                lambda x, y: x * ((z - y) / abs(z - y))**(-2),
                self.dblpoles, 1)
            smplpoles_contrib = reduce(
                lambda x, y: x * ((z - y) / abs(z - y))**(-1),
                self.smplpoles, 1)

        # If z is a pole return INF(nity)
        except ZeroDivisionError:
            return INF

        return phase * zero_contrib * dblpoles_contrib * smplpoles_contrib

    def clear(self):
        self.zeros = []
        self.smplpoles = []
        self.dblpoles = []

    def add_zero(self, z):
        self.zeros.append(complex(z))

    def add_dblpole(self, z):
        self.dblpoles.append(complex(z))

    def add_smplpole(self, z):
        self.smplpoles.append(complex(z))

    def _close_2dblpole(self, z, sensitivity):
        for x in self.dblpoles:
            if abs(z-x) < sensitivity:
                return True
        return False

    def _close_2smplpole(self, z, sensitivity):
        for x in self.smplpoles:
            if abs(z-x) < sensitivity:
                return True
        return False

    def close_2pole(self, z, sensitivity):
        return self._close_2smplpole(z, sensitivity=sensitivity) | \
            self._close_2dblpole(z, sensitivity=sensitivity)

    def integrate(self, trajectory):
        starting_point = trajectory[0]
        sqrt_monodromy = Monodromy(self(starting_point))

        value = 0
        for point in trajectory[1:]:
            dz = point - starting_point

            quad_value0 = sqrt_monodromy(self(starting_point))
            quad_value1 = sqrt_monodromy(self(point))
            avg_quad_value = (quad_value0 + quad_value1) / 2

            value += dz * avg_quad_value
            starting_point = point
        return value

    def intrinsic_length(self, trajectory):
        segment_starts = np.array(trajectory[:-1])
        segment_ends = np.array(trajectory[1:])

        starts_value = sqrt(abs(self(segment_starts)))
        ends_value = sqrt(abs(self(segment_ends)))

        average_value = (starts_value + ends_value) / 2

        dl = abs(segment_ends - segment_starts)
        integral_contributions = average_value * dl

        return np.sum(integral_contributions)

    def intrinsic_length_2(self, trajectory):
        starting_point = trajectory[0]
        value = 0
        for point in trajectory[1:]:
            dl = abs(point - starting_point)

            quad_value0 = sqrt(abs(self(starting_point)))
            quad_value1 = sqrt(abs(self(point)))
            avg_quad_value = (quad_value0 + quad_value1) / 2

            value += dl * avg_quad_value
            starting_point = point
        return value

    def save(self, path, name='quad_diff'):
        fname = os.path.join(path, name + '.json')

        data = {
            'zeros': [[x.real, x.imag] for x in self.zeros],
            'simple_poles': [[x.real, x.imag] for x in self.smplpoles],
            'double_poles': [[x.real, x.imag] for x in self.dblpoles],
            'phase': [self.phase.real, self.phase.imag]
        }

        with open(fname, 'w') as jsonfile:
            json.dump(data, jsonfile)

    @classmethod
    def from_file(cls, path, name='quad_diff'):
        fname = os.path.join(path, name + '.json')

        with open(fname, 'r') as jsonfile:
            data = json.load(jsonfile)

        quad_data = {
            key: [complex(*x) for x in value] 
            if key != 'phase' else complex(*value)
            for key, value in data.iteritems()
        }

        return cls(quad=quad_data)
