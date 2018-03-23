"""Module for Monodromy class."""

from cmath import *


class Monodromy(object):
    """Monodromy object class."""

    def __init__(self, point):
        # Internal state variables
        self.monodromy = False
        self.factor = 1
        self.patch1 = False
        self.patch2 = False

        # current point
        self.current = point

        # Initial update
        if point.real < -abs(point.imag):
            self.patch1 = point.imag > 0
            self.patch2 = point.imag < 0

    def update(self, target):

        if self.current.real <= 0 and target.real <= 0:
            if self.current.imag * target.imag < 0:
                self.monodromy = not self.monodromy

        if (self.current.real < -abs(self.current.imag) and
                target.real < 0):

            self.factor *= (-1)**self.monodromy * self.factor

            if target.imag > -target.real:
                self.patch1 = not self.patch1
            elif target.imag < target.real:
                self.patch2 = not self.patch2

        if (self.current.imag < 0 and
                self.current.imag < self.current.real and
                target.real < -abs(target.imag)):
            self.patch2 = not self.patch2

        if (self.current.imag > 0 and
                self.current.imag > -self.current.real and
                target.real < -abs(target.imag)):
            self.patch1 = not self.patch1

        if self.patch1 and self.patch2:
            self.patch1 = False
            self.patch2 = False

        self.current = target

    def dist(self, z):
        if self.patch1 != self.patch2:
            if self.patch1:
                return self.factor * sqrt(1j) * (sqrt(z * (-1j)))
            if self.patch2:
                return self.factor * sqrt(-1j) * (sqrt(z * (1j)))
        else:
            return self.factor * sqrt(z)
