# -*- coding: utf-8 -*-
#
#
#


class QuadraticDifferential:
    """Esta clase codifica una diferencial cuadratica en al esfera"""

    def __init__(self):
        self.zeros = []
        self.dblpoles = []
        self.smplpoles = []
        self.phase = complex(1,0) 

    def __call__(self, z):
        result = self.phase
        for x in self.zeros:
            result = result*(z-x)/abs(z-x)
        for x in self.dblpoles:
            result = result*((z-x)/abs(z-x))**-2
        for x in self.smplpoles:
            result = result*((z-x)/abs(z-x))**-1
        return result

    def add_zero(self, z):
        self.zeros.append(z)

    def add_dblpole(self, z):
        self.dblpoles.append(z)

    def add_smplpole(self, z):
        self.smplpoles.append(z)


class Monodromy:

    def __init__(self, first):
        self.monodromia = 0
        self.ma = 0
        self.patch1 = 0
        self.patch2 = 0
        self.source = first

    def __call__(self, target):
        if self.source.real <= 0 and \
           self.source.imag > 0 and \
           target.real <= 0 and \
           target.imag < 0:
            self.monodromia = (self.monodromia+1) % 2
            self.ma = (self.ma+1) % 2
        if self.source.real <= 0 and \
           self.source.imag < 0 and \
           target.real <= 0 and \
           target.imag > 0:
            self.monodromia = (self.monodromia+1) % 2
            self.ma = (self.ma+1) % 2
        if (-1)*self.source.real > abs(self.source.imag):
            if target.real < 0 and \
              target.imag > (-1)*target.real:
                self.patch1 = (self.patch1+1) % 2
                self.ma = 0
            if target.real < 0 and \
              target.imag < target.real:
                self.patch2 = (self.patch2+1) % 2
                self.ma = 0
        if self.source.imag < 0 and\
          self.source.imag < self.source.real:
            if (-1)*target.real > abs(target.imag):
                self.patch2 = (self.patch2+1) % 2
        if self.source.imag > 0 and\
          self.source.imag > (-1)*self.source.real:
            if (-1)*target.real > abs(target.imag):
                self.patch1 = (self.patch1+1) % 2
        if self.patch1 and self.patch2:
            self.patch1 = 0
            self.patch2 = 0
        self.source = target


