# -*- coding: utf-8 -*-
#
#
#
from intFormas import intfdz, intfdzCurve
from cmath import *

class QuadraticDifferential:
    """Esta clase codifica una diferencial cuadratica en el plano"""

    def __init__(self):
        self.zeros = []
        self.dblpoles = []
        self.smplpoles = []
        self.phase = complex(1,0) 
        self.saddles = []

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

    def compute_saddles(self):
        pass

    def DQ(self, z):
        result = 1.0
	for x in self.zeros:
		result = result*(z-x)
	for x in self.dblpoles:
		result = result*(z-x)**-2
	for x in self.smplpoles:
		result = result*(z-x)**-1
        if abs(result) > 10**4:
            print("mierda")
	return result

    def dqNot(self, z):
	result = self.phase
	for x in self.zeros:
            if x!= z:
		result = result*(z-x)/abs(z-x)
	for x in self.dblpoles:
		result = result*((z-x)/abs(z-x))**-2
	for x in self.smplpoles:
		result = result*((z-x)/abs(z-x))**-1
	return result.conjugate()


class Monodromy:
    """clase que implementa la monodromia de un punto a otro"""

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















def dist(z, patch1, patch2, monodromia, ma):
	if patch1 != patch2:
		if patch1:
			r=ri*(sqrt(z*(-1j)))
		if patch2:
			r=ric*(sqrt(z*(1j)))
	else:
		r=sqrt(z)
	return (-1)**(monodromia+ma)*r


def faseSilla(z,w,quad):
    f = lambda x: sqrt(quad.DQ(x))
    faseSilla = intfdz(f,z,w)
    return faseSilla/abs(faseSilla)

def faseSillaD(z, w, divisiones, quad):
    z, w = (w-z)*10**-3 +z, (z-w)*10**-3 +w
    incremento = (w-z)/divisiones
    quad.phase = 1
    inicio = quad(z)
    mon = Monodromy(inicio)
    fin = z+incremento
    if (-inicio.real) > abs(inicio.imag) and inicio.imag > 0:
    	mon.patch1=1
    if (-inicio.real) > abs(inicio.imag) and inicio.imag < 0:
    	mon.patch2=1
    inicio = z
    faseSilla=0
    for i in range(divisiones):
        f = lambda x: dist(quad.DQ(x), mon.patch1, mon.patch2, mon.monodromia, mon.ma)
        faseSilla += intfdz(f,inicio,fin)
	mon(quad(fin))
        inicio=fin
        fin+=incremento 
    return faseSilla

ri=sqrt(1j)
ric=sqrt(-1j)
