# -*- coding: utf-8 -*-
#
#
#
from intFormas import intfdz, intfdzCurve
from cmath import *
import numpy as np
from scipy.integrate import odeint
from pathos.multiprocessing import Pool
import sys

class QuadraticDifferential:
    """Esta clase codifica una diferencial cuadratica en el plano"""

    def __init__(self):
        self.zeros = []
        self.dblpoles = []
        self.smplpoles = []
        self.phase = complex(1,0) 
        self.saddles = []
        self.plotpoints = []
        self.trajectories = {}

    def __call__(self, z):
        result = self.phase
        for x in self.zeros:
            result = result*(z-x)/abs(z-x)
        for x in self.dblpoles:
            result = result*((z-x)/abs(z-x))**-2
        for x in self.smplpoles:
            result = result*((z-x)/abs(z-x))**-1
        return result

    def add_plotpoint(self, z):
        self.plotpoints.append(z)

    def add_zero(self, z):
        self.zeros.append(z)
        self.clear_trajectories()

    def add_dblpole(self, z):
        self.dblpoles.append(z)
        self.clear_trajectories()

    def add_smplpole(self, z):
        self.smplpoles.append(z)
        self.clear_trajectories()

    def compute_saddles(self):
        pass

    def clear_trajectories(self):
        self.trajectories={}
        self.saddles=[]

    def compute_trajectories(self):
        pool=Pool()
        result = {}
        for p in self.plotpoints:
            if (p, self.phase) not in self.trajectories:
                result[p] = pool.apply_async(self.compute_trajectory, (p,) )
        pool.close()
        pool.join()
        for p in self.plotpoints:
            if (p, self.phase) not in self.trajectories:
                self.trajectories[(p, self.phase)] = result[p].get()
        print('')
                

    def compute_trajectory(self, p):
        ts=TrajectorySolver(self, p)
        e=ts()
        t=Trajectory(ts)
        return t

    def QD(self, z):
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

    def qd_not(self, z):
	result = self.phase
	for x in self.zeros:
            if x!= z:
		result = result*(z-x)/abs(z-x)
	for x in self.dblpoles:
		result = result*((z-x)/abs(z-x))**-2
	for x in self.smplpoles:
		result = result*((z-x)/abs(z-x))**-1
	return result.conjugate()

    def close_2pole(self, z):
        for x in self.dblpoles:
            if abs(z-x) < 10**-2: return False
        return True

    def close_2smplpole(self, z):
        for x in self.smplpoles:
            if abs(z-x) < 10**-2: return False
        return True


class Monodromy:
    """clase que implementa la monodromia de un punto a otro"""

    def __init__(self, first):
        self.monodromia = 0
        self.ma = 0
        self.patch1 = 0
        self.patch2 = 0
        self.source = first
        if (-first.real)>abs(first.imag) and first.imag>0:
		self.patch1=1
	if (-first.real)>abs(first.imag) and first.imag<0:
		self.patch2=1

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

    def dist(self,z):
        if self.patch1 != self.patch2:
            if self.patch1:
                return (-1)**(self.monodromia+self.ma)*ri*(sqrt(z*(-1j)))
            if self.patch2:
                return (-1)**(self.monodromia+self.ma)*ric*(sqrt(z*(1j)))
        else:
            return (-1)**(self.monodromia+self.ma)*sqrt(z)


class Trajectory:

    def __init__(self, traj):
        self._phase = traj._phase
        self.plotpoint = traj.plotpoint
        self.first = traj.first
        self.last = traj.last
        self.first_mon = traj.first_mon
        self.last_mon = traj.last_mon
        self.coordinates = traj.coordinates


class TrajectorySolver:
    """clase que representa una trayectoria y los metodos para calcularla"""

    def __init__(self, quad, plotpoint, phase = None):
        self.qd = quad
        if phase: self._phase = phase
        else: self._phase = quad.phase
        self.plotpoint = plotpoint
        self.first = plotpoint
        self.last = plotpoint
        self.first_mon = Monodromy(quad(plotpoint).conjugate())
        self.last_mon = Monodromy(quad(plotpoint).conjugate())
        self.coordinates = np.array([[plotpoint.real, plotpoint.imag]])

    def __call__(self):
        if self.first == self.last:
            self.tp()
            self.tn()
            sys.stdout.write(".")
            sys.stdout.flush()
        return self.coordinates

    def traject(self, sign, mon, z):
        norma = 0.5
        coord = np.array([[z.real, z.imag]])
        rep = 0
        fin = z
        inicio = z
        ultimo = z
        start = z
        F = lambda y,t,mono: [sign*self.f(y,t,mono)[0],sign*self.f(y,t,mono)[1]]
        while (self.qd.close_2pole(fin) and self.qd.close_2smplpole(fin) and norma < lim and rep < maxreps):
            sol = odeint(F, [inicio.real, inicio.imag], t, mxstep = maxint, args=(mon,))
            fin = complex(sol[-1,0], sol[-1,1])
            mon(self.qd(fin).conjugate())
            if densidadPuntos < abs(fin-ultimo):
                coord = np.append(coord, np.array([[fin.real, fin.imag]]), axis = 0)
                ultimo = fin
            inicio = fin
            norma = abs(inicio)
            if 50 < rep and 0.01 > abs(fin-start):
                break
            rep += 1
        return coord

    def tp(self):
        self.coordinates = np.vstack((self.coordinates,  self.traject(1, self.last_mon, self.last)))
        self.last=complex(self.coordinates[-1][0], self.coordinates[-1][1])
        return

    def tn(self):
        self.coordinates = np.vstack((self.traject(-1, self.first_mon, self.first)[::-1], self.coordinates))
        self.first=complex(self.coordinates[0][0], self.coordinates[0][1])
        return

    def f(self, y, t, mono):
        x = y[0]
        y = y[1]
        z = mono.dist(self.qd(complex(x,y)).conjugate())
        z *= normav
        if abs(complex(x, y)) > 1:
            z *= abs(complex(x,y))
        return [z.real, z.imag]
   


lim = 30
maxreps = 50000
maxint = 100
t= np.linspace(0,.5,2) #intervalo temporal
densidadPuntos=0.05
normav = 0.001


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


