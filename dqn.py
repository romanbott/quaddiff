# -*- coding: utf-8 -*-
import numpy as np
import sys
from scipy.integrate import odeint
from scipy.integrate import quad
from cmath import *
from multiprocessing import Pool
from itertools import combinations
from intFormas import intfdz, intfdzCurve
from obj import QuadraticDifferential, Monodromy, faseSillaD, Trajectory
import cPickle as pickle


def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


def load_object(filename):
    with open(filename, 'rb') as input:
        obj=pickle.load(input)
    return obj

class QuadraticDrawer:
    def __init__(self, quad, figure):
        self.figure = figure
        self.point_type = 'c'
        self.quadratic_differential = None
        self.cid_press = None
        self.cid_click = None
        self.set_quad(quad)

    def __call__(self):
        self.figure.canvas.draw()

    def draw_critical(self):
        pass

    def draw_saddles(self):
        pass

    def set_quad(self, quad):
        if quad == self.quadratic_differential: return
        self.quadratic_differential = quad
        if self.cid_press:
            self.figure.canvas.mpl_disconnect(self.cid_press)
        self.cid_press = self.figure.canvas.mpl_connect('key_press_event', self.press)
        if self.cid_click:
            self.figure.canvas.mpl_disconnect(self.cid_click)
        self.cid_click= self.figure.canvas.mpl_connect('button_press_event', self.click)

    def press(self, event):
        if event.key == 'd':
            self.quadratic_differential.phase *= rect(1, -0.1) #hacer metodo en QuadraticDifferential
            print("- Rotando fase: "+str(self.quadratic_differential.phase))
        if event.key == 'a':
            self.quadratic_differential.phase *= rect(1, 0.1) #hacer metodo en QuadraticDifferential
            print("+ Rotando fase: "+str(self.quadratic_differential.phase))
        if event.key == 'r':
            print('redibujando...')
            self()
        if event.key == 'R':
            print('redibujando criticas...')
            self.draw_critical()
        if event.key == 'm':
            print('redibujando sillas...')
            self.draw_saddles()
        if event.key == 'e':
            plt.savefig('figura.pdf')
        if event.key == 'C':
            self.quadratic_differential.clear_zeros()
        if event.key == 'c': self.point_type = 'c'
        if event.key == 'p': self.point_type = 'p'
        if event.key == 'u': self.point_type = 'u'

    def click(self, event):
        if event.button == 1:
            self.quadratic_differential.add_plotpoint(complex(event.xdata, event.ydata))
            print('P')
        if event.button == 3:
            if self.point_type == 'c':
                self.quadratic_differential.add_zero(complex(event.xdata, event.ydata))
                plt.plot([event.xdata],[event.ydata], 'ro')
                print(self.point_type)

            if self.point_type == 'p':
                self.quadratic_differential.add_dblpole(complex(event.xdata, event.ydata))
                plt.plot([event.xdata],[event.ydata], 'go')
                print(self.point_type)
            if self.point_type == 'u':
                self.quadratic_differential.add_smplpole(complex(event.xdata, event.ydata))
                plt.plot([event.xdata],[event.ydata], 'bo')
                print(self.point_type)

import matplotlib as mpl
#mpl.use("pgf")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


#parametros
maxreps=50000 #repeticiones maximas de la rutina de integracion
maxint=100
fase=rect(1,0.0) #fase de la diferencial
normav=0.001 #determina la norma del campo
normamax = abs(complex(10,10)) #determina el area dentro de la cual se integra el campo
t= np.linspace(0,.5,2) #intervalo temporal
densidadPuntos=0.05 #determina cada cuantas repeticiones de la rutina de integracion se guarda un punto para su graficacion
colorLineas = 'k'
colorLineasCriticas = 'r'
colorLineasSillas= 'b'
colorLineasSillasD= 'g'
anchoLineas = 0.2
actualizaClick = 0#determina si se redibuja la figura cada click
dibujaEjes = 0
raizCubica= rect(1,2*pi/3)
lim = 30



fig = plt.figure(figsize = (9, 9))
ax = fig.add_subplot(111, xlim=(-3,3), ylim=(-3,3), autoscale_on=False)
offs= (0.0, 0.0)
coleccionTrayectorias=[]
col = LineCollection(coleccionTrayectorias,offsets=offs,color=colorLineas,linewidths=anchoLineas,antialiaseds=True)
colCriticas = LineCollection(coleccionTrayectorias,offsets=offs,color=colorLineasCriticas,linewidths=anchoLineas,antialiaseds=True)
colSillas= LineCollection(coleccionTrayectorias,offsets=offs,color=colorLineasSillas,linewidths=anchoLineas,antialiaseds=True)
colSillasD= LineCollection(coleccionTrayectorias,offsets=offs,color=colorLineasSillasD,linewidths=anchoLineas,antialiaseds=True)
ax.add_collection(col, autolim=True)
ax.add_collection(colCriticas, autolim=True)
ax.add_collection(colSillas, autolim=False)
ax.add_collection(colSillasD, autolim=False)





quad = QuadraticDifferential()
quad_drawer = QuadraticDrawer(quad, fig)


if __name__ == "__main__":
    plt.show()
else:
    plt.ion()
    plt.show()
