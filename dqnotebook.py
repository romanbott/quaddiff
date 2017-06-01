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
import dill as pickle
from matplotlib.collections import LineCollection



def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


def load_object(filename):
    with open(filename, 'rb') as input:
        obj=pickle.load(input)
    return obj

class QuadraticDrawer:
    def __init__(self, quad, figure, ancho=0.2):
        self.figure = figure
        self.point_type = 'c'
        self.quadratic_differential = None
        self.cid_press = None
        self.cid_click = None
        self.set_quad(quad)
        self.ax = figure.add_subplot(111, xlim=(-3,3), ylim=(-3,3), autoscale_on=False)
        self.trajectories = LineCollection([],offsets=offs,color='k',linewidths=ancho,antialiaseds=True)
        self.ax.add_collection(self.trajectories)

    def __call__(self):
        self.figure.canvas.draw()

    def draw_critical(self):
        pass

    def draw_saddles(self):
        pass

    def draw_trajectories(self):
        t=[]
        self.quadratic_differential.compute_trajectories()
        all_trajectories = self.quadratic_differential.trajectories
        for (x,phase) in all_trajectories:
            if phase == self.quadratic_differential.phase:
                t.append(all_trajectories[(x,phase)].coordinates)
        self.trajectories.set_segments(t)


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
            self.draw_trajectories()
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
        self()

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
offs= (0.0, 0.0)
anchoLineas = 0.2
