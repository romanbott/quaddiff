"""Base plotter module."""

from cmath import *
import numpy as np
import logging
import json
import os
from multiprocessing import Pool

from ..utils import MethodProxy
from ..core.trajectory import TrajectorySolver

class BasePlotter(object):
    """ Base Plotter class"""
    phase_points = 10

    def __init__(self, quad, solver_params=None, trajectories=None):
        self.qd = quad

        self.plotpoints = []
        self.trajectories = {} if trajectories is None else trajectories
        self.saddles = []
        self.phases = [
            exp(t*pi*2j) for t in np.linspace(0, 1, self.phase_points)]

        self.solver = TrajectorySolver(self.qd)
        if solver_params is not None:
            self.solver.parameters = solver_params

    def add_plotpoint(self, z):
        self.plotpoints.append(z)

    def make_mesh(self):
        pass

    def compute_saddles(self):
        pass

    def clear_trajectories(self):
        self.trajectories = {}
        self.saddles = []

    def compute_trajectories(self):
        arguments = [(point, phase)
                     for point in self.plotpoints
                     for phase in self.phases
                     if (point, phase) not in self.trajectories]
        # points, phases = zip(*arguments)
        logging.info('Computing {} new trayectories.'.format(len(arguments)))

        pickable_method = MethodProxy(self.solver, self.solver._calculate)

        pool = Pool()
        results = pool.map(pickable_method, arguments)
        pool.close()
        pool.join()

        for arg, res in zip(arguments, results):
            self.trajectories[arg] = res

        logging.info('done.')

    def save_trajectories(self, path, name='quad_diff'):
        def pkey(key):
            point = (key[0].real, key[0].imag)
            phase = (key[1].real, key[1].imag)
            return str((point, phase))
        # Save trajectories
        jsonable_trajectories = {
            pkey(key): value.tolist()
            for key, value in self.trajectories.iteritems()
        }

        fname = os.path.join(path, name + '_points.json')
        with open(fname, 'w') as jsonfile:
            json.dump(jsonable_trajectories, jsonfile)

    def from_file(self, path, name='quad_diff'):
        from ast import literal_eval
        # Trayectory files
        fname = os.path.join(path, name + '_points.json')
        with open(fname, 'r') as jsonfile:
            trajectories = json.load(jsonfile)

        def pkey(key):
            key = literal_eval(key)
            point = complex(*tuple(key[0]))
            phase = complex(*tuple(key[1]))
            return point, phase

        parsed_trajectories = {
            pkey(key): np.array(value)
            for key, value in trajectories.iteritems()
        }

        self.trajectories.update(parsed_trajectories)
        return parsed_trajectories
