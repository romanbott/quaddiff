"""Base plotter module."""

import numpy as np
import cmath as cm
import logging
import json
import os
from tqdm import tqdm
from multiprocessing import Pool

from ..utils import MethodProxy
from ..core.trajectory import TrajectorySolver


class BasePlotter(object):
    """ Base Plotter class"""
    name = 'Base'

    def __init__(self, quad, solver_params=None, trajectories=None):
        self.qd = quad

        self.plotpoints = []
        self.trajectories = {} if trajectories is None else trajectories
        self.saddles = []
        self.phases = [1 + 0j]

        self.solver = TrajectorySolver(self.qd)
        if solver_params is not None:
            self.solver.parameters = solver_params

    def add_plotpoint(self, z):
        self.plotpoints.append(z)

    def add_phase(self, phase):
        self.phases.append(phase)

    def make_mesh(self, N=6):
        self.plotpoints = [
            complex(x, y)
            for x in np.linspace(-5, 5, N)
            for y in np.linspace(-5, 5, N)]

    def make_phase_mesh(self, N=6):
        self.phases = [
            cm.rect(1, t*2*cm.pi) for t in np.linspace(0, 1, N)]

    def compute_saddles(self):
        pass

    def clear_trajectories(self):
        self.trajectories = {}
        self.saddles = []

    def calculate_trajectories(self):
        arguments = [(point, phase)
                     for point in self.plotpoints
                     for phase in self.phases
                     if (point, phase) not in self.trajectories]

        pickable_method = MethodProxy(self.solver, self.solver._calculate)

        pool = Pool()
        results = pool.imap(pickable_method, tqdm(arguments))
        pool.close()
        pool.join()

        for arg, res in zip(arguments, results):
            self.trajectories[arg] = res

    def get_trajectories(self, phase=None, plotpoint=None):
        if plotpoint is not None:
            if plotpoint not in self.plotpoints:
                msg = "Plotpoint selected is not in plotpoint list"
                raise ValueError(msg)

            selection = {
                phase: self.trajectories[(plotpoint, phase)]
                for phase in self.phases}
        elif phase is not None:
            if phase not in self.phases:
                msg = "Phase selected is not in phases list"
                raise ValueError(msg)

            selection = {
                point: self.trajectories[(point, phase)]
                for point in self.plotpoints
            }
        else:
            selection = self.trajectories
        return selection

    def save_trajectories(self, path, name='quad_diff'):
        def pkey(key):
            point = (key[0].real, key[0].imag)
            phase = (key[1].real, key[1].imag)
            return str((point, phase))
        # Save trajectories
        jsonable_trajectories = {
            pkey(key): [[x.real, x.imag] for x in value]
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
            pkey(key): [complex(*x) for x in value]
            for key, value in trajectories.iteritems()
        }

        self.trajectories.update(parsed_trajectories)
        return parsed_trajectories

    def get_phase_plot(self, phase, calculate=True):
        if calculate:
            if phase not in self.phases:
                self.add_phase(phase)
            self.calculate_trajectories()
        else:
            if phase not in self.phases:
                msg = "Phase not in phase list. Please"
                msg += " add phase to phase list and re-calculate"
                msg += " the trajectories."
                raise ValueError(msg)
        lines = self.get_trajectories(phase=phase)
        self.plot(lines)

    def __repr__(self):
        msg = '{}Plotter Object:\n'.format(self.name)
        msg += '\tPlotpoints: \n'
        for point in self.plotpoints:
            msg += '\t\t+ {} \n'.format(point)
        msg += '\tPhases: \n'
        for phase in self.phases:
            msg += '\t\t+ {} \n'.format(phase)
        return msg
