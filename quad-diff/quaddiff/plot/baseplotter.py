"""Base plotter module."""

import numpy as np
import cmath as cm
import logging
import json
import os

from ..core.constants import *

from ..core.trajectory import TrajectorySolver


class BasePlotter(object):
    """ Base Plotter class"""
    name = 'Base'
    distance_2line = DISTANCE_2LINE
    min_distance = MIN_DISTANCE

    def __init__(self, quad, solver_params=None, trajectories=None):
        self.qd = quad

        self.plotpoints = []
        self.trajectories = {} if trajectories is None else trajectories
        self.phases = []

        self.solver = TrajectorySolver(self.qd)
        if solver_params is not None:
            self.solver.parameters = solver_params

    def clear(self):
        self.plotpoints = []
        self.phases = []
        self.clear_trajectories()

    def add_plotpoint(self, z):
        self.plotpoints.append(z)

    def add_phase(self, phase):
        self.phases.append(phase)

    def add_point_mesh(self, N=6):
        self.plotpoints += [
            complex(x, y)
            for x in np.linspace(-5, 5, N)
            for y in np.linspace(-5, 5, N)]

    def add_phase_mesh(self, N=6):
        self.phases += [
            cm.rect(1, t*2*cm.pi) for t in np.linspace(0, 1, N)]

    def add_vecinity_points(self, point, N=10, epsilon=0.2):
        self.plotpoints += [
            point + cm.rect(epsilon, 2 * cm.pi * t)
            for t in np.linspace(0, 1, N)]

    def clear_trajectories(self):
        self.trajectories = {}

    def calculate_trajectories(self, progressbar=True):
        arguments = [(point, phase)
                     for point in self.plotpoints
                     for phase in self.phases
                     if (point, phase) not in self.trajectories]

        new_trajectories = self.solver.parallel_calculate(
            arguments,
            progressbar=progressbar)
        self.trajectories.update(new_trajectories)

    def get_trajectories(self, phase=None, plotpoint=None, simplify=True):
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

        if simplify:
            distance_2line = self.distance_2line
            min_distance = self.min_distance

            logging.info('Simplifying trayectories')
            selection = {
                key: value.simplify(
                    distance_2line=self.distance_2line,
                    min_distance=self.min_distance)
                for key, value in selection.iteritems()
            }
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

    def get_phase_plot(self, phase, calculate=True, **kwargs):
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
        self.plot(lines, **kwargs)

    def get_point_plot(self, plotpoint, calculate=True, **kwargs):
        if calculate:
            if plotpoint not in self.plotpoints:
                self.add_plotpoint(plotpoint)
            self.calculate_trajectories()
        else:
            if phase not in self.phases:
                msg = "Phase not in phase list. Please"
                msg += " add phase to phase list and re-calculate"
                msg += " the trajectories."
                raise ValueError(msg)
        lines = self.get_trajectories(plotpoint=plotpoint)
        self.plot(lines, **kwargs)

    def __repr__(self):
        msg = '{}Plotter Object:\n'.format(self.name)
        msg += '\tPlotpoints: \n'
        for point in self.plotpoints:
            msg += '\t\t+ {} \n'.format(point)
        msg += '\tPhases: \n'
        for phase in self.phases:
            msg += '\t\t+ {} \n'.format(phase)
        return msg


