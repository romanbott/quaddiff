"""Quadratic differential package for Python."""
import cmath as cm
import numpy as np

from .core.quaddiff import QuadraticDifferential
from .core.monodromy import Monodromy
from .core.trajectory import TrajectorySolver
from .core.trajectory import Trajectory
from .plot.baseplotter import BasePlotter
from .analyzer import Analyzer
from .plot import *
from .utils import INF


class QD(object):
    def __init__(self, qd=None, plotter=None):
        self.qd = QuadraticDifferential()
        self.plotter = MatplotlibPlotter(self.qd)
        self.analyzer = Analyzer(self.qd)

    def add_zero(self, zero):
        self.qd.add_zero(zero)
        self.plotter.clear_trajectories()

    def add_simple_pole(self, pole):
        self.qd.add_smplpole(pole)
        self.clear_trajectories()

    def add_double_pole(self, pole):
        self.qd.add_dblpole(pole)
        self.clear_trajectories()

    def clear(self):
        self.qd.clear()
        self.plotter.clear()

    def clear_trajectories(self):
        self.plotter.clear_trajectories()

    def add_plotpoint(self, point):
        self.plotter.add_plotpoint(point)

    def add_phase(self, phase):
        self.plotter.add_phase(phase)

    def phase_plot(self, phase):
        self.plotter.get_phase_plot(phase)

    def point_plot(self, point):
        self.plotter.get_point_plot(point)

    def animated_plot(self, **kwargs):
        self.plotter.animate(**kwargs)

    def get_phases(self):
        return self.plotter.phases

    def get_plotpoints(self):
        return self.plotter.plotpoints

    def add_phase_mesh(self, N=6):
        mesh = [cm.rect(1, 2 * cm.pi * t) for t in np.linspace(0, 1, N)]
        self.plotter.phases += mesh

    def add_points_mesh(self, N=6):
        xlim = self.plotter.xlim
        ylim = self.plotter.ylim

        mesh = [
            complex(x, y)
            for x in np.linspace(xlim[0], xlim[1], N)
            for y in np.linspace(ylim[0], ylim[1], N)]

        self.plotter.plotpoints += mesh
