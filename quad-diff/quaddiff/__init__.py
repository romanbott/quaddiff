"""Quadratic differential package for Python."""
from .core.quaddiff import QuadraticDifferential
from .core.monodromy import Monodromy
from .core.old_trajectory import Trajectory
from .core.trajectory import TrajectorySolver
from .plot.baseplotter import BasePlotter
from .plot import *
from .utils import INF

class QD(object):
    def __init__(self, qd=None, plotter=None):
        self.qd = QuadraticDifferential()
        self.plotter = MatplotlibPlotter()
