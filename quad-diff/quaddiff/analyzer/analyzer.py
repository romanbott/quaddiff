import cmath as cm

from .. import QuadraticDifferential
from .. import TrajectorySolver


class Analyzer(object):
    def __init__(self, quad):
        self.qd = quad

    def critical_trajectories_zero(self, zero, phase):
        pass