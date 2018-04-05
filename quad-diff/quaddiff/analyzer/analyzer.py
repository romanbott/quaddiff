import cmath as cm

from .. import QuadraticDifferential
from .. import TrajectorySolver, Trajectory

from ..core.constants import *  # pylint: disable=wildcard-import


class Analyzer(object):
    """
    Class for finding critical trajectories and other sutff
    Has methods critical_trajectories_zero, etc
    """
    epsilon = 1e-3
    factor = 3
    close_2pole = 1e-2
    max_step = .05
    lim = LIM
    distance_2limit = DISTANCE_2LIMIT
    stoping_distance = 1e-5

    def __init__(self, quad):
        self.qd = quad  # pylint: disable=invalid-name
        self.critical_trajectories = {}

    def critical_trajectories_zeros(self, phase, zeros=None):
        """Computes critical trajectories for all zeros"""
        if zeros is None:
            zeros = self.qd.zeros

        solver = TrajectorySolver(self.qd)
        solver.lim = self.lim
        solver.close_2pole = self.close_2pole
        solver.close_2zero = self.epsilon / self.factor
        solver.max_step = self.max_step

        trayectories = {}
        args = []
        for zero in zeros:
            args += self.args_critical_trajectories(zero, phase)
        new_trajectories = solver.parallel_calculate(args)
        trayectories.update(new_trajectories)

        return trayectories

    def args_critical_trajectories(self, zero, phase):
        """Returns the position of points to plot critical trajectories"""

        linearized_quaddiff = self.qd(zero, ignore_zero=True, normalize=True)
        phase_cbrt = cm.exp(cm.log(phase)/3)
        lqd_cbrt = cm.exp(cm.log(linearized_quaddiff)/3)
        unit_cbrt = cm.rect(1, 2 * cm.pi / 3.0)

        critical_phase = self.epsilon * unit_cbrt / (phase_cbrt * lqd_cbrt)

        points = [
            zero + critical_phase * (unit_cbrt)**j
            for j in range(3)]
        arguments = [(point, phase) for point in points]
        return arguments

    def saddle_trajectory(self, zero1, zero2, path):
        """Compute saddle trajectory from zero1 to zero2"""
        integration_path = Trajectory(path)
        integral = self.qd.integrate(
            integration_path.refine(max_distance=0.01))
        phase = integral/abs(integral)
        phase = phase**-2
        critical_trajectories = self.critical_trajectories_zeros(
            zeros=[zero1], phase=phase)
        return critical_trajectories
