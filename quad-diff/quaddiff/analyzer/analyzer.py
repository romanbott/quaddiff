import cmath as cm

from .. import QuadraticDifferential
from .. import TrajectorySolver


class Analyzer(object):
    epsilon = 1e-3

    def __init__(self, quad):
        self.qd = quad
        self.critical_trajectories = {}

    def critical_trajectories_zero(self, zero, phase):
        linearized_quaddiff = self.qd(zero, ignore_zero=True)
        phase_cbrt = cm.exp(cm.log(phase)/3)
        lqd_cbrt = cm.exp(cm.log(linearized_quaddiff)/3)
        unit_cbrt = cm.rect(1, 2 * cm.pi / 3.0)

        critical_phase = self.epsilon * unit_cbrt / (phase_cbrt * lqd_cbrt)

        solver = TrajectorySolver(self.qd)
        solver.close_2pole = 1e-2
        solver.max_step = .1

        points = [
            zero + critical_phase * (unit_cbrt)**j
            for j in range(3)]

        critical_trajectories = {
            (point, phase): solver.calculate(point, phase)
            for point in points}

        mini = 30
        for pole in self.qd.dblpoles:
            for traj in critical_trajectories.values():
                mini = min(mini, abs(traj[0] - pole), abs(traj[-1] - pole))

        return critical_trajectories

    def critical_trajectories_zeros(self, phase):
        trayectories = {}
        for zero in self.qd.zeros:
            new_trajectories = self.critical_trajectories_zero(zero, phase)
            trayectories.update(new_trajectories)

        return trayectories
