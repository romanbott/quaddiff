import cmath as cm

from .. import QuadraticDifferential
from .. import TrajectorySolver

from ..core.constants import *  # pylint: disable=wildcard-import


class Analyzer(object):
    epsilon = 1e-3
    distance_2limit = DISTANCE_2LIMIT
    stoping_distance = 1e-5

    def __init__(self, quad):
        self.qd = quad  # pylint: disable=invalid-name
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

    def critical_trajectories_zero_2(self, zero, phase):
        C = 2 * cm.pi / 3.0
        phase_intervals = [[t*C, (t+1)*C] for t in range(3)]

        solver = TrajectorySolver(self.qd)
        solver.close_2pole = 1e-2
        solver.max_step = 0.1
        solver.center = zero
        solver.lim = 1

        init_points = [zero + cm.rect(self.epsilon, p1) for p1, p2 in phase_intervals]
        end_points = [zero + cm.rect(self.epsilon, p2) for p1, p2 in phase_intervals]
        init_traj= solver.parallel_calculate(
            [(point, phase) for point in initial_points],
            progressbar=False)
        trajectories = [
            [init_traj[init_p], init_traj[end_p]]
            for init_p, end_p in zip(init_points, end_points)]

        critical = {}
        while True:
            for t1, t2 in trajectories:
                pass

            break

    def critical_trajectories_zeros(self, phase):
        trayectories = {}
        for zero in self.qd.zeros:
            new_trajectories = self.critical_trajectories_zero(zero, phase)
            trayectories.update(new_trajectories)

        return trayectories
