import cmath as cm

from .. import QuadraticDifferential
from .. import TrajectorySolver

from ..core.constants import *  # pylint: disable=wildcard-import


class Analyzer(object):
    """
    Class for finding critical trajectories and other sutff
    Has methods critical_trajectories_zero, etc
    """
    epsilon = 1e-3
    factor = 3
    close_2pole = 1e-2
    max_step = .01
    distance_2limit = DISTANCE_2LIMIT
    stoping_distance = 1e-5

    def __init__(self, quad):
        self.qd = quad  # pylint: disable=invalid-name
        self.critical_trajectories = {}

    def critical_trajectories_zero(self, zero, phase):
        """Computes critical trajectories for the given zero"""

        solver = TrajectorySolver(self.qd)
        solver.close_2pole = self.close_2pole
        solver.close_2zero = self.epsilon / self.factor
        solver.max_step = self.max_step

        new_trajectories = solver.parallel_calculate(
            self.args_critical_trajectories(zero, phase))
        critical_trajectories = {}
        critical_trajectories.update(new_trajectories)

        #mini = 30
        #for pole in self.qd.dblpoles:
        #    for traj in critical_trajectories.values():
        #        mini = min(mini, abs(traj[0] - pole), abs(traj[-1] - pole))

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
            [(point, phase) for point in init_points],
            progressbar=False)
        trajectories = [
            [init_traj[init_p], init_traj[end_p]]
            for init_p, end_p in zip(init_points, end_points)]

        critical = {}
        unsolved = [0, 1, 2]
        while True:
            # Check if solved
            for i in unsolved:
                t1, t2 = trajectories[i]
                if t1.converges(zero, distance_2limit=self.distance_2limit):
                    unsolved.remove(i)
                    critical[(t1.basepoint, phase)] = t1
                elif t2.converges(zero, distance_2limit=self.distance_2limit):
                    unsolved.remove(i)
                    critical[(t2.basepoint, phase)] = t2
                if abs(phase_intervals[i][0] - phase_intervals[i][1]) <= self.stoping_distance:
                    unsolved.remove(i)
                    critical[(t1.basepoint, phase)] = t1

            if len(unsolved) == 0:
                break

            new_args = []
            for i in unsolved:
                middle_phase = (phase_intervals[i]/2 +  phase_intervals[i][1])/2
                new_args.append((zero + cm.rect(self.epsilon, middle_phase), phase))

            new_trajectories = solver.parallel_calculate(new_args, progressbar=False)

            for i, arg in zip(unsolved, new_args):
                middle_trajectory = new_trajectories[arg]
                t1, t2 = trajectories[i]
                if same_zone(t1, middle_trajectory):
                    trajectories[i][0] = middle_trajectory
                    phase_intervals[i][0] = (phase_intervals[i]/2 +  phase_intervals[i][1])/2
                else:
                    trajectories[i][1] = middle_trajectory
                    phase_intervals[i][1] = (phase_intervals[i]/2 +  phase_intervals[i][1])/2

        return unsolved

    def critical_trajectories_zeros_2(self, phase):
        trajectories = {}
        for zero in self.qd.zeros:
            new_trajectories = self.critical_trajectories_zero_2(zero, phase)
            trajectories.update(new_trajectories)

        return trajectories

    def critical_trajectories_zeros(self, phase):
        """Computes critical trajectories for all zeros"""

        solver = TrajectorySolver(self.qd)
        solver.close_2pole = self.close_2pole
        solver.close_2zero = self.epsilon / self.factor
        solver.max_step = self.max_step
        trayectories = {}
        args = []
        for zero in self.qd.zeros:
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
        arguments = [(point, phase)
                     for point in points]
        return arguments

def same_zone(trajectory1, trajectory2, close=0.01):
    first1 = trajectory1[0]
    first2 = trajectory2[0]
    last1 = trajectory1[-1]
    last2 = trajectory2[-1]

    f1f2 = abs(first1 - first2) <= close
    f1l2 = abs(first1 - last2) <= close
    l1f2 = abs(last1 - first2) <= close
    l1l2 = abs(last1 - last2) <= close

    if f1f2:
        return l1l2
    elif f1l2:
        return l1f2
    return False

