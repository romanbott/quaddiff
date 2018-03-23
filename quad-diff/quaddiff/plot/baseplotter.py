"""Base plotter module."""

from multiprocessing import Pool

class BasePlotter(object):
    """ Base Plotter class"""

    def __init__(self, quadratic_differential):
        self._qd = quadratic_differential

        self.plotpoints = []
        self.trajectories = []
        self.saddles = []

    def add_plotpoint(self, z):
        self.plotpoints.append(z)

    def compute_saddles(self):
        pass

    def clear_trajectories(self):
        self.trajectories = {}
        self.saddles = []

    def compute_trajectories(self):
        pool = Pool()
        result = {}
        for p in self.plotpoints:
            if (p, self.phase) not in self.trajectories:
                result[p] = pool.apply_async(self.compute_trajectory, (p,))
        pool.close()
        pool.join()
        for p in self.plotpoints:
            if (p, self.phase) not in self.trajectories:
                self.trajectories[(p, self.phase)] = result[p].get()

    def compute_trajectory(self, p):
        ts = TrajectorySolver(self, p)
        e = ts()
        t = Trajectory(ts)
        return t