import os
import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib.collections import LineCollection
from matplotlib.animation import FuncAnimation

from .baseplotter import BasePlotter


class MatplotlibPlotter(BasePlotter):
    name = 'Matplotlib'
    linewidths = .1
    colors = ['#000000']
    linestyles = 'solid'
    cmap = 'jet'
    figsize = None
    xlim = [-5, 5]
    ylim = [-5, 5]
    zero_marker = 'o'
    smplpole_marker = 'x'
    dblpole_marker = '*'
    axis = 'off'
    animation_interval = 200
    format = 'png'

    def plot(self, lines, show=True, save=None, dir='.'):
        fig, ax = plt.subplots()

        ax.set_xlim(self.xlim[0], self.xlim[1])
        ax.set_ylim(self.ylim[0], self.ylim[1])
        self.plot_lines(lines, ax)
        self.plot_zeros()
        self.plot_smplpoles()
        self.plot_dblpoles()
        plt.legend()
        plt.axis(self.axis)

        if show:
            plt.show()
        if save is not None:
            path = os.path.join(dir, save)
            plt.savefig('{}.{}'.format(path, self.format))
        plt.close()

    def plot_lines(self, lines, ax):
        collection = LineCollection(
            tuple([[(z.real, z.imag) for z in line] for line in lines.values()]),
            linewidths=self.linewidths,
            colors=self.colors,
            linestyles=self.linestyles,
            cmap=self.cmap)
        ax.add_collection(collection)

    def animate(self, save=None, show=True):

        self.calculate_trajectories()
        frames = self.phases
        fig, ax = plt.subplots()

        def update(phase):
            lines = self.get_trajectories(phase=phase)
            ax.clear()
            ax.set_xlim(self.xlim[0], self.xlim[1])
            ax.set_ylim(self.ylim[0], self.ylim[1])
            self.plot_zeros()
            self.plot_smplpoles()
            self.plot_dblpoles()
            self.plot_lines(lines, ax)
            plt.legend()
            plt.axis(self.axis)

        interval = self.animation_interval
        anim = FuncAnimation(fig, update, frames=frames, interval=interval)

        if save is not None:
            anim.save(save + '.gif', dpi=80, writer='imagemagick')
        if show:
            plt.show()

    def plot_zeros(self):
        if self.qd.zeros:
            X, Y = self._complex2XY(self.qd.zeros)
            plt.plot(X, Y, self.zero_marker, label='zeros')

    def plot_smplpoles(self):
        if self.qd.smplpoles:
            X, Y = self._complex2XY(self.qd.smplpoles)
            plt.plot(X, Y, self.smplpole_marker, label='simple poles')

    def plot_dblpoles(self):
        if self.qd.dblpoles:
            X, Y = self._complex2XY(self.qd.dblpoles)
            plt.plot(X, Y, self.dblpole_marker, label='double poles')

    def _complex2XY(self, complex_list):
        X, Y = zip(*[(z.real, z.imag) for z in complex_list])
        return X, Y
