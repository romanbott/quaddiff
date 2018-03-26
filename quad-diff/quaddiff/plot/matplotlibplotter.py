import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib.collections import LineCollection
from matplotlib.animation import FuncAnimation

from .baseplotter import BasePlotter

class MatplotlibPlotter(BasePlotter):
    name = 'Matplotlib'
    linewidths = 1
    colors = ['#000000']
    linestyles = 'solid'
    cmap = None
    xlim = [-5, 5]
    ylim = [-5, 5]
    zero_marker = 'o'
    smplpole_marker = 'x'
    dblpole_marker = '*'
    axis = 'off'

    def plot(self, lines, show=True):
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
        return fig

    def plot_lines(self, lines, ax):
        collection = LineCollection(
            tuple([[(z.real, z.imag) for z in line] for line in lines.values()]),
            linewidths=self.linewidths,
            colors=self.colors,
            linestyles=self.linestyles,
            cmap=self.cmap)
        ax.add_collection(collection)

    def animate(self):
        fig, ax = plt.subplots()
        frames = self.phases

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

        anim = FuncAnimation(fig, update, frames=frames, interval=200)
        plt.show()

    def plot_zeros(self):
        X, Y = zip(*[(x.real, x.imag) for x in self.qd.zeros])
        plt.plot(X, Y, self.zero_marker, label='zeros')

    def plot_smplpoles(self):
        X, Y = zip(*[(x.real, x.imag) for x in self.qd.smplpoles])
        plt.plot(X, Y, self.smplpole_marker, label='simple poles')

    def plot_dblpoles(self):
        X, Y = zip(*[(x.real, x.imag) for x in self.qd.dblpoles])
        plt.plot(X, Y, self.dblpole_marker, label='double poles')
