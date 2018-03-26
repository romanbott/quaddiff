import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib.collections import LineCollection

from .baseplotter import BasePlotter

class MatplotlibPlotter(BasePlotter):
    name = 'Matplotlib'

    def plot(self, lines):
        fig = plt.figure()
        for line in lines:
            X, Y = zip(*[(z.real, z.imag) for z in line])
            plt.plot(X, Y)
        plt.show()
