"""Module for Monodromy class."""

import cmath as cm
import logging


class Monodromy(object):
    """Monodromy object class."""

    def __init__(self, point):
        # Internal state variables
        self.phase = cm.phase(point)
        self.point = point / abs(point)

    def update(self, target):
        new_point = target / abs(target)
        change = new_point / self.point
        self.point = new_point
        arg_change = cm.phase(change)
        #if abs(arg_change) > (cm.pi / 2):
        #    logging.info('Large steps can lead to erroneous monodromy')
        self.phase += arg_change

    def __call__(self, z):
        self.update(z)
        if (self.phase - cm.pi) % (4 * cm.pi) < (2 * cm.pi):
            factor = -1
        else:
            factor = 1
        return factor * cm.sqrt(z)

