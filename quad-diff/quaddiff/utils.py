"""Utils for Quadratic Differential Package."""
from core.constants import *


class Inf(object):
    """Infinity object class."""
    def __add__(self, other):
        if not isinstance(other, (complex, float, int)):
            msg = "Infinity cannot be operated with {}".format(type(other))
            raise TypeError(msg)
        return self

    def __add__(self, other):
        if not isinstance(other, (complex, float, int)):
            msg = "Infinity cannot be operated with {}".format(type(other))
            raise TypeError(msg)
        return self

    def __sub__(self, other):
        if not isinstance(other, (complex, float, int)):
            msg = "Infinity cannot be operated with {}".format(type(other))
            raise TypeError(msg)
        return self

    def __mul__(self, other):
        if not isinstance(other, (complex, float, int)):
            msg = "Infinity cannot be operated with {}".format(type(other))
            raise TypeError(msg)
        if other == 0:
            return 0
        else:
            return self

    def __rmul__(self, other):
        if not isinstance(other, (complex, float, int)):
            msg = "Infinity cannot be operated with {}".format(type(other))
            raise TypeError(msg)
        if other == 0:
            return 0
        else:
            return self

    def __repr__(self):
        return 'Infinity'

INF = Inf()


class MethodProxy(object):
    def __init__(self, obj, method):
        self.obj = obj
        if isinstance(method, basestring):
            self.methodName = method
        else:
            assert callable(method)
            self.methodName = method.func_name

    def __call__(self, *args, **kwargs):
        return getattr(self.obj, self.methodName)(*args, **kwargs)


def simplify_trajectory(
        trajectory,
        distance_2line=DISTANCE_2LINE,
        min_distance=MIN_DISTANCE):
    last = trajectory[0]
    reference_angle = (trajectory[1] - last)
    new_trajectory = [last]

    for i in range(1, len(trajectory) - 1):
        point = trajectory[i]
        component = orthogonal_component(reference_angle, point - last)
        if (component >= distance_2line and
                abs(point - last) >= min_distance):
            last = point
            new_trajectory.append(last)
            reference_angle = trajectory[i + 1] - point

    new_trajectory.append(trajectory[-1])
    return new_trajectory


def orthogonal_component(base, new):
    direction = base / abs(base)
    component = -(new.real * direction.imag) + (new.imag * direction.real)
    return abs(component)
