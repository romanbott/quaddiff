"""Utils for Quadratic Differential Package."""


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


# lim = 30
# maxreps = 500000
# maxint = 500
# t= np.linspace(0,.1,10) #intervalo temporal
# densidadPuntos=0.01
# normav = 0.001
# normav = 0.0002
# def fase_silla(z,w,quad):
#     longitud_silla=integrar(z,w,quad)
#     fase_silla=longitud_silla/abs(longitud_silla)
#     fase_silla=fase_silla**2
#     return fase_silla.conjugate()
# ri=sqrt(1j)
# ric=sqrt(-1j)
# def integrar(x,y,quad, pasos=10**4):
#     paso = (y-x)/pasos
#     medio_paso = paso/2
#     mon=Monodromy(quad(x))
#     integral = 0
#     for i in range(pasos):
#         #integral += mon.dist(quad(x+medio_paso))*paso*sqrt(abs(quad.QD(x+medio_paso)))
#         integral += paso*(sqrt(abs(quad.QD(x)))*mon.dist(quad(x))   +  sqrt(abs(quad.QD(x+paso)))*mon.dist(quad(x+paso))     )/2
#         x+=paso
#         mon(quad(x))
#     return integral