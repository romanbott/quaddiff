# -*- coding: utf-8 -*-
from scipy.integrate import quad
def intfdz(f,z,w):
    gama = lambda x: (w-z)*x+z
    integrandoReal = lambda x: ((w-z)*f(gama(x))).real
    integrandoImag = lambda x: ((w-z)*f(gama(x))).imag
    #raizI = lambda x: ((w-z)*sqrt(f(gama(x)))).imag
    #raizR = lambda x: ((w-z)*sqrt(f(gama(x)))).real
    integralReal = quad(integrandoReal, 0, 1)
    integralImag = quad(integrandoImag, 0, 1)
    return complex(integralReal[0], integralImag[0])
