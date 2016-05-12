# -*- coding: utf-8 -*-
from scipy.integrate import quad
from itertools import izip
def intfdz(f,z,w):
    gama = lambda x: (w-z)*x+z
    integrandoReal = lambda x: ((w-z)*f(gama(x))).real
    integrandoImag = lambda x: ((w-z)*f(gama(x))).imag
    #raizI = lambda x: ((w-z)*sqrt(f(gama(x)))).imag
    #raizR = lambda x: ((w-z)*sqrt(f(gama(x)))).real
    integralReal = quad(integrandoReal, 0, 1,limit=200)
    integralImag = quad(integrandoImag, 0, 1,limit=200)
    return complex(integralReal[0], integralImag[0])
def intfdzCurve(f,curva):
    resultado = 0
    for inicio, fin in izip(curva,curva[1:]):
       resultado += intfdz(f,inicio,fin) 
    return resultado
