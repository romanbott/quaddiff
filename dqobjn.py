# -*- coding: utf-8 -*-
import numpy as np
import sys
from scipy.integrate import odeint
from scipy.integrate import quad
from cmath import *
from multiprocessing import Pool
from itertools import combinations
from intFormas import intfdz, intfdzCurve
from obj import QuadraticDifferential, Monodromy, faseSillaD
from dqn import QuadraticDrawer
import cPickle as pickle


def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


def load_object(filename):
    with open(filename, 'rb') as input:
        obj=pickle.load(input)
    return obj






quad=QuadraticDifferential()

def dq(z):
    return quad(z).conjugate()
def DQ(z):
	evaluacion=1
	for x in ceros:
		evaluacion=evaluacion*(z-x)
	for x in polos:
		evaluacion=evaluacion*(z-x)**-2
	for x in polosSimples:
		evaluacion=evaluacion*(z-x)**-1
	return evaluacion
def dqnorm(z):
	evaluacion=1
	for x in polos:
		evaluacion=evaluacion*(abs(z-x))**-2
	return evaluacion
dqNot=quad.dqNot

###### generacion manual de los puntos
fase = quad.phase
puntos = []
ceros = quad.zeros
polos = quad.dblpoles
polosSimples = quad.smplpoles



#Constantes e inicializacion de variables
ri=sqrt(1j)
ric=sqrt(-1j)
pause = False

########### funciones #########
def dist(z, patch1, patch2, monodromia, ma):
	if patch1 != patch2:
		if patch1:
			r=ri*(sqrt(z*(-1j)))
		if patch2:
			r=ric*(sqrt(z*(1j)))
	else:
		r=sqrt(z)
	return (-1)**(monodromia+ma)*r

def f(y,t,patch1, patch2, monodromia, ma):
	xi=y[0]
	yi=y[1]
	#c=1
	z=dist(dq(complex(xi,yi)), patch1, patch2, monodromia, ma)
	#if abs(z)>10**4:
		#c=10.0**2/abs(z)
		#print c
	if z!=0:
		z=z*normav
	#	z=z/abs(z)*normav
	if abs(complex(xi,yi))>1:
		z=z*abs(complex(xi,yi))
	#z=z*c
	f0 = z.real
	f1 = z.imag
	return [f0, f1]

def fi(y,t, patch1, patch2, monodromia, ma):
	xi=y[0]
	yi=y[1]
	#c=1
	z=-dist(dq(complex(xi,yi)), patch1, patch2, monodromia, ma)
	#if abs(z)>10**4:
		#c=10.0**2/abs(z)
	if z!=0:
		z=z*normav
	#	z=z/abs(z)*normav
	if abs(complex(xi,yi))>1:
		z=z*abs(complex(xi,yi))
	#z=z*c
	f0 = z.real
	f1 = z.imag
	return [f0, f1]
def compute_trajectory(quad, phase, z, norm = 0.5):
    pass

def trayectp(z):
	norma=0.5
	coord=np.array([[z.real,z.imag]])
	rep=0
	inicio=dq(z)
	fin=z
        mon=Monodromy(dq(z))
	if (-inicio.real)>abs(inicio.imag) and inicio.imag>0:
		mon.patch1=1
	if (-inicio.real)>abs(inicio.imag) and inicio.imag<0:
		mon.patch2=1
	inicio = z
	ultimo = z
	start = z
	while ( cercaPolo(fin) and cercaPoloS(fin) and norma < lim and rep<maxreps):
                argumentos = ( mon.patch1, \
                        mon.patch2, \
                        mon.monodromia, \
                        mon.ma)
		sol = odeint(f,[inicio.real, inicio.imag],t,mxstep=maxint,args = argumentos)
		fin= complex(sol[-1,0],sol[-1,1])
		mon(dq(fin))
		if densidadPuntos < abs(fin-ultimo):
			coord=np.append(coord,np.array([[fin.real,fin.imag]]),axis=0)
			ultimo=fin
		inicio = fin
		norma = max(inicio.real,inicio.imag)
		norma = abs(inicio)
		if 50<rep and 0.01>abs(fin-start):
			break
	#	if rep%densidadPuntos==1:
	#		#x=(2*fin.real/(1+fin*fin.conjugate())).real
	#		#y=(2*fin.imag/(1+fin*fin.conjugate())).real
	#		#w=((1-fin*fin.conjugate())/(1+fin*fin.conjugate())).real
	#		if esfera==1:
	#			coord=np.append(coord,np.array([estereografica(fin)]),axis=0)
	#		else:
	#			coord=np.append(coord,np.array([[fin.real,fin.imag]]),axis=0)
			#coord=np.append(coord,np.array([[fin.real,fin.imag,0.0]]),axis=0)
			#np.append(coord,[[x,y,w]],axis=0)
		rep=rep+1
	#print('Trayectoria positiva, '+str(rep)+' repticiones')
	return coord  
def cercaPolo(z):
    global polos
    for x in polos:
        if abs(z-x)<10**-2: return False
    return True
def cercaPoloS(z):
    global polosSimpes
    for x in polosSimples:
        if abs(z-x)<10**-2: return False
    return True
def trayectn(z):
	norma=0.5
	#fin=z
	#x=(2*fin.real/(1+fin*fin.conjugate())).real
	#y=(2*fin.imag/(1+fin*fin.conjugate())).real
	#w=((1-fin*fin.conjugate())/(1+fin*fin.conjugate())).real
	coord=np.array([[z.real,z.imag]])
	#coord=np.array([[fin.real,fin.imag,0.0]])
	rep=0
	inicio=dq(z)
	fin=z
        mon=Monodromy(dq(z))
	if (-inicio.real)>abs(inicio.imag) and inicio.imag>0:
		mon.patch1=1
	if (-inicio.real)>abs(inicio.imag) and inicio.imag<0:
		mon.patch2=1
	inicio = z
	ultimo = z
	start = z
	while ( cercaPolo(fin) and cercaPoloS(fin) and norma < lim and rep<maxreps):
                argumentos = ( mon.patch1, \
                        mon.patch2, \
                        mon.monodromia, \
                        mon.ma)
		sol = odeint(fi,[inicio.real, inicio.imag],t,mxstep=maxint,args=argumentos)
		fin= complex(sol[-1,0],sol[-1,1])
		mon(dq(fin))
		if densidadPuntos< abs(fin-ultimo):
			coord=np.append(coord,np.array([[fin.real,fin.imag]]),axis=0)
			ultimo=fin
		inicio = fin
		norma = max(inicio.real,inicio.imag)
		norma = abs(inicio)
		if 50<rep and 0.01>abs(fin-start):
			break
	#	if rep%densidadPuntos==1:
	#		if esfera==1:
	#			coord=np.append(coord,np.array([estereografica(fin)]),axis=0)
	#		else:
	#			coord=np.append(coord,np.array([[fin.real,fin.imag]]),axis=0)
			#x=(2*fin.real/(1+fin*fin.conjugate())).real
			#y=(2*fin.imag/(1+fin*fin.conjugate())).real
			#w=((1-fin*fin.conjugate())/(1+fin*fin.conjugate())).real
			#coord=np.append(coord,np.array([[fin.real,fin.imag,0.0]]),axis=0)
		rep=rep+1
	#print('Trayectoria negativa, '+str(rep)+' repticiones')
	return coord
def trayectoria(z):
	tn=trayectn(z)
	tp=trayectp(z)
	total=np.vstack((tn[::-1],tp))
	#x=tn[0]
	#x.reverse()
	#x.extend(tp[0])
	#y=tn[1]
	#y.reverse()
	#y.extend(tp[1])
	#z=tn[2]
	#z.reverse()
	#z.extend(tp[2])
	#print('.'),
	sys.stdout.write(".")
	sys.stdout.flush()
	return total


def trayectorias(listaPuntos):
	pool=Pool()
	tray=[]
	tray = pool.map(trayectoria,listaPuntos)
	pool.close()
	pool.join()
	print('!')
	return np.array(tray)












import matplotlib as mpl
#mpl.use("pgf")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


#parametros
maxreps=50000 #repeticiones maximas de la rutina de integracion
maxint=100
fase=rect(1,0.0) #fase de la diferencial
normav=0.001 #determina la norma del campo
normamax = abs(complex(10,10)) #determina el area dentro de la cual se integra el campo
t= np.linspace(0,.5,2) #intervalo temporal
densidadPuntos=0.05 #determina cada cuantas repeticiones de la rutina de integracion se guarda un punto para su graficacion
colorLineas = 'k'
colorLineasCriticas = 'r'
colorLineasSillas= 'b'
colorLineasSillasD= 'g'
anchoLineas = 0.2
actualizaClick = 0#determina si se redibuja la figura cada click
dibujaEjes = 0
raizCubica= rect(1,2*pi/3)
lim = 30



fig = plt.figure(figsize = (9, 9))
ax = fig.add_subplot(111, xlim=(-3,3), ylim=(-3,3), autoscale_on=False)
offs= (0.0, 0.0)
coleccionTrayectorias=[]
col = LineCollection(coleccionTrayectorias,offsets=offs,color=colorLineas,linewidths=anchoLineas,antialiaseds=True)
colCriticas = LineCollection(coleccionTrayectorias,offsets=offs,color=colorLineasCriticas,linewidths=anchoLineas,antialiaseds=True)
colSillas= LineCollection(coleccionTrayectorias,offsets=offs,color=colorLineasSillas,linewidths=anchoLineas,antialiaseds=True)
colSillasD= LineCollection(coleccionTrayectorias,offsets=offs,color=colorLineasSillasD,linewidths=anchoLineas,antialiaseds=True)
ax.add_collection(col, autolim=True)
ax.add_collection(colCriticas, autolim=True)
ax.add_collection(colSillas, autolim=False)
ax.add_collection(colSillasD, autolim=False)


quad_drawer = QuadraticDrawer(quad, fig)


if __name__ == "__main__":
    plt.show()
else:
    plt.ion()
    plt.show()
