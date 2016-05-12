# -*- coding: utf-8 -*-
import numpy as np
import sys
from scipy.integrate import odeint
from scipy.integrate import romberg 
from scipy.integrate import quad
from cmath import *
from multiprocessing import Pool
from itertools import combinations
from intFormas import intfdz
########################################
#este archivo trabaja en paralelo y guarda los datos en hoyo.npy
#Definicion de la diferencial cuadratica
#def dq(z):
	#return (fase*((1/z**4+z-1+1/(z-complex(0.3,0.7))**3)*((z+1.2)**3))).conjugate() #diferencial cuadratica
	#return (fase.conjugate()*(1-z)*(z+1)*(z-complex(1,1))*(z-complex(1,-1))*(z-complex(-1,1))*(z-complex(-1,-1))*(z-complex(0,-1.0))*(z-complex(0,1.0))).conjugate()
#	return (fase.conjugate()*(z)*(z-complex(1,1))*(z-complex(1,-1))).conjugate()
	#return (fase*(z-complex(0,1))*(z-complex(0,-1))).conjugate()# raices en las tres raices de la unidad, una coneccion silla
	#return (fase*z*(z-rect(1,0))*(z-rect(1,8*pi/5))*(z-rect(1,6*pi/5))*(z-rect(1,2*pi/5))*(z-rect(1,4*pi/5))).conjugate()# raices en las tres raices de la unidad, una coneccion silla
	#return (fase*(z)).conjugate() #diferencial cuadratica
	#return (fase*(complex(-1,0)/z**2+(1-z)**5)).conjugate() #diferencial cuadratica

def dq(z):
	evaluacion=fase.conjugate()
	for x in ceros:
		evaluacion=evaluacion*(z-x)/abs(z-x)
	for x in polos:
		evaluacion=evaluacion*((z-x)/abs(z-x))**-2
	for x in polosSimples:
		evaluacion=evaluacion*((z-x)/abs(z-x))**-1
	return evaluacion.conjugate()
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
def dqNot(z):
	evaluacion=fase.conjugate()
	for x in ceros:
            if x!= z:
		evaluacion=evaluacion*(z-x)/abs(z-x)
	for x in polos:
		evaluacion=evaluacion*((z-x)/abs(z-x))**-2
	for x in polosSimples:
		evaluacion=evaluacion*((z-x)/abs(z-x))**-1
	return evaluacion.conjugate()
###### generacion manual de los puntos
puntos = []
ceros = []
polos= []
polosSimples= []
tipoPunto='c'
#for x in np.linspace(-3,3.9,200):
#	for y in np.linspace(-3,3,20):
#		puntos.append(complex(-x,y))
#	puntos.append(complex(x,0)*rect(1,pi/3+0.1))
#	puntos.append(complex(x,0)*rect(1,-pi/3+0.1))


#parametros
maxreps=50000 #repeticiones maximas de la rutina de integracion
maxint=100
fase=rect(1,0.0) #fase de la diferencial
normav=0.001 #determina la norma del campo
normamax = abs(complex(10,10)) #determina el area dentro de la cual se integra el campo
t= np.linspace(0,1,5) #intervalo temporal
densidadPuntos=0.05 #determina cada cuantas repeticiones de la rutina de integracion se guarda un punto para su graficacion
colorLineas = 'k'
colorLineasCriticas = 'r'
colorLineasSillas= 'b'
anchoLineas = 0.2
actualizaClick = 0#determina si se redibuja la figura cada click
dibujaEjes = 0
cuadros=2
rotacion = 2*pi/cuadros
raizCubica= rect(1,2*pi/3)
lim = 30
esfera=0
####matplotlib
import matplotlib as mpl
#mpl.use("pgf")
import matplotlib.pyplot as plt
#from matplotlib import animation
from matplotlib.collections import LineCollection
#from mpl_toolkits import mplot3d as m3d
###Animacion


#Constantes e inicializacion de variables
ri=sqrt(1j)
ric=sqrt(-1j)
monodromia=0
ma=0
patch1=0
patch2=0
pause = False

########### funciones #########
def estereografica(z):
	return [(2*z.real/(1+z*z.conjugate())).real,(2*z.imag/(1+z*z.conjugate())).real,((1-z*z.conjugate())/(1+z*z.conjugate())).real]
def mon(inicio, fin):
	global monodromia
	global ma
	global patch1
	global patch2
	if inicio.real<=0 and inicio.imag>0 and fin.real<=0 and fin.imag<0:
		monodromia=(monodromia + 1)%2
		ma=(ma+1)%2
		#print('m')
		#print(str(patch1)+str(patch2))
	if inicio.real<=0 and inicio.imag<0 and fin.real<=0 and fin.imag>0:
		monodromia=(monodromia + 1)%2
		ma=(ma+1)%2
		#print('-m')
		#print(str(patch1)+str(patch2))
	if (-1)*inicio.real>abs(inicio.imag):
		if fin.real<0 and fin.imag>(-1)*fin.real:
			patch1=(patch1+1)%2
			ma=0
		if fin.real<0 and fin.imag<fin.real:
			patch2=(patch2+1)%2
			ma=0
	if inicio.imag<0 and inicio.imag<inicio.real:
		if (-1)*fin.real>abs(fin.imag):
			patch2=(patch2+1)%2
	if inicio.imag>0 and inicio.imag>(-1)*inicio.real:
		if (-1)*fin.real>abs(fin.imag):
			patch1=(patch1+1)%2
	if patch1 and patch2:
		#print("igu")
		patch1=0
		patch2=0
	return
def dist(z):
	if patch1 != patch2:
		if patch1:
			r=ri*(sqrt(z*(-1j)))
		if patch2:
			r=ric*(sqrt(z*(1j)))
	else:
		r=sqrt(z)
	return (-1)**(monodromia+ma)*r

def f(y,t):
	xi=y[0]
	yi=y[1]
	#c=1
	z=dist(dq(complex(xi,yi)))
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

def fi(y,t):
	xi=y[0]
	yi=y[1]
	#c=1
	z=-dist(dq(complex(xi,yi)))
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


def trayectp(z):
	global monodromia
	global patch1
	global patch2
	global ma
	monodromia, patch1, patch2, ma, = 0, 0, 0, 0
	norma=0.5
	#fin=z
	#x=(2*fin.real/(1+fin*fin.conjugate())).real
	#y=(2*fin.imag/(1+fin*fin.conjugate())).real
	#w=((1-fin*fin.conjugate())/(1+fin*fin.conjugate())).real
	if esfera==1:
		coord=np.array([estereografica(z)])
	else:
		coord=np.array([[z.real,z.imag]])
	#coord=np.array([[fin.real,fin.imag,0.0]])
	rep=0
	inicio=dq(z)
	fin=z
	if (-inicio.real)>abs(inicio.imag) and inicio.imag>0:
		patch1=1
	if (-inicio.real)>abs(inicio.imag) and inicio.imag<0:
		patch2=1
	inicio = z
	ultimo = z
	start = z
	while ( cercaPolo(fin) and cercaPoloS(fin) and norma < lim and rep<maxreps):
		sol = odeint(f,[inicio.real, inicio.imag],t,mxstep=maxint)
		fin= complex(sol[-1,0],sol[-1,1])
		mon(dq(inicio),dq(fin))
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
	global monodromia
	global patch1
	global patch2
	global ma
	monodromia, patch1, patch2, ma, = 0, 0, 0, 0
	norma=0.5
	#fin=z
	#x=(2*fin.real/(1+fin*fin.conjugate())).real
	#y=(2*fin.imag/(1+fin*fin.conjugate())).real
	#w=((1-fin*fin.conjugate())/(1+fin*fin.conjugate())).real
	if esfera==1:
		coord=np.array([estereografica(z)])
	else:
		coord=np.array([[z.real,z.imag]])
	#coord=np.array([[fin.real,fin.imag,0.0]])
	rep=0
	inicio=dq(z)
	fin=z
	if (-inicio.real)>abs(inicio.imag) and inicio.imag>0:
		patch1=1
	if (-inicio.real)>abs(inicio.imag) and inicio.imag<0:
		patch2=1
	inicio = z
	ultimo = z
	start = z
	while ( cercaPolo(fin) and cercaPoloS(fin) and norma < lim and rep<maxreps):
		sol = odeint(fi,[inicio.real, inicio.imag],t,mxstep=maxint)
		fin= complex(sol[-1,0],sol[-1,1])
		mon(dq(inicio),dq(fin))
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

#interfaz
#def onpress(event):
	#global fase
	#global pause
	#if event.key=='d':
		#fase = fase*rect(1,-0.1)
		#col.set_segments(trayectorias())
		#fig.canvas.draw()
		#print('detecto a')
		#print(fase)
	#if event.key=='a':
		#print(str(len(puntos))+' puntos seleccionados')
		#pause = True
		##anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
		#anim.save('basic_animation.mp4',bitrate=3000)
	#if event.key=='r':
		#print('redibujando...')
		#col.set_segments(coleccionTrayectorias)
		#fig.canvas.draw()
	#if event.key=='e':
		#plt.savefig('figura.pgf')

def trayectorias(listaPuntos):
	pool=Pool()
	tray=[]
	tray = pool.map(trayectoria,listaPuntos)
	pool.close()
	pool.join()
	print('!')
	return np.array(tray)
frame=0



#inicializacion de la rutina en paralelo
pool=Pool()
#resTray = pool.map(trayectoria,puntos)
#pool.close()
#pool.join()
#for curvaPrueb in resTray:
	#coleccionTrayectorias.append(list(zip(curvaPrueb[0],curvaPrueb[1],curvaPrueb[2])))
##for p in puntos:
	##curvaPrueb=trayectoria(p)
	##coleccionTrayectorias.append(list(zip(curvaPrueb[0],curvaPrueb[1],curvaPrueb[2])))
###########
#col.set_segments(coleccionTrayectorias)

#def calculaCuadro(i):
	#fase=rect(1,i*rotacion)
	#print('cuadro '+str(i))
	#return trayectorias(puntos)
#pool=Pool()
#anim=pool.map(calculaCuadro, range(0,cuadros))
#pool.close()
#pool.join()

#sys.stdout.write("\n Escribiendo resultados\n")
#resultado=np.array(anim)
#np.save("anim.npy",resultado)

#inicializacion de la figura
fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111, xlim=(-3,3), ylim=(-3,3), autoscale_on=False)
offs= (0.0, 0.0)
ax.get_xaxis().set_visible(dibujaEjes)
ax.get_yaxis().set_visible(dibujaEjes)
coleccionTrayectorias=[]
col = LineCollection(coleccionTrayectorias,offsets=offs,color=colorLineas,linewidths=anchoLineas,antialiaseds=True)
colCriticas = LineCollection(coleccionTrayectorias,offsets=offs,color=colorLineasCriticas,linewidths=anchoLineas,antialiaseds=True)
colSillas= LineCollection(coleccionTrayectorias,offsets=offs,color=colorLineasSillas,linewidths=anchoLineas,antialiaseds=True)
ax.add_collection(col, autolim=True)
ax.add_collection(colCriticas, autolim=True)
ax.add_collection(colSillas, autolim=False)
def onpress(event):
	global fase
        global puntos
        global tipoPunto
	if event.key=='d':
		fase = fase*rect(1,-0.1)
		col.set_segments(trayectorias(puntos))
		fig.canvas.draw()
		print('detecto d')
		print(fase)
	if event.key=='a':
		fase = fase*rect(1,0.1)
		col.set_segments(trayectorias(puntos))
		fig.canvas.draw()
		print('detecto a')
		print(fase)
	if event.key=='r':
		print('redibujando...')
		col.set_segments(trayectorias(puntos))
		fig.canvas.draw()
	if event.key=='R':
                puntosExtra=[]
		print('redibujando criticas...')
                for x in ceros:
                    for i  in range(3):
                        z=x+(raizCubica**i*dqNot(x)**(1/3.0))*0.1
                        puntosExtra.append(z) 
		colCriticas.set_segments(trayectorias(puntosExtra))
		fig.canvas.draw()
	if event.key=='e':
		plt.savefig('figura.pgf')
	if event.key=='C':
            puntos=[]
	if event.key=='c':
            tipoPunto=event.key
	if event.key=='p':
            tipoPunto=event.key
	if event.key=='u':
            tipoPunto=event.key
        if event.key=='m':
            trayectoriasSilla=[]
            puntosSilla=[]
	    print('redibujando sillas...')
            faseOld = fase
            for x in combinations(ceros,2):
                z=x[0]
                w=x[1]
                faseS=faseSilla(z,w)
                #fase=faseS.conjugate()
                fase=faseS**2
                puntosSilla.append(z-faseS**2*0.01)
                #tray=trayectoria(z+faseS*0.1)
                #trayectoriasSilla.append(list(zip(tray[0],tray[1])))
                #trayectoriasSilla.append(trayectoria(z-faseS**(1/3.0)*0.3))
                for i in range(3):
                    puntoInicial=z+raizCubica**i*dqNot(z)**(1.0/3.0)*0.1
                    trayectoriasSilla.append(trayectoria(puntoInicial))
                    #plt.plot([puntoInicial.real],[puntoInicial.imag],'bo')
                print("fase  %g+i%g, z %g+i%g" % (fase.real,fase.imag,z.real,z.imag))
	    #colSillas.set_segments(trayectorias(puntosSilla))
            #fase = faseOld
            colSillas.set_segments(trayectoriasSilla)
	    fig.canvas.draw()
                





def faseSilla(z,w):
    f = lambda x: sqrt(DQ(x))
    faseSilla = intfdz(f,z,w)
    return faseSilla/abs(faseSilla)
#	periodoQI = quad(raizI, 0, 1)
#	periodoQR = quad(raizR, 0, 1)


def onclick(event):
	global puntos
        global tipoPunto
	if event.button == 1:
		puntos.append(complex(event.xdata,event.ydata))
		print('P')
		if actualizaClick ==1:
			col.set_segments(coleccionTrayectorias)
			fig.canvas.draw()
	if event.button == 3:
                if tipoPunto=='c':
	        	global ceros 
	        	ceros.append(complex(event.xdata,event.ydata))
	        	plt.plot([event.xdata],[event.ydata],'ro')
	        	print('C')
	        	usarCeros = True
                if tipoPunto=='p':
	        	polos.append(complex(event.xdata,event.ydata))
	        	plt.plot([event.xdata],[event.ydata],'go')
	        	print('C')
                if tipoPunto=='u':
	        	global polosSimples 
	        	polosSimples.append(complex(event.xdata,event.ydata))
	        	plt.plot([event.xdata],[event.ydata],'bo')
	        	print('s')
        if event.button == 2:
		polos.append(complex(event.xdata,event.ydata))
		plt.plot([event.xdata],[event.ydata],'go')
		print('P')

fig.canvas.mpl_connect('key_press_event', onpress)
fig.canvas.mpl_connect('button_press_event', onclick)
##mostrar el ambiente

#def animate(i):
	#global fase, frame
	#if pause:
		#fase = fase*rect(1,rotacion)	
		#col.set_segments(trayectorias(puntos))
		#frame = frame+1
		#print(' Cuadro '+str(frame)+' de 50 completado')
	#return col,

##anim = animation.FuncAnimation(fig, animate, init_func=init,
                               ##frames=200, interval=50, blit=True)
#pool=Pool(1)
#resTray = pool.map(trayectoria,puntos)
#pool.close()
#pool.join()
#lista=[]
#for x in puntos:
	#lista.append(trayectoria(x))

#inicializacion de la figura
#np.save("hoyo.npy",np.array(resTray))
#fig = plt.figure()
#if esfera==1:
	#ax = fig.add_subplot(111, xlim=(-1,1), ylim=(-1,1), zlim=(-1,1), projection='3d')
	#col = m3d.art3d.Line3DCollection(resTray,color=colorLineas,linewidths=anchoLineas,antialiaseds=True)
	#ax.add_collection3d(col)
#else:
#plt.savefig('figure.pgf')

def dqq(z):
	return (fase.conjugate()*(z)*(z-complex(1,1))*(z-complex(1,-1)))
	#return (fase*(1-z)*(z+1)*(z-complex(1,1))*(z-complex(1,-1))*(z-complex(-1,1))*(z-complex(-1,-1))*(z-complex(0,-1.0))*(z-complex(0,1.0)))
fases = np.linspace(0,2*pi,500)
#gama = lambda x: complex(x,-x)
#raizI = lambda x: (complex(1,-1)*sqrt(dqq(gama(x)))).imag
#raizR = lambda x: (complex(1,1)*sqrt(dq(gama(x)))).real
#for alfa in fases:
#	fase=rect(1,alfa)
#	#periodoR = romberg(raiz, 0, 1, show=True)
#	periodoQI = quad(raizI, 0, 1)
#	periodoQR = quad(raizR, 0, 1)
#	print("periodo: %g + i%g fase: %g" % (periodoQR[0],periodoQI[0],alfa))
#	if 10**-3>abs(periodoQI[0]):
#		print("la fase buscada es %g+i%g" % (fase.real,fase.imag))
#		break
plt.show()
