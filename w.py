from __future__ import division
import os
import numpy as np
import matplotlib.pyplot as plt
#import h5py
import math
import scipy.integrate as integrate
import cmath
from matplotlib.pyplot import cm
from scipy.integrate import dblquad
from sympy import *
import traceback
import time


try:
	start_time=time.time()
	L=12;
	nx=2;
	ny=2;

	#plot limits
	angle=16 #nPi/angle, loops go from 0 to 2*angle; determines how fine directions selection is
	sigma=5
	limit=L

	step=1
	X_max=L+1
	X_min=0
	Y_max=X_max
	Y_min=X_min
	total=((X_max-X_min+step-1)*(Y_max-Y_min+step-1))/step**2


	#calculation constants
	norm_constant= 1/(sigma*math.sqrt(2*math.pi))
	k2= np.array([-1,1])#arguments of the wavefunction
	k1= np.array([1,1])
	
	#definitions
	def ko(theta):#posistion in k space determined by angle
		return np.array([math.cos(math.pi*theta/angle),math.sin(math.pi*theta/angle)])
	def r(x,y):#integration variable in r space
		return np.array([x,y])

	# def wf(xo,yo):#wave function pb
	# 	return math.sqrt(2/L)*math.sin(nx*math.pi*xo/L)*math.sin(ny*math.pi*yo/L)
	# def wf(xo,yo):#wave function of th epure momemntum state
	# 	return math.cos(np.dot(k1,np.array([xo,yo])))+math.cos(np.dot(k2,np.array([xo,yo])))
	def wf(xo,yo):#other pure momentum
		return np.conj(cmath.exp(1j*np.dot(k1,np.array([xo,yo])))+cmath.exp(-1j*np.dot(k2,np.array([xo,yo]))))
	
	def gauss(theta,x,y,xo,yo):#coherent state conjugated
		return norm_constant*cmath.exp((-(x-xo)**2-(y-yo)**2)/(4*sigma**2))*cmath.exp( 1j*np.dot(ko(theta),r(x,y)))

	def cint(theta,xo,yo):

		a= dblquad(lambda y, x: np.real(gauss(theta,x,y,xo,yo)*wf(xo,yo)), 0, limit, lambda x: 0, lambda x: limit)#integration of real
		b= dblquad(lambda y, x: np.imag(gauss(theta,x,y,xo,yo)*wf(xo,yo)), 0, limit, lambda x: 0, lambda x: limit)#integration of imaginary
		return (np.linalg.norm(a[0]+1j*b[0]))**2

	def kazsa(xo,yo,n):
		bx = plt.axes()
		if n==1:
			n=np.amax(z)
		if n==0:
			n=np.amax(z[xo,yo])

		for t in range (0,2*angle):
			if z[xo,yo,t]<10**-30:
				z[xo,yo,t]=0
			bx.arrow(xo,yo,
				z[xo,yo,t]*math.cos(math.pi*t/angle)/(n),
				z[xo,yo,t]*math.sin(math.pi*t/angle)/(n),
				head_width=.05, head_length=0.1, fc='k', ec='k')

		bx.set_xlim(X_min-step,X_max+step)
		bx.set_ylim(Y_min-step,Y_max+step)
		plt.show()
	


	os.system('clear')
	counter=0
	z=np.zeros([X_max,Y_max,2*angle])#array that stores |<gauss|wf>|^2 for each k(theta)

	#calcutions and also makes arrow objects
	ax = plt.axes()
	for xo in range (X_min,X_max,step):#where the gauss is centered?
		for yo in range (Y_min,Y_max,step):
			for t in range (0,2*angle):#angle that determines posistion in k space
				c=dblquad(lambda y, x: gauss(t,x,y,xo,yo)*wf(xo,yo),0, limit, lambda x: 0, lambda x: limit)
				z[xo,yo,t]=(np.linalg.norm(c[0]))**2#magnitude squared of the complex number from integrating
				# z[xo,yo,t]=cint(t,xo,yo)
			os.system('clear')
			counter +=1#progress bar counter
			print (xo,yo), "\t %2d%%" %(100*counter/total)
			
	os.system('clear')

	for xo in range (X_min,X_max,step):
		for yo in range (Y_min,Y_max,step):
			for t in range (0,2*angle):
				if z[xo,yo,t]<10**-30:
					z[xo,yo,t]=0

				ax.arrow(xo,yo,
					z[xo,yo,t]*math.cos(math.pi*t/angle)/(2*np.amax(z)),
					z[xo,yo,t]*math.sin(math.pi*t/angle)/(2*np.amax(z)),
					head_width=.05, head_length=0.1, fc='k', ec='k')






	# displays arrows
	os.system('play --no-show-progress --null --channels 1 synth %s sine %f' % ( .5, 523))
	os.system('play --no-show-progress --null --channels 1 synth %s sine %f' % ( .5, 880))
	os.system('play --no-show-progress --null --channels 1 synth %s sine %f' % ( .5, 699))

	print("Computation time: %.7s seconds" % (time.time() - start_time))

	ax.set_xlim(X_min-step,X_max+step)
	ax.set_ylim(Y_min-step,Y_max+step)
	plt.show()



except Exception as x:
	print("Computation time: %s seconds" % (time.time() - start_time))
	os.system('play --no-show-progress --null --channels 1 synth %s sine %f' % ( .5, 523))
	os.system('play --no-show-progress --null --channels 1 synth %s sine %f' % ( .5, 493.88))
	os.system('play --no-show-progress --null --channels 1 synth %s sine %f' % ( 1, 466.16))
	traceback.print_exc()
	print "whoops..."