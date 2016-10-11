import numpy as np
import matplotlib.pyplot as plt
#import h5py
import math
import scipy.integrate as integrate
import cmath
from matplotlib.pyplot import cm
from scipy.integrate import dblquad



alpha=1
beta=0.8
sigma=3
norm_constant= 1/(math.sqrt(2*math.pi))

angle=10 #nPi/angle, loops go from 0 to 2*angle; determines how fine directions selection is
X_max=2
X_min=1
Y_max=X_max
Y_min=X_min

riemann_delta_x=0.25
riemann_delta_y=riemann_delta_x
riemann_interval_x=4 #multiples of four
riemann_interval_y=riemann_interval_x



k_o= lambda theta: np.array([math.cos(math.pi*theta/angle),math.sin(math.pi*theta/angle)])
k_1= np.array([1,1])
k_2= np.array([1,-1])
r= lambda x,y: np.array([x,y])



def gauss_wf(xo,yo,variables):

	x=variables[0]
	y=variables[1]
	theta=variables[2]

	momentum_state= alpha*math.cos(k_1[0]*r(x,y)[0]+k_1[1]*r(x,y)[1])+beta*math.sin(k_2[0]*r(x,y)[0]+k_2[1]*r(x,y)[1])
	
	r_minus_ro= math.sqrt((x-xo)**2+(y-yo)**2)

	gauss_conjugate=norm_constant*cmath.exp((r_minus_ro/(4*sigma*2))+1j*(k_o(theta)[0]*r(x,y)[0]+k_o(theta)[1]*r(x,y)[1]))#only one with theta dependence

	return gauss_conjugate*momentum_state

def midpoint_riemann_sum(x,y,variables):
	S=0 #Sum
	
	for xo in range(0,int(riemann_interval_x*(1/riemann_delta_x))):
		for yo in range(0,int(riemann_interval_y*(1/riemann_delta_y))):

			x_1=xo*riemann_delta_x +riemann_delta_x/2
			y_1=yo*riemann_delta_y+riemann_delta_y/2

			S=S+gauss_wf(x_1,y_1,variables)*riemann_delta_x*riemann_delta_y
	
	S_magnitude=math.sqrt(S.real**2+S.imag**2)
	print S.real**2+S.imag**2#accidently build a unit circle? all magnitudes are the same...

	return S_magnitude**2


ax = plt.axes()
for x in range (X_min,X_max):
	for y in range (Y_min,Y_max):
		for theta in range(0,2*angle):

			variables=np.array([x,y,theta])
			S=midpoint_riemann_sum(x,y,variables)
			ax.arrow(x,y,S*math.cos(math.pi*theta/angle),S*math.sin(math.pi*theta/angle), head_width=0, head_length=0,linewidth=.5)

			print (x,y,theta),"\t",S

plt.ylim([-10,10])
plt.xlim([-10,10])
#plt.show()