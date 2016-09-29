import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import math
import scipy.integrate as integrate
import cmath
from decimal import Decimal as D

from matplotlib.pyplot import cm

from scipy.integrate import quad
from plotly.tools import FigureFactory as FF

pi=math.pi

sigma=2#deltak/k~35%?

angle=16 #nPi/angle, loops go from 0 to 2*angle; determines how fine directions selection is
Rrange=10
Crange=10
step=10
a=1
b=.8



Mag_at_angle = np.zeros(shape=(Rrange,Crange,angle*2+1)) #the value after you process a direction
Mag_in_x=np.zeros(shape=(Rrange,Crange,angle*2+1),dtype=object)#component x of mag_at_angle
Mag_in_y=np.zeros(shape=(Rrange,Crange,angle*2+1),dtype=object)#component y of mag_at_angle
resultant_mag=np.zeros(shape=(Rrange,Crange,4),dtype=object)#sum of Mags in each quadrant which correspond to test vectors[1,1],[-1,1],etc

def holy_hand_gernade(x,y,theta,xo,yo):	#Husimi process 

	r=np.array([x,y])#position in space
	ro=np.array([xo,yo])#used in summation loop

	#Acos(k1 dot r)+Bcos(k2 dot r)

	k1=np.array([1,1])
	k2=np.array([-1,1])
	ko=np.array([math.cos(math.pi*theta/angle),math.sin(math.pi*theta/angle)])

	r_minus_ro_squared=(r[0]-ro[0])**2 + (r[1]-ro[1])**2 #|[(x-xo),(y-yo)]|^2=(r-ro)^2. appears in the exponential of the coherent state
	ko_dot_ro= ko[0]*ro[0]+ko[1]*ro[1] #appears in the exponential of the coherent state

	holy= np.conjugate(a*np.cos(np.dot(k1,r))+b*np.cos(np.dot(k2,r)))#complex conjugate of the wavefunction

	grail = (1/(2*math.pi*sigma))*math.exp(r_minus_ro_squared/(4*sigma**2))*cmath.exp(1j*ko_dot_ro)#gaussian

	holygrail= abs(holy*grail)# psi* gauss

	return holygrail**2


#summation
temp=0
for row in range (0,Rrange,step):#position
	for column in range(0,Crange,step):#position
		for theta in range (0,2*angle):#r direction
			for xo in range(0,5):#xo
				for yo in range(0,5):#yo
					temp=temp+holy_hand_gernade(row,column,theta,xo,yo)

			Mag_at_angle[row,column,theta]=temp
			Mag_in_x[row,column,theta]=Mag_at_angle[row,column,theta]*np.cos((math.pi*theta/angle))#split up the mag so arrows can be plotted
			Mag_in_y[row,column,theta]=Mag_at_angle[row,column,theta]*np.sin((math.pi*theta/angle))
			temp=0

ax = plt.axes()
for row in range (0,Rrange,step):
	for column in range (0,Crange,step):
		for theta in range (0,2*angle):
			print Mag_in_x[row,column,theta],Mag_in_y[row,column,theta]
			ax.arrow(row,column,Mag_in_x[row,column,theta],Mag_in_y[row,column,theta], head_width=0, head_length=0,linewidth=.5)


# #take the dot prodcts of the husumi vectors (that lie between 0,pi/2; pi/2,pi; pi,3pi/2; 3pi/2,2pi) with the test vectors ( 1,1; -1,1; -1,-1; 1,-1)  
# for row in range (0,Rrange,step):
# 	for column in range(0,Crange,step):#add some ifs
# 		for theta in range (0,angle/2):
# 			resultant_mag[row,column,0]=resultant_mag[row,column,0]+np.dot(np.array((Mag_in_x[row,column,theta],Mag_in_y[row,column,theta])),np.array((1,1)))
# 		for theta in range (angle/2,angle):

# 			resultant_mag[row,column,1]=resultant_mag[row,column,1]+np.dot(np.array((Mag_in_x[row,column,theta],Mag_in_y[row,column,theta])),np.array((-1,1)))

# 		for theta in range (angle,3*angle/2):
# 			resultant_mag[row,column,2]=resultant_mag[row,column,2]+np.dot(np.array((Mag_in_x[row,column,theta],Mag_in_y[row,column,theta])),np.array((-1,-1)))
# 		for theta in range (3/2*angle,2*angle):
# 			resultant_mag[row,column,3]=resultant_mag[row,column,3]+np.dot(np.array((Mag_in_x[row,column,theta],Mag_in_y[row,column,theta])),np.array((1,-1)))


# #subtract the dot products from the test vectors and then multiply it times the test vectors?
# for row in range (0,Rrange,step):
# 	for column in range(0,Crange,step):

# 		resultant_mag[row,column,0]=np.subtract(np.array((1,1)),((resultant_mag[row,column,0])*np.array((1,1))))
# 		resultant_mag[row,column,1]=np.subtract(np.array((-1,1)),((resultant_mag[row,column,1])*np.array((-1,1))))
# 		resultant_mag[row,column,2]=np.subtract(np.array((-1,-1)),((resultant_mag[row,column,2])*np.array((-1,-1))))
# 		resultant_mag[row,column,3]=np.subtract(np.array((1,-1)),((resultant_mag[row,column,3])*np.array((1,-1))))

# #plot that shit
# #ax.arrow(x position, y position, x length, y length)			
# for row in range (0,Rrange,step):
# 	for column in range(0,Crange,step):
# 		print row, column
# 		print resultant_mag[row,column,:]
# 		ax = plt.axes()
# 		ax.arrow(row, column, resultant_mag[row,column,0][0],resultant_mag[row,column,0][1], head_width=1, head_length=1)
# 		ax.arrow(row, column, resultant_mag[row,column,1][0],resultant_mag[row,column,1][1], head_width=1, head_length=1)
# 		ax.arrow(row, column, resultant_mag[row,column,2][0],resultant_mag[row,column,2][1], head_width=1, head_length=1)
# 		ax.arrow(row, column, resultant_mag[row,column,3][0],resultant_mag[row,column,3][1], head_width=1, head_length=1)

plt.ylim([-1,Crange])
plt.xlim([-1,Rrange])
plt.show()
