import numpy as np
import math
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.lines as lin
from minuit2 import Minuit2 as minuit
from numpy  import *

castor_inner_octant_radius = 41
beampipe_r = 57/2 #in mm
beampipe_x = 0 #we shift beampipe instead of castor halfs
beampipe_y = 0

def distance_to_beampipe( x , y ):
   "calculate distance for any given point to beam pipe outer circle"
   global beampipe_x
   global beampipe_y
   global beampipe_r
   r = math.sqrt((x-beampipe_x)**2 + (y-beampipe_y)**2)
   return math.fabs(beampipe_r - r);

class castor_half():
    def __init__(self, name,sensor_angles,sensor_pos,rnew):
        self.nsensors = len(sensor_angles)
        assert  self.nsensors == len(rnew), '#angles must be #rnew'
        assert  self.nsensors == len(sensor_pos), '#angles must be #positions'
        if name.find("far") != -1:
            self.isFarHalf = True
        elif name.find("near") != -1:
            self.isFarHalf = False
        else:
            sys.stderr.write("Please choose name with either \"near\" or \"far\" in it\n")
            exit(1)
        self.sensor_angles = sensor_angles #x direction is 0. counter-clockwise
        self.sensor_pos = sensor_pos #in ideal geometry
        self.rnew = rnew
        self.verbosity = 1

    def setVerbosity(self, verbosity):
        #0: no output, 1:essential results, 2:everything
        self.verbosity = verbosity

    def fit_pos(self):
        def f(x,y):
            #chi2 function. returns distance to beam pipe circle for all vector(sensorpos)-vector(r) for given (x,y) shift
            chi2 = 0
            sigma = 1
            for i in range(0,self.nsensors):
                pos = array(self.sensor_pos[i])
                angle = self.sensor_angles[i]
                r = self.rnew[i]
                pointingat = pos - array([r * math.sin(radians(angle)), r * math.cos(radians(angle))])
                chi2 += distance_to_beampipe(pointingat[0]+x,pointingat[1]+y) ** 2 / sigma ** 2
            return chi2
        m = minuit(f)
        if verbosity > 1 : m.printMode = 1
        m.migrad()
        self.x = m.values["x"]
        self.y = m.values["y"]
        chi2 = m.fval / self.nsensors
        if verbosity > 0 :print "far half: " if self.isFarHalf else "near half: ","x={0:.2f} y={1:.2f} chi2={2:.2f}".format(self.x,self.y,chi2)
        return

verbosity = 1 
far_angles = [22.5,-67.5]
far_pos = []
for angle in far_angles:
    far_pos.append([math.cos(math.radians(angle))*castor_inner_octant_radius , math.sin(math.radians(angle))*castor_inner_octant_radius])
if verbosity > 0:
    print "\nPosition of far side sensors is: "
    for pos in far_pos: print "{0[0]:.3f},{0[1]:.3f}".format(pos)

far_r_old = [8.68439,20.2883] #top,bottom
far_r_old_error = [0,0.00508586] #top,bottom
farside_old = castor_half("farside_old",far_angles,far_pos,far_r_old)
farside_old.setVerbosity(verbosity)
farside_old.fit_pos()

far_r_new = [12.6556,22.2788] #top,bottom
far_r_new_error = [0,0.00172011] #top,bottom
farside_new = castor_half("farside_new",far_angles,far_pos,far_r_new)
farside_new.setVerbosity(verbosity)
farside_new.fit_pos()


near_angles = [180-22.5,180+67.5]
near_pos = []
for angle in near_angles:
    near_pos.append([math.cos(math.radians(angle))*castor_inner_octant_radius , math.sin(math.radians(angle))*castor_inner_octant_radius])
if verbosity > 0:
    print "\nPosition of near side sensors is: "
    for pos in near_pos: print "{0[0]:.3f},{0[1]:.3f}".format(pos)

near_r_old = [15.2248,25.8462] #top,bottom
near_r_old_error = [0.00119269,0.0232433] #top,bottom
nearside_old = castor_half("nearside_old",near_angles,near_pos,near_r_old)
nearside_old.setVerbosity(verbosity)
nearside_old.fit_pos()

near_r_new = [32.9261,32.6869] #top,bottom
near_r_new_error = [2.44163e-15,0.00781898] #top,bottom
nearside_new = castor_half("nearside_new",near_angles,near_pos,near_r_new)
nearside_new.setVerbosity(verbosity)
nearside_new.fit_pos()



def printShift(old,new):
    dx = new.x-old.x
    dy = new.y-old.y
    dr = math.sqrt(dx**2+dy**2)
    print "Shift due to magnetic field (","far  half" if old.isFarHalf else "near half", "): dx={0:.3f} dy={1:.3f} dr={2:.3f}".format(dx,dy,dr)

print("\n\n")
printShift(farside_old,farside_new)
printShift(nearside_old,nearside_new)

def draw(old,new):

    circle1=plt.Circle((beampipe_x,beampipe_y),beampipe_r,color='0.8')
    fig = plt.figure(figsize=[4,4])
    fig.gca().add_artist(circle1)
                   
    assert old.nsensors == new.nsensors
    data_old = []
    for i in range(0,old.nsensors):
        data_old.append([old.sensor_pos[i][0]+old.x,old.sensor_pos[i][1]+old.y])
        
    if verbosity>1: print "drawing sensors:" << data_old ,"\n", [row[0] for row in data_old], [row[1] for row in data_old]
    plt.scatter([row[0] for row in data_old], [row[1] for row in data_old], s=5, alpha=0.5, color='r') #s=area
                    
    l = lin([0,100],[0,10])                                    
    plt.axis.add_line(l)

    plt.xlabel('x [mm]')
    plt.ylabel('y [mm]')
    plt.title('IR Sensor [Jan 07 -> Jan 17]')
    #plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
    plt.axis([-50, 50, -50, 50])
    
    plt.savefig("test.png")

draw(farside_old,farside_new)

# e = np.e
# X, Y = np.meshgrid(np.linspace(0, castor_inner_octant_radius, 100), np.linspace(0, castor_inner_octant_radius, 100))
# F = X ** Y
# G = Y ** X

# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1)
# plt.axis([-100, 100, -100, 100])
# plt.contour(X, Y, (F - G), [0])
# plt.plot([e], [e], 'g.', markersize=20.0)
# #plt.show()
# plt.savefig("test.png",bbox_inches=0.100)
