from __future__ import division #always float divisions for python version < 3.0
import numpy as np
import math
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.lines as lin
import matplotlib.patches as pat
import matplotlib.legend as pyleg
import matplotlib.patches as patches
from matplotlib.path import Path
from minuit2 import Minuit2 as minuit #minuit2.so needed. Please compile yourself
from numpy  import *
from scipy.interpolate import UnivariateSpline
from scipy import arange, array, exp

castor_inner_octant_radius = 40.6507 #cos(22.5) * 44
beampipe_r = (57+0.5)/2 #in mm #+white paper = 0.5mm
beampipe_x = 0#-0.7 #hauke's value
beampipe_y = 0#-1.4

def rotatePoint(centerPoint,point,angle):
    """Rotates a point around another centerPoint. Angle is in degrees.
    Rotation is counter-clockwise"""
    angle = math.radians(angle)
    temp_point = point[0]-centerPoint[0] , point[1]-centerPoint[1]
    temp_point = ( temp_point[0]*math.cos(angle)-temp_point[1]*math.sin(angle) , temp_point[0]*math.sin(angle)+temp_point[1]*math.cos(angle))
    temp_point = temp_point[0]+centerPoint[0] , temp_point[1]+centerPoint[1]
    return temp_point

def distance_to_beampipe( x , y , xe=0 , ye=0 ):
   """calculate distance for any given point to beam pipe outer circle"""
   global beampipe_x
   global beampipe_y
   global beampipe_r
   r = math.sqrt((x-beampipe_x)**2 + (y-beampipe_y)**2)
   distance = math.fabs(beampipe_r - r);
   ddisdx = 1 / 2 / r * 2 * (x-beampipe_x)
   ddisdy = 1 / 2 / r * 2 * (y-beampipe_y)
   error = math.sqrt( xe**2 * ddisdx**2 + ye**2 * ddisdy**2)
   return distance , error

class castor_half():
    """Class that takes for a castor half that takes sensor information and can shift/fit position so that sensors point to beam pipe"""
    def __init__(self, name,sensors):
        self.nsensors = len(sensors)
        if name.find("far") != -1:
            self.isFarHalf = True
        elif name.find("near") != -1:
            self.isFarHalf = False
        else:
            sys.stderr.write("Please choose name with either \"near\" or \"far\" in it\n")
            exit(1)
        self.sensors = sensors #x direction is 0. counter-clockwise
        self.verbosity = 1

    def setVerbosity(self, verbosity):
        #0: no output, 1:essential results, 2:everything
        self.verbosity = verbosity

    def fit_pos(self):
        def f(x,y):
            #chi2 function. returns distance to beam pipe circle for all vector(sensorpos)-vector(r) for given (x,y) shift
            chi2 = 0
            for iSen in self.sensors:
                r = iSen.GetCalibratedDist()
                pointingat = iSen.pos - array([r * math.cos(radians(iSen.angle)), r * math.sin(radians(iSen.angle))])

                error_r = sqrt(1**2 + iSen.GetDistError()**2) #1mm sys + stat error
                error_theta = radians(1)

                dpxdr     = -math.cos(radians(iSen.angle))
                dpxdtheta = r * math.sin(radians(iSen.angle))
                dpydr     = -math.sin(radians(iSen.angle))
                dpydtheta = -r * math.cos(radians(iSen.angle))

                xe = sqrt( error_r**2 * dpxdr**2 + error_theta**2 * dpxdtheta**2)#x and y is initial sensor positions from drawings. maybe no uncertainty
                ye = sqrt( error_r**2 * dpydr**2 + error_theta**2 * dpydtheta**2)#x and y is initial sensor positions from drawings. maybe no uncertainty

                delta,sigma = distance_to_beampipe(pointingat[0]+x,pointingat[1]+y, xe, ye)
                if self.verbosity>1: print "xe=",xe,"ye=",ye, " --> delta=",delta,"sigma=",sigma
                chi2 += delta ** 2 / sigma ** 2
                
            return chi2
        m = minuit(f)
        if verbosity > 1 : m.printMode = 1
        m.migrad()
        m.minos()
        self.x = m.values["x"]
        self.y = m.values["y"]
        self.xeu = m.merrors["x", 1]
        self.yeu = m.merrors["y", 1]
        self.xel = m.merrors["x", -1]
        self.yel = m.merrors["y", -1]
        chi2 = m.fval / self.nsensors
        if verbosity > 0 :print "far half: " if self.isFarHalf else "near half: "," fitted position (x,y)=({0:.2f}{1:+.2f}{2:+.2f},{3:.2f}{4:+.2f}{5:+.2f}) with chi2={6:.2f}".format(self.x,self.xeu,self.xel,self.y,self.yeu,self.yel,chi2)
        return

class sensor():
    """stores information about position and angle about each sensor"""        
    def __init__(self, pos, angle):
        """angle counter clock-wise in deg [-180,180] starting from 3pm looking from IP"""
        self.pos = pos
        self.angle = angle
        self.meas_r = 0
        self.meas_r_err = 0
        self.cal_meas = []
        self.cal_true = []
        self.verbosity = 1

        assert len(self.pos) == 2, 'pos must be provided as list [x,y]'
        assert -180 < self.angle < 180, 'error with angle for sensor'

        if self.verbosity > 0:
            print "             nominal position (x,y)=({0[0]:.2f},{0[1]:.2f})".format(self.pos)
    @classmethod
    def fromsensor(cls, sensor): #copy constructor
        new = cls(sensor.pos,sensor.angle)
        new.SetCalibrationData(sensor.cal_meas,sensor.cal_true)
        new.SetDist(sensor.meas_r,sensor.meas_r_err)
        return new


    def SetCalibrationData(self,meas,true):
        """set measured calibration data. like [-1.2, 10.1, 20.4] and [0, 10, 20]"""
        self.cal_meas = meas
        self.cal_true = true
        assert len(self.cal_meas) == len(self.cal_true), 'error in calibration data'
        self.cal_spline = UnivariateSpline(self.cal_meas,self.cal_true,k=2)
        
    def GetCalibratedDist(self): #apply zero shift and linearity from measurement
        """apply calibration and return distance. if outside of calibration data it will be extrapolated"""
        if (self.cal_spline):
            return self.cal_spline(self.meas_r)[0]
        else:
            if self.verbosity > 1:
                print "No calibration data set. Using uncalibrated"
            return self.meas_r
    
    def GetDist(self):
        return self.meas_r

    def GetDistError(self):
        return self.meas_r_err

    def SetDist(self,r,r_err):
        self.meas_r = r
        self.meas_r_err = r_err
        
####FITTING####
verbosity = 1
offcenter = -6.35 #mm distance because the sending device is not in the middle of the housing

####FITTING FAR SIDE####
#parameters are [angle, off-center unrotated laser y-direction in mm, zeroshift means distance when touching the metal in mm]
if verbosity > 0:
    print "\nPositions of -far side- sensors are: "

sensor_fartop = sensor(rotatePoint([0,0],[castor_inner_octant_radius,offcenter], 180-67.5), 180-67.5) #offcenter applied to y position unrotated. this is correct (clockwise)
sensor_farbot = sensor(rotatePoint([0,0],[castor_inner_octant_radius,offcenter], -180+22.5), -180+22.5)
                     
sensor_fartop.SetDist(8.68439-1,0)
sensor_farbot.SetDist(20.2883-1,0.00508586)

sensor_fartop.SetCalibrationData([-2.,9.8,19.1], [0,10,20])
sensor_farbot.SetCalibrationData([0.5,10.1,20.2], [0,10,20])

farside_old = castor_half("farside_old",[sensor_fartop,sensor_farbot])
farside_old.setVerbosity(verbosity)
print "Before B field"
farside_old.fit_pos()

sensor_fartop = sensor.fromsensor(sensor_fartop)
sensor_farbot = sensor.fromsensor(sensor_farbot)
sensor_fartop.SetDist(12.6556-1,0)
sensor_farbot.SetDist(22.2788-1,0.00172011)

farside_new = castor_half("farside_new",[sensor_fartop,sensor_farbot])
farside_new.setVerbosity(verbosity)
print "After B field"
farside_new.fit_pos()

####FITTING NEAR SIDE####
if verbosity > 0:
    print "\nPositions of -near side- sensors are: "

sensor_neartop = sensor(rotatePoint([0,0],[castor_inner_octant_radius,offcenter], 67.5), 67.5) #offcenter applied to y position unrotated. this is correct (clockwise)
sensor_nearbot = sensor(rotatePoint([0,0],[castor_inner_octant_radius,offcenter], -22.5), -22.5)

sensor_neartop.SetCalibrationData([-2.,9.8,19.1], [0,10,20])
sensor_nearbot.SetCalibrationData([0.5,10.1,20.2], [0,10,20])

sensor_neartop.SetDist(15.2248-2,0.00119269)
sensor_nearbot.SetDist(25.8462-2,0.0232433)


nearside_old = castor_half("nearside_old",[sensor_neartop,sensor_nearbot])
nearside_old.setVerbosity(verbosity)
print "Before B field"
nearside_old.fit_pos()

sensor_neartop = sensor.fromsensor(sensor_neartop)
sensor_nearbot = sensor.fromsensor(sensor_nearbot)
sensor_neartop.SetDist(32.9261-2,2.44163e-15)
sensor_nearbot.SetDist(32.6869-2,0.00781898)

nearside_new = castor_half("nearside_new",[sensor_neartop,sensor_nearbot])
nearside_new.setVerbosity(verbosity)

print "After B field"
nearside_new.fit_pos()

def printShift(old,new):
    dx = new.x-old.x
    dy = new.y-old.y
    dr = math.sqrt(dx**2+dy**2)
    print "Shift due to magnetic field ({half}): dx={dx:.2f} dy={dy:.2f} dr={dr:.2f}".format(half="far  half" if old.isFarHalf else "near half",dx=dx,dy=dy,dr=dr)

print("\n\n")
printShift(farside_old,farside_new)
printShift(nearside_old,nearside_new)



#####DRAWING####
fig = plt.figure(figsize=[8,8])
ax = fig.gca()

circle1=plt.Circle((beampipe_x,beampipe_y),beampipe_r,color='black',fill=False,label="beampipe",lw=2)
ax.add_artist(circle1)

leglables=["beampipe"]
legpointers=[circle1]

def draw(fig,old,new):
    ax = fig.gca()
    assert old.nsensors == new.nsensors

    def drawSensor(label,color,pos,pointingat):
       sightline = lin.Line2D( [pos[0],pointingat[0]] , [pos[1],pointingat[1]], color=color, label=label)
       ax.add_artist(sightline)
       sensor=plt.Circle((pos[0],pos[1]),3,color=color,fill=True, alpha=0.5)
       ax.add_artist(sensor)
       if(label): legpointers.append(sensor)
       if(label): leglables.append(label)

    for i in range(0,old.nsensors):
        pos = array(old.sensors[i].pos)
        angle = old.sensors[i].angle
        rold = old.sensors[i].GetCalibratedDist()
        shift = array([old.x,old.y])
        pointingat = pos - array([rold * math.cos(radians(angle)), rold * math.sin(radians(angle))])

        #draw ideal position
        verts = [
            (castor_inner_octant_radius, -25.4/2), # left, bottom
            (castor_inner_octant_radius, 25.4/2), # left, top
            (castor_inner_octant_radius+15.25, 25.4/2), # right, top
            (castor_inner_octant_radius+15.25, -25.4/2), # right, bottom
            (0., 0.), # ignored
            ]

        codes = [Path.MOVETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.CLOSEPOLY,
                 ]
        for j in range(len(verts)):
            verts[j] = rotatePoint([0,0],verts[j],angle)
        path = Path(verts, codes)
        rec = patches.PathPatch(path, facecolor='0.9', lw=2)
        ax.add_patch(rec)
        drawSensor("nominal position" if i==0 else "","0.5",pos,pointingat); #no name skips label for i>1, 0.5 is grey colour
        #draw fitted old position
        drawSensor("fitted w/o B-field" if i==0 else "","r",pos+shift,pointingat+shift)
        posnew = array(new.sensors[i].pos)
        anglenew = new.sensors[i].angle
        rnew = new.sensors[i].GetCalibratedDist()
        shiftnew = array([new.x,new.y])
        
        pointingatnew = posnew - array([rnew * math.cos(radians(angle)), rnew * math.sin(radians(angle))])
        #draw fitted new position
        drawSensor("fitted w/ B-field" if i==0 else "","g",posnew+shiftnew,pointingatnew+shiftnew);
        #arrows and text
        ax.annotate("",
                xy=(posnew[0]+shiftnew[0],posnew[1]+shiftnew[1]), xycoords='data',
                xytext=(pos[0]+shift[0],pos[1]+shift[1]), textcoords='data',
                arrowprops=dict(arrowstyle="fancy", #linestyle="dashed",
                                color="0.1",
                                shrinkA=0, shrinkB=0,
                                patchA=None,
                                patchB=None,
                                connectionstyle="arc3,rad=0.2",
                                ),
                )


        text = ("Far" if new.isFarHalf else "Near") + " (without B field)\nx={0:.2f}+-{1:.2f}\ny={2:.3f}+-{3:.2f}".format(old.x,(old.xeu-old.xel)/2,old.y,(old.yeu-old.yel)/2)
        an1 = ax.annotate(text, xy=(0.02 if new.isFarHalf else 0.59,0.97), xycoords="axes fraction",
                  va="top", ha="left" if new.isFarHalf else "right",
                  bbox=dict(boxstyle="round", fc="w"))

        from matplotlib.text import OffsetFrom
        offset_from = OffsetFrom(an1, (0.5, 0))
        text = ("Far" if new.isFarHalf else "Near") + " (with B field)\nx={0:.2f}+-{1:.2f}\ny={2:.3f}+-{3:.2f}".format(new.x,(new.xeu-new.xel)/2,new.y,(new.yeu-new.yel)/2)
        an2 = ax.annotate(text, xy=(0.1, 0.1), xycoords="data",
                          xytext=(0, -10), textcoords=offset_from,
                          # xytext is offset points from "xy=(0.5, 0), xycoords=at"
                          va="top", ha="center",
                          bbox=dict(boxstyle="round", fc="w"))

plt.xlabel('x [mm]')
plt.ylabel('y [mm]')
plt.title('IR Sensors (IP side) [Jan 08 -> Jan 17]')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')

draw(fig,farside_old,farside_new)
leg = ax.legend(legpointers,leglables,loc='upper right', fancybox=True)
draw(fig,nearside_old,nearside_new)

plt.axis([-70, 70, -60, 80]) #make it square

plt.savefig("ir_IP_pos.png",bbox_inches="tight")
plt.savefig("ir_IP_pos.pdf",bbox_inches="tight")

