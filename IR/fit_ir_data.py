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
from scipy import arange, array, exp, optimize
from ROOT import *

create_chi2_plot = True
castor_inner_octant_radius = 40.6507 #cos(22.5) * 44
beampipe_r = (57+0.5)/2 #in mm #+white paper = 0.5mm
beampipe_x_shift = 0#-0.7 #hauke's value
beampipe_y_shift = 0#-1.4
activeSensingAreaShift = -6.35 #mm distance because the sending device is not in the middle of the housing - only for IP side!

""" Calibration data from ROOT - one for each sensor"""
"""
GeometricCorrectionFile = TFile.Open("GeometricCorrection.root")
Correction_IP_NEAR_TOP = GeometricCorrectionFile.Get("IP_NEAR_TOP")
Correction_IP_NEAR_BOTTOM = GeometricCorrectionFile.Get("IP_NEAR_BOTTOM")
Correction_IP_FAR_TOP = GeometricCorrectionFile.Get("IP_FAR_TOP")
Correction_IP_FAR_BOTTOM = GeometricCorrectionFile.Get("IP_FAR_BOTTOM")

Correction_NONIP_NEAR_TOP = GeometricCorrectionFile.Get("NONIP_NEAR_TOP")
Correction_NONIP_NEAR_CENTER = GeometricCorrectionFile.Get("NONIP_NEAR_CENTER")
Correction_NONIP_NEAR_BOTTOM = GeometricCorrectionFile.Get("NONIP_NEAR_BOTTOM")
Correction_NONIP_FAR_TOP = GeometricCorrectionFile.Get("NONIP_FAR_TOP")
Correction_NONIP_FAR_CENTER = GeometricCorrectionFile.Get("NONIP_FAR_CENTER")
Correction_NONIP_FAR_BOTTOM = GeometricCorrectionFile.Get("NONIP_FAR_BOTTOM")
"""
""" ++++++++++++++++++++++++++++++++++++++++++++++ """

""" ++++++++++++++++ """
""" Helper Functions """
def rotatePoint(centerPoint,point,angle):
    """Rotates a point around another centerPoint. Angle is in degrees.
    Rotation is counter-clockwise"""
    angle = math.radians(angle)
    temp_point = point[0]-centerPoint[0] , point[1]-centerPoint[1]
    temp_point = ( temp_point[0]*math.cos(angle)-temp_point[1]*math.sin(angle) , temp_point[0]*math.sin(angle)+temp_point[1]*math.cos(angle))
    temp_point = temp_point[0]+centerPoint[0] , temp_point[1]+centerPoint[1]
    return temp_point

def distanceTwoPoints(x1,y1,x2,y2):
    return sqrt( (x1-x2)**2 + (y1-y2)**2 )
""" ++++++++++++++++++++++++++++++++++++++++++++++ """

class castor_half(object):
    """Class for a castor half that takes sensor information and can shift/fit position so that sensors point to beam pipe"""
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

    def SetVerbosity(self, verbosity):
        #0: no output, 1:essential results, 2:everything
        self.verbosity = verbosity

    def GetChi2(self,x,y):
        """chi2 function. loops over sensors and adds up chi2. Return chi2,ndf"""
        chi2 = 0
        for iSen in self.sensors:
            chi2 += iSen.GetChi2(x,y)
        ndf = len(self.sensors);
        return chi2,ndf

class castor(object):
    """Class object to store both halves and sensors (e.g. potentiometers) that are connected to both halves"""
    def __init__(self, name, halves):
        self.nsensors = 0
        self.sensors = []
        self.verbosity = 1
        self.half = {}
        for iHalf in halves:
            self.half["far" if iHalf.isFarHalf else "near"] = iHalf
        
        self.checkConsistency()

    def SetVerbosity(self,verb):
        self.verbosity = verb
        for key,iHalf in self.half.iteritems():
            iHalf.SetVerbosity(self.verbosity)

    def checkConsistency(self):
        for iSen in self.sensors:
            assert type(iSen) == openingSensor, "only opening sensors should be added to castor instance"
        assert len(self.half)==2,"Two halfes needed in castor instance"

    def addSensors(self, sensors):
        """add list of sensors"""
        self.nsensors = len(sensors)
        self.sensors = sensors
        self.checkConsistency()

    def GetChi2(self,x_far,y_far,x_near,y_near):
        near_half = self.half["near"]
        far_half = self.half["far"]
        chi2 = 0.
        ndf = 0

        c,n = near_half.GetChi2(x_near,y_near)
        chi2 += c
        ndf += n
        c,n = far_half.GetChi2(x_far,y_far)
        chi2 += c
        ndf += n

        for iSen in self.sensors:
            c = iSen.GetChi2(x_near,y_near,x_far)
            chi2 += c
            ndf += 1
#            assert ndf > 4, "at least one sensor needed for fit? huh?"
#        if y_near >0 : chi2 += 30*abs(y_near) #remove later. constrain to negative y
        return chi2

    def fit(self):
        self.checkConsistency()

        near_half = self.half["near"]
        far_half = self.half["far"]

        def f(x_far,y_far,x_near,y_near):
            chi2 = 0.
            ndf = 0

            c,n = near_half.GetChi2(x_near,y_near)
            chi2 += c
            ndf += n
            c,n = far_half.GetChi2(x_far,y_far)
            chi2 += c
            ndf += n

            for iSen in self.sensors:
                c = iSen.GetChi2(x_near,y_near,x_far)
                chi2 += c
                ndf += 1
#            assert ndf > 4, "at least one sensor needed for fit? huh?"
#            if y_near >0 : chi2 += 30*abs(y_near) #remove later. constrain to negative y
            return chi2

        m = minuit(f)
        m.values["x_far"] = -5 #starting values
        m.values["y_far"] = -10
        m.values["x_near"] = 15
        m.values["y_near"] = -10
#        m.tol = 10000 #how close result must be to estimated minumum


        if verbosity > 1 : print "Starting values: ", m.values
        if verbosity > 1 : m.printMode = 1
        m.migrad()
        m.minos()

        near_half.x = m.values["x_near"] #maybe one should add
        near_half.y = m.values["y_near"]
        near_half.xeu = m.merrors["x_near", 1]
        near_half.yeu = m.merrors["y_near", 1]
        near_half.xel = m.merrors["x_near", -1]
        near_half.yel = m.merrors["y_near", -1]

        far_half.x = m.values["x_far"]
        far_half.y = m.values["y_far"]
        far_half.xeu = m.merrors["x_far", 1]
        far_half.yeu = m.merrors["y_far", 1]
        far_half.xel = m.merrors["x_far", -1]
        far_half.yel = m.merrors["y_far", -1]

        chi2 = m.fval
        if verbosity > 0 :
            print "Fitted position (chi2/ndf={chi2:.2f})".format(chi2=chi2)
            print "Near half: (x,y)=({0:.2f}{1:+.2f}{2:+.2f},{3:.2f}{4:+.2f}{5:+.2f})".format(self.half["near"].x,self.half["near"].xeu,self.half["near"].xel,self.half["near"].y,self.half["near"].yeu,self.half["near"].yel)
            print "Far half: (x,y)=({0:.2f}{1:+.2f}{2:+.2f},{3:.2f}{4:+.2f}{5:+.2f})".format(self.half["far"].x,self.half["far"].xeu,self.half["far"].xel,self.half["far"].y,self.half["far"].yeu,self.half["far"].yel)
            print "Opening: {op:.2f} mm".format(op = near_half.x - far_half.x)
        return


""" ++++++++++++++++++++++++++ """
""" SENSOR classes """
""" ++++++++++++++++++++++++++ """


class sensor(object):
    """ABSTRACT sensor class. stores information about position and angle about each sensor"""
    def __init__(self, pos, angle, name):
        """angle counter clock-wise in deg [-180,180] starting from 3pm looking from IP"""
        self.pos = pos
        self.angle = angle
        self.meas_r = 0
        self.meas_r_err = 0
        self.cal_meas = []
        self.cal_true = []
        self.cal_spline = None
        self.verbosity = 1
        self.name = name

        assert len(self.pos) == 2, 'pos must be provided as list [x,y]'
        assert -180.5 < self.angle < 180.5, 'error with angle for sensor'

        if self.verbosity > 2:
            print "             nominal position (x,y)=({0[0]:.2f},{0[1]:.2f})".format(self.pos)

    @classmethod
    def fromsensor(cls, sensor):
        """copy constructor"""
        new = cls(sensor.pos,sensor.angle,sensor.name)
        new.verbosity = sensor.verbosity
        new.SetCalibrationData(sensor.cal_meas,sensor.cal_true)
        new.SetDist(sensor.meas_r,sensor.meas_r_err)
        return new

    def SetCalibrationData(self,meas,true):
        """set measured calibration data. like [-1.2, 10.1, 20.4] and [0, 10, 20]"""
        if not meas or not true: return

        self.cal_meas = meas
        self.cal_true = true
        assert len(self.cal_meas) == len(self.cal_true), 'error in calibration data'
        self.cal_spline = UnivariateSpline(self.cal_meas,self.cal_true,k=2)

    def DrawCalibration(self):
        """Draws plot for calibration data by SetData()"""
        assert self.cal_spline, "!cannot draw calibration curve since calibration data not set"
        fig = plt.figure(figsize=[8,8])
        ax = fig.gca()
        plt.plot(self.cal_meas, self.cal_true, 'bs')
        xcal = linspace(-5, 40, 1000)
        ycal = self.cal_spline(xcal)
        plt.plot([-5,40], [-5,40],"-k")
        plt.plot(xcal, ycal,"-b",lw=2)
        plt.xlabel('measured [mm]')
        plt.ylabel('truth [mm]')
        plt.title('Sensors (IP side) {name}'.format(name=self.name.replace("_"," ")))
        plt.savefig("ir_sens_calib_{name}.png".format(name=self.name.replace(" ","_")))

    def GetCalibratedDist(self):
        """apply calibration and return distance. if outside of calibration data it will be extrapolated"""
        if (self.cal_spline):
            return self.cal_spline(self.meas_r)[0] #for other scipy version array is returned. [0] is needed
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

    def _GetPointingAt(self,x,y):
        assert False, "please define pure virtual method"

    def _AwayFromTarget(self):
        assert False, "please define pure virtual method"

    def GetChi2(self,x,y):
        """returns distance to beam pipe circle for all vector(sensorpos)-vector(r) for given (x,y) shift"""
        delta,sigma = self._AwayFromTarget(self._GetPointingAt(x,y))
        return delta ** 2 / sigma ** 2

class infraredBeamPipeSensor(sensor):
    """subclass of sensor that measures the distance to mount or opening"""
    def __init__(self,pos,angle,name,orientationSwap=False):
        sensor.__init__(self,pos,angle,name)
        sensor.__seenBeampipe = None
        self.swapped = orientationSwap

    def SetGeometryCorrection(self,CalGraph):
        """ Assign the right Correction Graph to the Sensor """
        self.CorrectionGraph = CalGraph

    def GetCorrectedDist(self,x,y):
        """apply calibration and return distance. using the 3D calibration for large angles, x and y are the current sensor position"""
	if self.CorrectionGraph:
	    real_r = self.CorrectionGraph.Eval(self.meas_r,rotatePoint([0,0],array([x,y]), -self.angle)[1])
	else:
	    real_r = self.meas_r
        if real_r == 0 :
            if self.verbosity > 1:
                print "No calibration data set. Using uncalibrated"
            return self.meas_r #if 3D calibration fails return uncalibrated distance
        else:
            return real_r

    def distance_to_beampipe(self, x , y , xe=0 , ye=0 ):
        """calculate distance for any given point to beam pipe outer circle"""
        global beampipe_x_shift
        global beampipe_y_shift
        global beampipe_r
        projP = rotatePoint([0,0],[x,y], -self.angle) #now incident angle of line of sight is parallel to x-axis coming from near side
        x = projP[0]
        y = projP[1]
        if abs(y)<beampipe_r:
            beampipe_intersect_x = math.sqrt((beampipe_r**2)-(y**2))
        else:
            beampipe_intersect_x = 0
        distance = abs(beampipe_intersect_x - x)
        error = sqrt(xe**2+ye**2);
        return distance, error

    def _AwayFromTarget(self,pointing_at):
        r = self.meas_r
        error_r = sqrt(2**2 + self.GetDistError()**2) #2 mm sys + stat error
        error_theta = radians(1.5) #deg systematic uncertainty

        dpxdr     = -math.cos(radians(self.angle))
        dpxdtheta = r * math.sin(radians(self.angle))
        dpydr     = -math.sin(radians(self.angle))
        dpydtheta = -r * math.cos(radians(self.angle))

        xe = sqrt( error_r**2 * dpxdr**2 + error_theta**2 * dpxdtheta**2)#x and y is initial sensor positions from drawings. maybe no uncertainty
        ye = sqrt( error_r**2 * dpydr**2 + error_theta**2 * dpydtheta**2)#x and y is initial sensor positions from drawings. maybe no uncertainty

        if self.verbosity>1: print "xe=",xe,"ye=",ye, " --> delta=",delta,"sigma=",sigma
        delta,sigma = self.distance_to_beampipe(pointing_at[0],pointing_at[1], xe, ye)
        return delta,sigma

    def _GetPointingAt(self,x,y):
        """calculates r in direction of angle. returns pointing absolute position, x and y are the CASTOR shifts"""
        r = self.GetCorrectedDist(self.pos[0]+x, self.pos[1]+y)
        return array(self.pos) - array([r * math.cos(radians(self.angle)), r * math.sin(radians(self.angle))]) + array([x,y])


    def drawNominalPosition(self,label,ax,leglabels,legpointers):
       
        #draw ideal position
        projP=rotatePoint([0,0],self.pos,-self.angle)
        x=projP[0]
        y=projP[1]

        if not self.swapped:
            verts = [
                (x,y-6.2), # left, bottom
                (x,y+19.2), # left, top
                (x+15.25,y+19.2), # right, top
                (x+15.25,y-6.2), # right, bottom
                (0., 0.), # ignored
                ]

            """     old definition from Colin, only for IP side possible       
            verts = [
                (castor_inner_octant_radius, -25.4/2), # left, bottom
                (castor_inner_octant_radius, 25.4/2), # left, top
                (castor_inner_octant_radius+15.25, 25.4/2), # right, top
                (castor_inner_octant_radius+15.25, -25.4/2), # right, bottom
                (0., 0.), # ignored
                ]
            """
          
        else:
            verts = [
                (x,y+6.2), # left, bottom
                (x,y-19.2), # left, top
                (x+15.25,y-19.2), # right, top
                (x+15.25,y+6.2), # right, bottom
                (0., 0.), # ignored
                ]

        codes = [Path.MOVETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.CLOSEPOLY,
                 ] #this draws the rectanle
        for j in range(len(verts)):
            verts[j] = rotatePoint([0,0],verts[j],self.angle) #and rotated according to angle of sensor
        path = Path(verts, codes)
        rec = patches.PathPatch(path, facecolor='0.9', lw=2)
        ax.add_patch(rec)

        sensor=plt.Circle((self.pos[0],self.pos[1]),3,color="0.5",fill=True, alpha=0.5)
        ax.add_artist(sensor)
        if(label): legpointers.append(sensor)
        if(label): leglables.append(label)


    def drawSensor(self,label,color,x,y,ax,leglabels,legpointers):
        pointing_at = array(self.pos) + array([x,y]) - array([self.GetCorrectedDist(self.pos[0]+x, self.pos[1]+y) * math.cos(radians(self.angle)), self.GetCorrectedDist(self.pos[0]+x, self.pos[1]+y) * math.sin(radians(self.angle))])
        sightline = lin.Line2D( [self.pos[0]+x,pointing_at[0]] , [self.pos[1]+y,pointing_at[1]], color=color, label=label)
        if color != "0.5": ax.add_artist(sightline)
        sensor=plt.Circle((self.pos[0]+x,self.pos[1]+y),3,color=color,fill=True, alpha=0.5)
        ax.add_artist(sensor)
        if(label): legpointers.append(sensor)
        if(label): leglables.append(label)


class openingSensor(sensor):
    """subclass of sensor that measures the opening distance. """
    def __init__(self,pos_y,name):
        """pos = at (0,y)"""
        sensor.__init__(self,(0,pos_y),angle=0,name=name)

    @classmethod
    def fromsensor(cls, sensor):
        """copy constructor"""
        new = cls(sensor.pos[0],sensor.name)
        new.verbosity = sensor.verbosity
        new.SetCalibrationData(sensor.cal_meas,sensor.cal_true)
        new.SetDist(sensor.meas_r,sensor.meas_r_err)
        return new

    def GetChi2(self,x_near,y_near,x_far):
        """overwriting GetChi2 from base class to accept x from other half"""
        delta,sigma = self._AwayFromTarget(self._GetPointingAt(x_near,y_near),x_far)
        return delta ** 2 / sigma ** 2

    def _AwayFromTarget(self,pointing_at,x_other_half):
        delta = abs(pointing_at[0]-x_other_half)
        sigma=self.GetDistError()
        return delta,sigma

    def _GetPointingAt(self,x,y):
        """overwriting GetChi2 from base class to avoid geometric calibration correction"""
        r = self.meas_r
        return array(self.pos) - array([r * math.cos(radians(self.angle)), r * math.sin(radians(self.angle))]) + array([x,y])


class fixationPotentiometer(sensor):
    """subclass of potentiometer sensors that measure the distance to the beam pipe fixation at the non IP side"""


########FITTING########
verbosity = 1
fitNonIP = True

### First: IP side ######
######FITTING OLD CASTOR####
####Defining sensors####
sensor_fartop = infraredBeamPipeSensor(rotatePoint([0,0],[castor_inner_octant_radius,activeSensingAreaShift], 180-67.5), 180-67.5, "far top") #activeSensingAreaShift applied to y position unrotated. this is correct (clockwise)
sensor_farbot = infraredBeamPipeSensor(rotatePoint([0,0],[castor_inner_octant_radius,activeSensingAreaShift], -180+22.5), -180+22.5, "far bot")
sensor_fartop.SetDist(8.68439,0)
sensor_farbot.SetDist(20.2883,0.18)
sensor_farbot.SetGeometryCorrection(None)
sensor_fartop.SetGeometryCorrection(None)
#sensor_fartop.SetCalibrationData([0.5,10.1,20.2], [0,10,20])
#sensor_farbot.SetCalibrationData([-2.,9.8,19.1], [0,10,20])
#sensor_fartop.DrawCalibration()
#sensor_farbot.DrawCalibration()

sensor_neartop = infraredBeamPipeSensor(rotatePoint([0,0],[castor_inner_octant_radius,activeSensingAreaShift], 67.5), 67.5, "near top") #activeSensingAreaShift applied to y position unrotated. this is correct (clockwise)
sensor_nearbot = infraredBeamPipeSensor(rotatePoint([0,0],[castor_inner_octant_radius,activeSensingAreaShift], -22.5), -22.5, "near bot")
#sensor_neartop.SetCalibrationData([0.8,10.7,19.9], [0,10,20])
#sensor_nearbot.SetCalibrationData([0.7,11.1,20.5], [0,10,20])
#sensor_neartop.DrawCalibration()
#sensor_nearbot.DrawCalibration()
sensor_neartop.SetDist(15.2248,0.04)
sensor_nearbot.SetDist(25.8462,0.82)
sensor_nearbot.SetGeometryCorrection(None)
sensor_neartop.SetGeometryCorrection(None)

#potis are fixed to far side. castor flange 200mm minus ~10mm
sensor_pot_top = openingSensor(190,"potentiometer top")
sensor_pot_top.SetDist(8.3,2.) #2 mm uncertainty
sensor_pot_bottom = openingSensor(-190,"potentiometer bottom")
sensor_pot_bottom.SetDist(14.1,2.)

####Defining Halves#####
nearside_old = castor_half("nearside_old",[sensor_neartop,sensor_nearbot])
farside_old = castor_half("farside_old",[sensor_fartop,sensor_farbot])

####Defining Whole Castor Instance####
castor_old = castor("castor at old position",[nearside_old,farside_old])
castor_old.addSensors([sensor_pot_bottom,sensor_pot_top])
castor_old.SetVerbosity(verbosity)

####Fitting####
print "Before B field"
castor_old.fit()

######FITTING NEW CASTOR#######
####Defining Sensors####
sensor_fartop = infraredBeamPipeSensor.fromsensor(sensor_fartop)
sensor_farbot = infraredBeamPipeSensor.fromsensor(sensor_farbot)
sensor_fartop.SetDist(12.6556,0)
sensor_farbot.SetDist(22.2788,0.04)
sensor_farbot.SetGeometryCorrection(None)
sensor_fartop.SetGeometryCorrection(None)

sensor_neartop = infraredBeamPipeSensor.fromsensor(sensor_neartop)
sensor_nearbot = infraredBeamPipeSensor.fromsensor(sensor_nearbot)
sensor_neartop.SetDist(32.9261,4e-6)
sensor_nearbot.SetDist(32.6869,0.18)
sensor_nearbot.SetGeometryCorrection(None)
sensor_neartop.SetGeometryCorrection(None)

sensor_pot_top = openingSensor.fromsensor(sensor_pot_top)
sensor_pot_top.SetDist(19.1,2.) #2 mm uncertainty
sensor_pot_bottom = openingSensor.fromsensor(sensor_pot_bottom)
sensor_pot_bottom.SetDist(22.1,2.)

####Defining Halves####
nearside_new = castor_half("nearside_new",[sensor_neartop,sensor_nearbot])
farside_new = castor_half("farside_new",[sensor_fartop,sensor_farbot])

####Defining Whole Castor Instance####
castor_new = castor("castor at new position",[nearside_new,farside_new])
castor_new.addSensors([sensor_pot_bottom,sensor_pot_top])
castor_new.SetVerbosity(verbosity)

####Fitting####
print "After B field"
castor_new.fit()

def printShift(old,new):
    dx = new.x-old.x
    dy = new.y-old.y
    dr = math.sqrt(dx**2+dy**2)
    print "Shift of the IP side due to magnetic field ({half}): dx={dx:.2f} dy={dy:.2f} dr={dr:.2f}".format(half="far  half" if old.isFarHalf else "near half",dx=dx,dy=dy,dr=dr)

print("\n\n")
printShift(castor_old.half['near'],castor_new.half['near'])
printShift(castor_old.half['far'],castor_new.half['far'])



#### NON-IP side #####
if fitNonIP :

    ######FITTING OLD CASTOR####
    ####Defining sensors####
    sensor_nonip_fartop = infraredBeamPipeSensor(rotatePoint([0,0],[castor_inner_octant_radius,1.238], 180-67.5), 180-67.5, "far top") #activeSensingAreaShift applied to y position unrotated. this is correct (clockwise)
    sensor_nonip_farcenter = infraredBeamPipeSensor(rotatePoint([0,0],[castor_inner_octant_radius,-1.6], 180), 180, "far center", True)
    sensor_nonip_farbot = infraredBeamPipeSensor(rotatePoint([0,0],[castor_inner_octant_radius,2.238], -180+67.5), -180+67.5, "far bot", True)
    sensor_nonip_fartop.SetDist(6.27352,0.08)
    sensor_nonip_farcenter.SetDist(13.5765,0.13)
    sensor_nonip_farbot.SetDist(18.7727,0.06)
    sensor_nonip_fartop.SetGeometryCorrection(None)
    sensor_nonip_farcenter.SetGeometryCorrection(None)
    sensor_nonip_farbot.SetGeometryCorrection(None)
    #sensor_nonip_fartop.SetCalibrationData([0.5,10.1,20.2], [0,10,20])
    #sensor_nonip_farbot.SetCalibrationData([-2.,9.8,19.1], [0,10,20])
    #sensor_nonip_fartop.DrawCalibration()
    #sensor_nonip_farbot.DrawCalibration()

    sensor_nonip_neartop = infraredBeamPipeSensor(rotatePoint([0,0],[castor_inner_octant_radius,2.238], 67.5), 67.5, "near top", True) #activeSensingAreaShift applied to y position unrotated. this is correct (clockwise)
    sensor_nonip_nearcenter = infraredBeamPipeSensor(rotatePoint([0,0],[castor_inner_octant_radius,-1.6], 0), 0, "near center", True) #activeSensingAreaShift applied to y position unrotated. this is correct (clockwise)
    sensor_nonip_nearbot = infraredBeamPipeSensor(rotatePoint([0,0],[castor_inner_octant_radius,-1.2], -67.5), -67.5, "near bot")
    #sensor_nonip_neartop.SetCalibrationData([0.8,10.7,19.9], [0,10,20])
    #sensor_nonip_nearbot.SetCalibrationData([0.7,11.1,20.5], [0,10,20])
    #sensor_nonip_neartop.DrawCalibration()
    #sensor_nonip_nearbot.DrawCalibration()
    sensor_nonip_neartop.SetDist(13.3061,0.03)
    sensor_nonip_nearcenter.SetDist(26.8516,0)
    sensor_nonip_nearbot.SetDist(21.9653,0.06)
    sensor_nonip_neartop.SetGeometryCorrection(None)
    sensor_nonip_nearcenter.SetGeometryCorrection(None)
    sensor_nonip_nearbot.SetGeometryCorrection(None)
    """
    #potis are fixed to far side. castor flange 200mm minus ~10mm
    sensor_nonip_pot_top = openingSensor(190,"potentiometer top")
    sensor_nonip_pot_top.SetDist(8.3,2.) #2 mm uncertainty
    sensor_nonip_pot_bottom = openingSensor(-190,"potentiometer bottom")
    sensor_nonip_pot_bottom.SetDist(14.1,2.)
    """

    ####Defining Halves#####
    nearside_nonip_old = castor_half("nearside_nonip_old",[sensor_nonip_neartop,sensor_nonip_nearcenter,sensor_nonip_nearbot])
    farside_nonip_old = castor_half("farside_nonip_old",[sensor_nonip_fartop,sensor_nonip_farcenter,sensor_nonip_farbot])

    ####Defining Whole Castor Instance####
    castor_nonip_old = castor("castor at old position",[nearside_nonip_old,farside_nonip_old])
    #castor_nonip_old.addSensors([sensor_pot_bottom,sensor_pot_top])
    castor_nonip_old.SetVerbosity(verbosity)

    ####Fitting####
    print "Before B field"
    castor_nonip_old.fit()

    ######FITTING NEW CASTOR#######
    ####Defining Sensors####
    sensor_nonip_fartop = infraredBeamPipeSensor.fromsensor(sensor_nonip_fartop)
    sensor_nonip_farcenter = infraredBeamPipeSensor.fromsensor(sensor_nonip_farcenter)
    sensor_nonip_farbot = infraredBeamPipeSensor.fromsensor(sensor_nonip_farbot)
    sensor_nonip_fartop.SetDist(10.7459,0)
    sensor_nonip_farcenter.SetDist(13.5471,0.01)
    sensor_nonip_farbot.SetDist(14.717,0.25)
    sensor_nonip_fartop.SetGeometryCorrection(None)
    sensor_nonip_farcenter.SetGeometryCorrection(None)
    sensor_nonip_farbot.SetGeometryCorrection(None)

    sensor_nonip_neartop = infraredBeamPipeSensor.fromsensor(sensor_nonip_neartop)
    sensor_nonip_nearcenter = infraredBeamPipeSensor.fromsensor(sensor_nonip_nearcenter)    
    sensor_nonip_nearbot = infraredBeamPipeSensor.fromsensor(sensor_nonip_nearbot)
    sensor_nonip_neartop.SetDist(15.9262,0.08)
    sensor_nonip_nearcenter.SetDist(23.8852,0.06)
    sensor_nonip_nearbot.SetDist(18.8557,0.06)
    sensor_nonip_neartop.SetGeometryCorrection(None)
    sensor_nonip_nearcenter.SetGeometryCorrection(None)
    sensor_nonip_nearbot.SetGeometryCorrection(None)
    """
    sensor_nonip_pot_top = openingSensor.fromsensor(sensor_pot_top)
    sensor_nonip_pot_top.SetDist(19.1,2.) #2 mm uncertainty
    sensor_nonip_pot_bottom = openingSensor.fromsensor(sensor_pot_bottom)
    sensor_nonip_pot_bottom.SetDist(22.1,2.)
    """
    ####Defining Halves####
    nearside_nonip_new = castor_half("nearside_nonip_new",[sensor_nonip_neartop,sensor_nonip_nearcenter,sensor_nonip_nearbot])
    farside_nonip_new = castor_half("farside_nonip_new",[sensor_nonip_fartop,sensor_nonip_farcenter,sensor_nonip_farbot])

    ####Defining Whole Castor Instance####
    castor_nonip_new = castor("castor at new position",[nearside_nonip_new,farside_nonip_new])
    #castor_nonip_new.addSensors([sensor_nonip_pot_bottom,sensor_nonip_pot_top])
    castor_nonip_new.SetVerbosity(verbosity)

    ####Fitting####
    print "After B field"
    castor_nonip_new.fit()

    def printShift(old,new):
        dx = new.x-old.x
        dy = new.y-old.y
        dr = math.sqrt(dx**2+dy**2)
        print "Shift of the non-ip side due to magnetic field ({half}): dx={dx:.2f} dy={dy:.2f} dr={dr:.2f}".format(half="far  half" if old.isFarHalf else "near half",dx=dx,dy=dy,dr=dr)

    print("\n\n")
    printShift(castor_nonip_old.half['near'],castor_nonip_new.half['near'])
    printShift(castor_nonip_old.half['far'],castor_nonip_new.half['far'])



####Creating Chi2 Plot####
if create_chi2_plot:
    fig = plt.figure(figsize=[8,8])
    ax = fig.gca()
    xg ,yg = np.mgrid [-25:45:50j , -40:40:50j] #range and steps
    f=[]
    print "Creating Chi2 plot"
    for i in range(len(xg)):
        g=[]
        for j in range(len(xg[i])):
            g.append(castor_new.GetChi2(-9.2,-4.45,xg[i][j],yg[i][j]))#-6.5,-6.5,xg[i][j],yg[i][j]))
        f.append(g)
    #print xg, "\n yg ", yg, "\n f" , f
    from matplotlib.colors import LogNorm
    p = plt.pcolor(xg,yg,np.array(f),norm=LogNorm())
    cb = fig.colorbar(p, ax=ax)
    plt.savefig("chi2_fixed_far_half.png")


#####DRAWING###### (most of the function should acutally be member functions of the objects itself. please fix)
fig_IP = plt.figure(figsize=[8,8])
ax_IP = fig_IP.gca()
circle1=plt.Circle((beampipe_x_shift,beampipe_y_shift),beampipe_r,color='black',fill=False,label="beampipe",lw=2)
ax_IP.add_artist(circle1)

leglables=["beampipe"]
legpointers=[circle1]

def draw(fig,old,new,leglabels,legpointers):
    ax = fig.gca()
    assert old.nsensors == new.nsensors

    for i in range(0,old.nsensors):
        old.sensors[i].drawNominalPosition("nominal position" if i==0 else "",ax,leglabels,legpointers); #no name skips label for i>1
        #draw fitted old position
        old.sensors[i].drawSensor("fitted w/o B-field" if i==0 else "","r",old.x,old.y,ax,leglabels,legpointers)
        #draw fitted new position
        new.sensors[i].drawSensor("fitted w/ B-field" if i==0 else "","g",new.x,new.y,ax,leglabels,legpointers);
        #arrows and text
        """
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
                ) #arrow for shift due to magnetic field
        """

        text = ("Far" if new.isFarHalf else "Near") + " (without B field)\nx={0:.2f}+-{1:.2f}\ny={2:.3f}+-{3:.2f}".format(old.x,(old.xeu-old.xel)/2,old.y,(old.yeu-old.yel)/2)
        an1 = ax.annotate(text, xy=(0.02 if new.isFarHalf else 0.59,0.97), xycoords="axes fraction",
                  va="top", ha="left" if new.isFarHalf else "right",
                  bbox=dict(boxstyle="round", fc="w")) #box with values 1

        from matplotlib.text import OffsetFrom
        offset_from = OffsetFrom(an1, (0.5, 0))
        text = ("Far" if new.isFarHalf else "Near") + " (with B field)\nx={0:.2f}+-{1:.2f}\ny={2:.3f}+-{3:.2f}".format(new.x,(new.xeu-new.xel)/2,new.y,(new.yeu-new.yel)/2)
        an2 = ax.annotate(text, xy=(0.1, 0.1), xycoords="data",
                          xytext=(0, -10), textcoords=offset_from,
                          # xytext is offset points from "xy=(0.5, 0), xycoords=at"
                          va="top", ha="center",
                          bbox=dict(boxstyle="round", fc="w")) #box with values 2

plt.xlabel('x [mm]')
plt.ylabel('y [mm]')
plt.title('IR Sensors (IP side) [Jan 08 -> Jan 17]')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')

draw(fig_IP,farside_old,farside_new,leglables,legpointers)
leg = ax_IP.legend(legpointers,leglables,loc='upper right', fancybox=True)
draw(fig_IP,nearside_old,nearside_new,leglables,legpointers)

plt.axis([-70, 70, -50, 90]) #make it square

plt.savefig("ir_IP_pos.png",bbox_inches="tight")
plt.savefig("ir_IP_pos.pdf",bbox_inches="tight")

if fitNonIP :
    leglables=["beampipe"]
    legpointers=[circle1]
    fig_nonIP = plt.figure(figsize=[8,8])
    ax_nonIP = fig_nonIP.gca()
    circle2=plt.Circle((beampipe_x_shift,beampipe_y_shift),beampipe_r,color='black',fill=False,label="beampipe",lw=2)
    ax_nonIP.add_artist(circle2)
    plt.xlabel('x [mm]')
    plt.ylabel('y [mm]')
    plt.title('IR Sensors (non IP side) [Jan 08 -> Jan 17]')
    draw(fig_nonIP,farside_nonip_old,farside_nonip_new,leglables,legpointers)
    leg = ax_nonIP.legend(legpointers,leglables,loc='upper right', fancybox=True)
    draw(fig_nonIP,nearside_nonip_old,nearside_nonip_new,leglables,legpointers)

    plt.axis([-70, 70, -50, 90]) #make it square

    plt.savefig("ir_nonIP_pos.png",bbox_inches="tight")
    plt.savefig("ir_nonIP_pos.pdf",bbox_inches="tight")



""" Backup of old code
    def distance_to_beampipe_old(self, x , y , xe=0 , ye=0 ):
        calculate distance for any given point to beam pipe outer circle
        global beampipe_x_shift
        global beampipe_y_shift
        global beampipe_r
        r = math.sqrt((x-beampipe_x_shift)**2 + (y-beampipe_y_shift)**2)
        distance = math.fabs(beampipe_r - r);
        ddisdx = 1 / 2 / r * 2 * (x-beampipe_x_shift) #derivative d(distance)/d(x)
        ddisdy = 1 / 2 / r * 2 * (y-beampipe_y_shift)
        error = math.sqrt( xe**2 * ddisdx**2 + ye**2 * ddisdy**2)
        return distance , error

    def distance_to_beampipe_corr(self, x , y , xe=0 , ye=0 ):
        calculate distance for any given point to beam pipe outer circle but with correction from geometric attenuation of the light
        global beampipe_x_shift
        global beampipe_y_shift
        global beampipe_r
        projP = rotatePoint([0,0],[x,y], -self.angle) #now incident angle of line of sight is parallel to x-axis coming from near side
        x = projP[0]
        y = projP[1]
        beampipe_intersect_x = beampipe_r - self.__getSeenBeampipe(y)
        distance = abs(beampipe_intersect_x - x)
        error = sqrt(xe**2+ye**2);
        return distance, error

    def distance_to_beampipe_proj(self, x , y , xe=0 , ye=0 ):
        calculate distance for any given point to beam pipe outer circle
        global beampipe_x_shift
        global beampipe_y_shift
        global beampipe_r
        
        projP = rotatePoint([0,0],[x,y], -self.angle) #now incident angle of line of sight is parallel to x-axis coming from near side
        #print "projP=",projP
        
        activeSensing = activeSensingAreaShift / 2.

        # if abs(projP[1])+abs(activeSensing) < beampipe_r:
        #     areaOutsideOfCircle = 0
        # else:
        #     areaOutsideOfCircle = (abs(projP[1])+abs(activeSensing)-beampipe_r)/(abs(activeSensing)*2) # how much of the sensing area covers space outside the beampipe circle
        # print "area outside:",areaOutsideOfCircle , math.tanh(areaOutsideOfCircle * math.pi) * beampipe_r
        # if areaOutsideOfCircle > 1: return -1,xe
        
        sensingY = [projP[1]]#np.linspace(abs(projP[1])-abs(activeSensing), min(abs(projP[1])+abs(activeSensing),beampipe_r),num=10)
        sensingDist = []
        for y in sensingY:
            beampipe_intersect_y = y
            if abs(y) < beampipe_r:
                beampipe_intersect_x = sqrt(beampipe_r**2-beampipe_intersect_y**2) #r^2=x^2+y^2 -> solve for x
                sensingDist.append(abs(beampipe_intersect_x - projP[0]))  #line of sight distance
            elif y > beampipe_r:
                sensingDist.append(distanceTwoPoints(projP[0],projP[1], 0,beampipe_r))
            else: #y < beampipe_r
                sensingDist.append(distanceTwoPoints(projP[0],projP[1], 0,-beampipe_r))
            
            
        distance = sum(sensingDist) / len(sensingDist)
        #if abs(projP[1]) > beampipe_r:
        #    distance += (abs(projP[1])-beampipe_r)*(abs(projP[1])-beampipe_r)**2#( math.tanh( (abs(projP[1])/beampipe_r-1) * 5 * math.pi ) +1 ) * beampipe_r #add the beampipe radius for overhanging area = contrain fit
        #print "Distance:", sensingY, sensingDist, distance
        error = sqrt(xe**2+ye**2);
        #ddisdx = 1 / 2 / r * 2 * (x-beampipe_x_shift) #derivative d(distance)/d(x)
        #ddisdy = 1 / 2 / r * 2 * (y-beampipe_y_shift)
        #error = math.sqrt( xe**2 * ddisdx**2 + ye**2 * ddisdy**2)
        return distance , error

"""
