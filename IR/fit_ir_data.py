import numpy
import math
import os
import matplotlib as ml
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
    def __init__(self, isFarHalf,sensor_angles,sensor_pos,rnew):
        assert len(sensor_angles) == len(rnew), '#angles must be #rnew'
        assert len(sensor_angles) == len(sensor_pos), '#angles must be #positions'
        self.far = isFarHalf
        self.sensor_angles = sensor_angles
        self.sensor_pos = sensor_pos #in ideal geometry
        self.rnew = rnew
        self.verbosity = 1

    def setVerbosity(self, verbosity):
        #0: no output, 1:essential results, 2:everything
        self.verbosity = verbosity

    def fit_pos(self):
        n = len(self.sensor_angles)
        def f(x,y):
            chi2 = 0
            sigma = 1
            for i in range(0,n):
                pos = array(self.sensor_pos[i])
                angle = self.sensor_angles[i]
                r = self.rnew[i]
                pointingat = pos - array([r * math.sin(radians(angle)), r * math.cos(radians(angle))])
                chi2 += distance_to_beampipe(pointingat[0]+x,pointingat[1]+y) ** 2 / sigma ** 2
            return chi2
        m = minuit(f)
        if verbosity > 1 : m.printMode = 1
        m.migrad()
        x = m.values["x"]
        y = m.values["y"]
        chi2 = m.fval / n
        if verbosity > 0 :print("far half: " if self.far else "near half: ","x=",x," y=",y, "chi2: ", chi2)
        return

verbosity = 1 
far_angles = [22.5,-67.5]
far_pos = []
for angle in far_angles:
    far_pos.append([math.cos(math.radians(angle))*castor_inner_octant_radius , math.sin(math.radians(angle))*castor_inner_octant_radius])
print("Position of far side sensors is: ", far_pos)

far_r_old = [10,15]
farside_old = castor_half(1,far_angles,far_pos,far_r_old)
farside_old.setVerbosity(verbosity)
farside_old.fit_pos()

far_r_new = [10,20]
farside_new = castor_half(1,far_angles,far_pos,far_r_new)
farside_new.setVerbosity(verbosity)
farside_new.fit_pos()

