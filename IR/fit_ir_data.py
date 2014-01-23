import numpy
import math
import os
import matplotlib as ml
from minuit2 import Minuit2 as minuit

verbosity = 2 #0: no output 1:essential results 2:everything

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
    def __init__(self, isFarHalf,sensor_angles,dr):
        assert len(sensor_angles) == len(dr), '#angle must be #dr'
        self.far = isFarHalf
        self.sensor_angles = sensor_angles
        self.dr = dr
    def fit_pos(self):
        n = len(self.sensor_angles)
        def f(x,y):
            chi2 = 0
            sigma = 1
            for i in range(0,n):
                chi2 += distance_to_beampipe(x,y) ** 2 / sigma ** 2
            return chi2
        m = minuit(f)
        global verbosity
        if verbosity > 1 : m.printMode = 1
        m.migrad()
        print(m.values["x"], m.values["y"], m.fval)

farside = castor_half(1,[-22.5,67.5],[10,13])
farside.fit_pos()

nearside = castor_half(0,[-22.5,67.5],[10,13])
#nearside.fit_pos()


