from __future__ import division
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt

data=np.array([
      [663,80],
      [672,151],
      [678,217],
      [684,249],
      [696,294],
      [716,337],
      [741,364],
      [778,372],
      [805,368],
      [827,357],
      [844,335],
      [858,309],
      [873,270],
      [880,235],
      [886,182],
      [899,88]])
data-=[778,372]
data=data/4
data2=np.array(data)
data.T[0]*=-1
data=np.array(data.tolist()+data2.tolist()) #poor man's version or mirror data to get symmetric fit

x=data.T[0]
y=-data.T[1]


r=27.5

fitfunc = lambda p, x: p[0]*x**2 + p[1]*x**4 + p[2]*x**6 # Target function (split in 2 functions to draw)
errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
p0 = [1,1,1]
p1, success = optimize.leastsq(errfunc, p0[:], args=(x, y))
if success: print "Fitted successfully"

xs=np.linspace(x.min(),x.max())
xc=np.linspace(-r,r,200)
yc = r-r*np.sin(np.arccos(xc/r))

fig=plt.figure(figsize=[8,8])
plt.plot(x,y,"kx")
plt.plot(xs,fitfunc(p1, xs),color="r")
plt.plot(xc,yc,ls="dashed")
plt.savefig("special_calibration_rough_1.png")

fig=plt.figure(figsize=[8,4])
plt.plot(xc,yc-fitfunc(p1,xc))
plt.savefig("special_calibration_rough_2.png")
