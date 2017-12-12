# -*- coding: utf-8 -*-
import numpy as np
#datadir = './TannerData/'
import matplotlib.pyplot as p
from scipy.interpolate import interp1d
import scipy.optimize as opt
import material_est as mm
p.ion()

reload(mm)

V = {'emod':1000.,
       'sa':0.,
       'C':0.,
       'C1':100.,
       'C20':0.1,
       'C3':0.,
       'C4':0.,
       'C5':0.,
       'C6':0.,
       'C7':0.,
       'C8':0.,
       'C9':0.,
       'm':2.,
       'n':3.,
       'q':0.,
       'xi':1.,
       'rate0':1.,
       'sig0':100.,
       'integrate':'implicit',
       'tempformat':'Kelvin',
       'timeformat':'Total'}
       
       
rate = 1e-3
temp = 1123.
times = [0.]
ttime = 0.
strains0 = np.linspace(0,1.,51)
straincyc = np.r_[strains0,strains0[-2::-1],-strains0[1:],-strains0[-2::-1]]
strains = np.r_[straincyc,straincyc[1:],straincyc[1:],straincyc[1:],strains0[1::]]
strainsff = np.r_[straincyc,strains0[1::]]
displ = np.exp(strains)-1
for i in range(strains.size-1):
  ttime+=np.abs(strains[i+1]-strains[i])/rate
  times+=[ttime]

stress,strain = mm.stresscalc(displ,times[1],temp,**V)

p.close(1)
p.figure(1)
p.plot(strain,stress,'k-')
p.title('Dislocation Based')
