# -*- coding: utf-8 -*-
import numpy as np
#datadir = './TannerData/'
import matplotlib.pyplot as p
from scipy.interpolate import interp1d
import scipy.optimize as opt
import material_mtdb as mm
p.ion()

reload(mm)

V = {'mu0':85000.,
       'd0':8300.,
       't0':200.,
       'nu':0.33,
       'sa':35.,
       's0':50.,
       'Y0':1.,
       'C1':100.,
       'C20':10.,
       'si':0.,
       'a0e':2.,
       'a02':0.55,
       'a0i':0.9,
       'rate0':1e7,
       'rate02':1e10,
       'qe':1.2,
       'pe':0.5,
       'qi':1.5,
       'pi':2./3,
       'tempformat':'Celcius'}
       

strains0 = np.linspace(0,1.,51)
straincyc = np.r_[strains0,strains0[-2::-1],-strains0[1:],-strains0[-2::-1]]
strains = np.r_[straincyc,straincyc[1:],straincyc[1:]]#,straincyc[1:]]
displ = np.exp(strains0)-1

rate = 1.
dtime = strains[1]/rate

stressT0,strainT0 = mm.stresscalc(displ,strains[1]/1,25.,**V)
stressT1,strainT1 = mm.stresscalc(displ,strains[1]/1,50.,**V)
stressT2,strainT2 = mm.stresscalc(displ,strains[1]/1,150.,**V)
stressT3,strainT3 = mm.stresscalc(displ,strains[1]/1,250.,**V)
stressT4,strainT4 = mm.stresscalc(displ,strains[1]/1,350.,**V)

p.close(1)
p.figure(1)
p.plot(strainT1,stressT1,'k:',label=r'50$^o$C')
p.plot(strainT2,stressT2,'k-.',label=r'150$^o$C')
p.plot(strainT3,stressT3,'k--',label=r'250$^o$C')
p.plot(strainT4,stressT4,'k-',label=r'350$^o$C')
p.legend(loc='lower right')
p.title(r'Temperature dependence at $\dot{\varepsilon}=1^{-1}$')


stressR0,strainR0 = mm.stresscalc(displ,strains[1]/0.0001,25.,**V)
stressR1,strainR1 = mm.stresscalc(displ,strains[1]/10000.,25.,**V)
stressR2,strainR2 = mm.stresscalc(displ,strains[1]/10,25.,**V)
stressR3,strainR3 = mm.stresscalc(displ,strains[1]/0.01,25.,**V)
stressR4,strainR4 = mm.stresscalc(displ,strains[1]/0.00001,25.,**V)

p.close(2)
p.figure(2)
p.plot(strainR1,stressR1,'k:',label=r'10000$^{-1}$')
p.plot(strainR2,stressR2,'k-.',label=r'10$^{-1}$')
p.plot(strainR3,stressR3,'k--',label=r'0.01$^{-1}$')
p.plot(strainR4,stressR4,'k-',label=r'0.00001$^{-1}$')
p.legend(loc='lower right')
p.title(r'Rate dependence at 25$^o$C')

#change in temperature:
temps=np.ones(strains0.size)*50.
temps[np.where(strains0<0.3)[0]]=150.
temps[np.where(strains0>0.6)[0]]=250.
stressTC,strainTC = mm.stresscalc(displ,strains[1]/1,temps,**V)

p.close(3)
p.figure(3)
p.plot(strainT1,stressT1,'k:')
p.plot(strainT2,stressT2,'k:')
p.plot(strainT3,stressT3,'k:')
p.plot(strainTC,stressTC,'k-',label=r'150$^o$C,250.,50')
p.legend(loc='lower right')
p.title(r'Temperature Change')

dtimes = np.ones(strains0.size)*strains[1]/0.01
dtimes[np.where(strains0<0.3)[0]]=strains[1]/10000.
dtimes[np.where(strains0>0.6)[0]]=strains[1]/10.
stressRC,strainRC = mm.stresscalc(displ,dtimes,25.,**V)

p.close(4)
p.figure(4)
p.plot(strainR1,stressR1,'k:',label=r'10000$^{-1}$')
p.plot(strainR2,stressR2,'k:',label=r'10$^{-1}$')
p.plot(strainR3,stressR3,'k:',label=r'0.01$^{-1}$')
p.plot(strainRC,stressRC,'k-',label=r'10,0.01,10000.$^{-1}$,')
p.legend(loc='lower right')
p.title(r'Rate Change')

