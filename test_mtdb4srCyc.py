# -*- coding: utf-8 -*-
import numpy as np
#datadir = './TannerData/'
import matplotlib.pyplot as p
from scipy.interpolate import interp1d
import scipy.optimize as opt
import material_mtdb4statcyc as mm
from os import system,path
p.ion()

reload(mm)

savefigs = False
fignames = 'MTDBcyc'

figdir = 'figures/'+fignames+'/'
if path.exists('figures/')==False:
  system('mkdir figures/')
if path.exists(figdir)==False:
  system('mkdir '+figdir)

V = {'mu0':100000.,
       'd0':10000.,
       't0':200.,
       'nu':0.3,
       'sa':10.,
       's0':100.,
       'Y0':1.,
       'C1':100.,
       'C20':10.,
       'lam0':0.,
       'cD':10.,
       'C0':0,#10.,#100
       'r':1.,
       'si':0.,
       'a0e':2.,
       'a02':.5,
       'a0i':0.9,
       'rate0':1e7,
       'rate02':1e10,
       'qe':1.2,
       'pe':0.5,
       'qi':1.5,
       'pi':2./3,
       'r0':0.1, #0.1
       'beta':1,
       'U0':5000.,
       'Aback':.00,
       'Qback':0.,#30.,
       'mback':0.5,
       'Cback1':0,#100.,
       'Cback2':0,#1.,
       'C5':100.,
       'C6':0.,       
       'tempformat':'Celcius'}
       
V={'mu0':100000.,
     'd0':10000.,
     't0':200.,
     'nu':0.3,
     'sa':10.,
     'sig0':100.,
     'a0e':2.,
     'C0':10.,
     'cD':10.,
     'C0m':1.,
     'C1':100.,
     'C20':10.,
     'a02':2.,
     'rate0':1e7,
     'rate02':1e10,
     'qe':1.,
     'pe':0.5,
     'C30':0.5,
     'C3m':0.1,
     'C3U':5000.,
     'C4':0,#100.,
     'C50':0.5,
     'C5U':5000.,
     'qC1':0,#50.,
     'qC2':0.,
     'tempformat':'Celcius',
     'Y0':1.,
     'lam0':0.}

# strains0 = np.linspace(0,1.,51)
# straincyc = np.r_[strains0,strains0[-2::-1],-strains0[1:],-strains0[-2::-1]]
# strains = np.r_[straincyc,straincyc[1:],straincyc[1:]]#,straincyc[1:]]
# displ = np.exp(strains0)-1
#
# rate = 1.
# dtime = strains[1]/rate
#
# stressT0,strainT0 = mm.stresscalc(displ,strains[1]/1,25.,**V)
# stressT1,strainT1 = mm.stresscalc(displ,strains[1]/1,25.,**V)
# stressT2,strainT2 = mm.stresscalc(displ,strains[1]/1,100.,**V)
# stressT3,strainT3 = mm.stresscalc(displ,strains[1]/1,200.,**V)
# stressT4,strainT4 = mm.stresscalc(displ,strains[1]/1,300.,**V)
# stressT5,strainT5 = mm.stresscalc(displ,strains[1]/1,400.,**V)
# stressT6,strainT6 = mm.stresscalc(displ,strains[1]/1,500.,**V)
#
# p.close(1)
# p.figure(1)
# p.plot(strainT1,stressT1,'ks',label=r'25$^o$C',ms=4.)
# p.plot(strainT2,stressT2,'k^',label=r'100$^o$C',ms=4.)
# p.plot(strainT3,stressT3,'kv',label=r'200$^o$C',ms=4.)
# p.plot(strainT4,stressT4,'k<',label=r'300$^o$C',ms=4.)
# p.plot(strainT5,stressT5,'k>',label=r'400$^o$C',ms=4.)
# p.plot(strainT6,stressT6,'kd',label=r'500$^o$C',ms=4.)
# p.xticks(np.array(range(11))/10.)
# ymax = int(p.ylim()[1]/100)+1
# p.yticks(np.array(range(ymax))*100.)
# p.grid(True)
# p.legend(loc='lower right')
# p.xlabel('Strain')
# p.ylabel('Stress [MPa]')
# p.title(r'Temperature dependence at $\dot{\varepsilon}=1s^{-1}$')
# if savefigs:
#   p.savefig(figdir+'f01_Temperatures.png')
#
#
# stressR0,strainR0 = mm.stresscalc(displ,strains[1]/0.0001,25.,**V)
# stressR1,strainR1 = mm.stresscalc(displ,strains[1]/10000.,25.,**V)
# stressR2,strainR2 = mm.stresscalc(displ,strains[1]/10,25.,**V)
# stressR3,strainR3 = mm.stresscalc(displ,strains[1]/0.01,25.,**V)
# stressR4,strainR4 = mm.stresscalc(displ,strains[1]/0.001,25.,**V)
# stressR5,strainR5 = mm.stresscalc(displ,strains[1]/0.0001,25.,**V)
# stressR6,strainR6 = mm.stresscalc(displ,strains[1]/0.00001,25.,**V)
#
# p.close(2)
# p.figure(2)
# p.plot(strainR1,stressR1,'ks',label=r'10000$s^{-1}$',ms=4.)
# p.plot(strainR2,stressR2,'k^',label=r'10$s^{-1}$',ms=4.)
# p.plot(strainR3,stressR3,'kv',label=r'0.01$s^{-1}$',ms=4.)
# #p.plot(strainR4,stressR4,'k<',label=r'0.001$s^{-1}$',ms=4.)
# #p.plot(strainR5,stressR5,'k>',label=r'0.0001$s^{-1}$',ms=4.)
# p.plot(strainR6,stressR6,'kd',label=r'0.00001$s^{-1}$',ms=4.)
# p.xticks(np.array(range(11))/10.)
# ymax = int(p.ylim()[1]/100)+1
# p.yticks(np.array(range(ymax))*100.)
# p.grid(True)
# p.xlabel('Strain')
# p.ylabel('Stress [MPa]')
# p.legend(loc='lower right')
# p.title(r'Rate dependence at 25$^o$C')
# if savefigs:
#   p.savefig(figdir+'f02_Rates025.png')
#
#
# stressR00,strainR0 = mm.stresscalc(displ,strains[1]/0.0001,250.,**V)
# stressR01,strainR01 = mm.stresscalc(displ,strains[1]/10000.,250.,**V)
# stressR02,strainR02 = mm.stresscalc(displ,strains[1]/10,250.,**V)
# stressR03,strainR03 = mm.stresscalc(displ,strains[1]/0.01,250.,**V)
# stressR04,strainR04 = mm.stresscalc(displ,strains[1]/0.001,250.,**V)
# stressR05,strainR05 = mm.stresscalc(displ,strains[1]/0.0001,250.,**V)
# stressR06,strainR06 = mm.stresscalc(displ,strains[1]/0.00001,250.,**V)
#
# p.close(3)
# p.figure(3)
# p.plot(strainR01,stressR01,'ks',label=r'10000$s^{-1}$',ms=4.)
# p.plot(strainR02,stressR02,'k^',label=r'10$s^{-1}$',ms=4.)
# p.plot(strainR03,stressR03,'kv',label=r'0.01$s^{-1}$',ms=4.)
# p.plot(strainR04,stressR04,'k<',label=r'0.001$s^{-1}$',ms=4.)
# p.plot(strainR05,stressR05,'k>',label=r'0.0001$s^{-1}$',ms=4.)
# p.plot(strainR06,stressR06,'kd',label=r'0.00001$s^{-1}$',ms=4.)
# p.xticks(np.array(range(11))/10.)
# ymax = int(p.ylim()[1]/100)+1
# p.yticks(np.array(range(ymax))*100.)
# p.grid(True)
# p.xlabel('Strain')
# p.ylabel('Stress [MPa]')
# p.legend(loc='lower right')
# p.title(r'Rate dependence at 250$^o$C')
# if savefigs:
#   p.savefig(figdir+'f03_Rates250.png')
#
#
#
# #change in temperature:
# temps=np.ones(strains0.size)*25.
# temps[np.where(strains0<0.3)[0]]=200.
# temps[np.where(strains0>0.6)[0]]=500.
# stressTC,strainTC = mm.stresscalc(displ,strains[1]/1,temps,**V)
#
# p.close(4)
# p.figure(4)
# p.plot(strainT1,stressT1,'k-',alpha=0.3)
# p.plot(strainT6,stressT6,'k-',alpha=0.3)
# p.plot(strainT3,stressT3,'k-',alpha=0.3)
# p.plot(strainTC,stressTC,'ks',label=r'150$^o$C,250.,50',ms=4.)
# p.xticks(np.array(range(11))/10.)
# ymax = int(p.ylim()[1]/100)+1
# p.yticks(np.array(range(ymax))*100.)
# p.grid(True)
# p.xlabel('Strain')
# p.ylabel('Stress [MPa]')
# #p.legend(loc='lower right')
# p.title(r'200$^o$C then 25$^o$C then 500$^o$C at $\dot{\varepsilon}=1s^{-1}$')
# if savefigs:
#   p.savefig(figdir+'f04_change_temp.png')
#
# dtimes = np.ones(strains0.size)*strains[1]/0.00001
# dtimes[np.where(strains0<0.3)[0]]=strains[1]/10000.
# dtimes[np.where(strains0>0.6)[0]]=strains[1]/10.
# stressRC,strainRC = mm.stresscalc(displ,dtimes,25.,**V)
#
# p.close(5)
# p.figure(5)
# p.plot(strainR1,stressR1,'k-',label=r'10000$s^{-1}$',alpha=0.3)
# p.plot(strainR2,stressR2,'k-',label=r'10$s^{-1}$',alpha=0.3)
# p.plot(strainR6,stressR6,'k-',label=r'0.01$s^{-1}$',alpha=0.3)
# p.plot(strainRC,stressRC,'ks',label=r'10,0.01,10000.$s^{-1}$,',ms=4.)
# p.xticks(np.array(range(11))/10.)
# ymax = int(p.ylim()[1]/100)+1
# p.yticks(np.array(range(ymax))*100.)
# p.grid(True)
# p.xlabel('Strain')
# p.ylabel('Stress [MPa]')
# #p.legend(loc='lower right')
# p.title(r'$\dot{\varepsilon}=10000s^{-1}$ then $0.00001s^{-1}$ then $10s^{-1}$ at 25$^o$C')
# if savefigs:
#   p.savefig(figdir+'f05_change_rate.png')
#
#
# ltp5 = np.where(strains0<=0.5)[0]
# mtp5 = np.where(strains0>0.5)[0]
# lst = np.r_[ltp5,mtp5+1]
# strainrelax = np.r_[strains0[ltp5],0.5,strains0[mtp5]]
# displh = np.exp(strainrelax)-1.
# rate = 0.0001
# dtime = strains0[1]/rate
# dtimes0 = np.r_[np.ones(ltp5.size)*dtime,np.ones(mtp5.size)*dtime]
# dtimes1 = np.r_[np.ones(ltp5.size)*dtime,60.,np.ones(mtp5.size)*dtime]
# dtimes2 = np.r_[np.ones(ltp5.size)*dtime,20*60.,np.ones(mtp5.size)*dtime]
# dtimes3 = np.r_[np.ones(ltp5.size)*dtime,60*60.*24,np.ones(mtp5.size)*dtime]
# temps0 = np.r_[np.ones(ltp5.size)*25,np.ones(mtp5.size)*250]
# temps = np.r_[np.ones(ltp5.size)*25,250.,np.ones(mtp5.size)*250]
# stress0min,strain0min = mm.stresscalc(displ,dtimes0,temps0,**V)
# stress1min,strain1min = mm.stresscalc(displh,dtimes1,temps,**V)
# stress20min,strain20min = mm.stresscalc(displh,dtimes2,temps,**V)
# stress60min,strain60min = mm.stresscalc(displh,dtimes3,temps,**V)
# stress0u,strain0u = mm.stresscalc(displ,strains0[1]/0.0001,25.,**V)
# stress0l,strain0l = mm.stresscalc(displ,strains0[1]/0.0001,250.,**V)
# p.close(6)
# p.figure(6)
# p.plot(strain0u,stress0u,'k-',alpha=0.3)
# p.plot(strain0l,stress0l,'k-',alpha=0.3)
# #p.plot(strain0min,stress0min,'k',label='no hold')
# p.plot(strain1min[lst],stress1min[lst],'ks',label='1 min hold',ms=4.)
# p.plot(strain20min[lst],stress20min[lst],'k^',label='20 min hold',ms=4)
# p.plot(strain60min[lst],stress60min[lst],'kv',label='24 hrs hold',ms=4.)
# p.title(r'25$^o$C then 250$^o$C at $\dot{\varepsilon}=0.0001s^{-1}$')
# p.xticks(np.array(range(11))/10.)
# ymax = int(p.ylim()[1]/100)+1
# p.yticks(np.array(range(ymax))*100.)
# p.grid(True)
# p.xlabel('Strain')
# p.ylabel('Stress [MPa]')
# p.legend(loc='lower right')
# if savefigs:
#  p.savefig(figdir+'f06_holdtime.png')



strains0 = np.linspace(0,0.5,26)
strainslunload = np.r_[strains0,strains0[-2::-1]]
straincyc = np.r_[strains0,strains0[-2::-1],-strains0[1:],-strains0[-2::-1]]#,strains0[1:]]
strains = np.r_[straincyc,straincyc[1:],strains0[1:]]#,straincyc[1:]]
displ = np.exp(strainslunload)-1
displ = np.exp(strains)-1


stresscyc1,straincyc1 = mm.stresscalc(displ,strains[1]/1,25.,**V)
# stresscyc2,straincyc2 = mm.stresscalc(displ,strains[1]/0.0001,250.,**V)
# stresscyc3,straincyc3 = mm.stresscalc(displ,strains[1]/1000,450.,**V)
p.close(7)
p.figure(7)
p.plot(straincyc1,stresscyc1,'k-',alpha=0.3)
p.plot(straincyc1,stresscyc1,'ks',ms=4.)
p.xlabel('Strain')
p.ylabel('Stress [MPa]')
p.title(r'Cyclic behaviour at $\dot{\varepsilon}=1s^{-1}$ 25$^o$C')
#p.xticks(np.linspace(-0.6,0.6,13))
ymax = int(p.ylim()[1]/100)
#p.yticks((np.array(range(2*ymax+1))-ymax)*100.)
p.grid(True)
if savefigs:
 p.savefig(figdir+'f07_cycle1.png')
# p.close(8)
# p.figure(8)
# p.plot(straincyc2,stresscyc2,'k-',alpha=0.3)
# p.plot(straincyc2,stresscyc2,'ks',ms=4.)
# p.xlabel('Strain')
# p.ylabel('Stress [MPa]')
# p.title(r'Cyclic behaviour at $\dot{\varepsilon}=0.0001s^{-1}$ 200$^o$C')
# p.xticks(np.linspace(-0.6,0.6,13))
# ymax = int(p.ylim()[1]/100)
# #p.yticks((np.array(range(2*ymax+1))-ymax)*100.)
# p.grid(True)
# if savefigs:
#  p.savefig(figdir+'f08_cycle2.png')
# p.close(9)
# p.figure(9)
# p.plot(straincyc3,stresscyc3,'k-',alpha=0.3)
# p.plot(straincyc3,stresscyc3,'ks',ms=4.)
# p.xlabel('Strain')
# p.ylabel('Stress [MPa]')
# p.title(r'Cyclic behaviour at $\dot{\varepsilon}=1000s^{-1}$ 450$^o$C')
# p.xticks(np.linspace(-0.6,0.6,13))
# ymax = int(p.ylim()[1]/100)
# #p.yticks((np.array(range(2*ymax+1))-ymax)*100.)
# p.grid(True)
# if savefigs:
#  p.savefig(figdir+'f09_cycle3.png')