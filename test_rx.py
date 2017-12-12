import numpy as np
import matplotlib.pyplot as p
#from scipy.interpolate import interp1d
#import scipy.optimize as opt
import material_rx as mm
#import material_recrystallise as mm
from os import system,path
p.ion()

reload(mm)

savefigs = True
fignames = 'MTDB_rx'

figdir = 'figures/'+fignames+'/'
if path.exists('figures/')==False:
  system('mkdir figures/')
if path.exists(figdir)==False:
  system('mkdir '+figdir)

V = {
       # Elastic Properties:
       'E_mu0':1e5,
       'E_d0':1e3,
       'E_t0':200.,
       'E_nu':0.3,
       # Reference stress values
       'sigma_a':1.,
       'sigma_0':100.,
       # ISV's at time=0 
       '1/Lg_0':1.e-10,
       'rho_0':1.,
       # Evolution of 1/Lg:      
       'c/grainsize':40.,
       'Lm':0.7, # Exponent
       'grainsize_0':60., # Initial grain size
       # Stage IV:   
       'C00':1000., # Constant
       # Storage
       'C1':100.,
       # Dynamic recovery:
       'C20':10., # Constant
       'a02':0.5, # \frac{k_B}{g02 b^3} related to the activation energy
       'rate02':1e10,
       # Thermal recovery
       'C30':0.1, # Constant
       'C3M':2., # Exponent
       'C3U':5e3, # \frac{Q/R} related to the activation energy
       # Recrystallisation
       'C40':1.e2, # Constant
       'C4T':1.5e4, # Temperature dependence
       'C4L':10., # Misorientation Constant
       'C4D':1., # Misorientation Scaling Theta/Theta_max = min([C4D/L,1])
       'C4M':5., # Misorientation exponent
       'C4A':0.22, # Interfacial Area Exponent "a"
       'C4B':1.3, # Interfacial Area Exponent "b"
       'C4C':20., # Interfacial Area Constant "c"
       'RX_nr':10, # Number of recrystalised volume fractions to track
       'RX_sizef':1., # Recrystallised grain size fraction
       # Kinetic equation / scaling function
       'a0e':2., # \frac{k_B}{g0e b^3} related to the activation energy
       'rate0':1e7,
       'qe':1.2,
       'pe':0.5,
       # Input temperature (specify Celcius or Kelvin is assumed)
       'tempformat':'Celcius',
       # Return additional information (Evolution of Recrystallised volume fractions, 1/Lg, rho and Approximate average grain size )
       'ADDinfo':False,
       # move recrystallised volume fractions if Xv>0.999 (Xvf2-->Xvf1-->Xvf0) 
       'moveXv':False,
       # use Finite Difference Gradients?
       'useFD':True    
       }

strains = np.linspace(0,1.,51)

rate = 1.
dtime = strains[1]/rate

V['ADDinfo'] = True
stressT300,Xv300,Y300,L300,P300,D300,eqp300 = mm.stresscalc(strains,strains[1]/0.0001,400.,**V)
Y1 = Y300[:,0]*(1.-Xv300[:,1])
Y2 = Y300[:,1]*(Xv300[:,1]-Xv300[:,2])
Y3 = Y300[:,2]*(Xv300[:,2]-Xv300[:,3])
Y4 = Y300[:,3]*(Xv300[:,3]-Xv300[:,4])
Y5 = Y300[:,4]*(Xv300[:,4]-Xv300[:,5])
nrRX = Xv300.shape[1]-1
dXv300 = Xv300[:,0:nrRX]-Xv300[:,1:]
Ytotal = np.sum(dXv300*Y300,1)

fnum = 1

p.close(fnum)
p.figure(fnum)
p.plot(strains,Ytotal,'k-',label=r'$\varrho$',lw=2.,alpha=0.5)
p.plot(strains,Y1,'k-s',label=r'$\varrho_{1-X_1}(1-X_1)$',ms=4.)
p.plot(strains,Y2,'k-^',label=r'$\varrho_{X_1-X_2}(X_1-X_2)$',ms=4.)
p.plot(strains,Y3,'k-v',ms=4.)
p.plot(strains,Y4,'k-<',ms=4.)
p.plot(strains,Y5,'k->',ms=4.)
p.xlabel('Strain')
p.ylabel('Volume-Averaged Dislocation Density Ratio')
p.ylim([0,70])
p.legend(loc='upper right')
if savefigs:
  p.savefig(figdir+'f00_T300_R0001_VADens.png')
  fnum -= 1
fnum+=1


p.close(fnum)
p.figure(fnum)
p.plot(strains,Xv300[:,1],'k-s',label=r'$X_1$',ms=4.)
p.plot(strains,Xv300[:,2],'k-^',label=r'$X_2$',ms=4.)
p.plot(strains,Xv300[:,3],'k-v',label=r'$X_3$',ms=4.)
p.plot(strains,Xv300[:,4],'k-<',label=r'$X_4$',ms=4.)
p.xlabel('Strain')
p.ylabel('Recrystallised Volume Fraction')
#p.legend(loc='upper left')
if savefigs:
  p.savefig(figdir+'f00_T300_R0001_RecVolFrac.png')
  fnum -= 1
fnum+=1

p.close(fnum)
p.figure(fnum)
#lstnrs=np.where((Xv300[:,1]<0.999)&(Xv300[:,0]>0.0001))[0]
p.plot(strains,Y300[:,0],'k-s',label=r'$\varrho_{1-X_1}$',ms=4.)
#lstnrs=np.where((Xv300[:,2]<0.999)&(Xv300[:,1]>0.0001))[0]
lstnrs =np.where(Y300[:,1]>1.)[0]
if lstnrs.size>2:
  lstnrs=np.r_[lstnrs[0]-1,lstnrs]
  p.plot(strains[lstnrs],Y300[lstnrs,1],'k-^',label=r'$\varrho_{X_1-X_2}$',ms=4.)
  lstnrs =np.where(Y300[:,2]>1.)[0]
  if lstnrs.size>2:
    lstnrs=np.r_[lstnrs[0]-1,lstnrs]
    p.plot(strains[lstnrs],Y300[lstnrs,2],'k-v',label=r'$\varrho_{X_2-X_3}$',ms=4.)
    #lstnrs=np.where((Xv300[:,3]<0.999)&(Xv300[:,2]>0.0001))[0]
    lstnrs =np.where(Y300[:,3]>1.)[0]
    if lstnrs.size>2:
      lstnrs=np.r_[lstnrs[0]-1,lstnrs]
      p.plot(strains[lstnrs],Y300[lstnrs,3],'k-<',label=r'$\varrho_{X_3-X_4}$',ms=4.)
      lstnrs =np.where(Y300[:,4]>1.)[0]
      #lstnrs=np.where((Xv300[:,4]<0.999)&(Xv300[:,3]>0.0001))[0]
      if lstnrs.size>2:
        lstnrs=np.r_[lstnrs[0]-1,lstnrs]
        p.plot(strains[lstnrs],Y300[lstnrs,4],'k->',label=r'$\varrho_{X_4-X_5}$',ms=4.)
        lstnrs =np.where(Y300[:,5]>1.)[0]
        #lstnrs=np.where((Xv300[:,5]<0.999)&(Xv300[:,4]>0.0001))[0]
        if lstnrs.size>2:
          lstnrs=np.r_[lstnrs[0]-1,lstnrs]
          p.plot(strains[lstnrs],Y300[lstnrs,5],'k-d',label=r'$\varrho_{X_5-X_6}$',ms=4.)
p.xlabel('Strain')
p.ylabel('Recrystallised Volume Fraction Dislocation Density Ratio')
p.legend(loc='upper left')
if savefigs:
  p.savefig(figdir+'f00_T300_R0001_Dens.png')
  fnum -= 1
fnum+=1

p.close(fnum)
p.figure(fnum)
p.plot(strains,L300[:,0],'k-s',label=r'$\frac{1}{Lg}_{1-X_1}$',ms=4.)
lstnrs =np.where(Y300[:,1]>1.)[0]
if lstnrs.size>2:
  lstnrs=np.r_[lstnrs[0]-1,lstnrs]
  p.plot(strains[lstnrs],L300[lstnrs,1],'k-^',label=r'$\frac{1}{Lg}_{X_1-X_2}$',ms=4.)
  lstnrs =np.where(Y300[:,2]>1.)[0]
  if lstnrs.size>2:
    lstnrs=np.r_[lstnrs[0]-1,lstnrs]
    p.plot(strains[lstnrs],L300[lstnrs,2],'k-v',label=r'$\frac{1}{Lg}_{X_2-X_3}$',ms=4.)
    lstnrs =np.where(Y300[:,3]>1.)[0]
    if lstnrs.size>2:
      lstnrs=np.r_[lstnrs[0]-1,lstnrs]
      p.plot(strains[lstnrs],L300[lstnrs,3],'k-<',label=r'$\frac{1}{Lg}_{X_3-X_4}$',ms=4.)
      lstnrs =np.where(Y300[:,4]>1.)[0]
      if lstnrs.size>2:
        lstnrs=np.r_[lstnrs[0]-1,lstnrs]
        p.plot(strains[lstnrs],L300[lstnrs,4],'k->',label=r'$\frac{1}{Lg}_{X_4-X_5}$',ms=4.)
        lstnrs =np.where(Y300[:,5]>1.)[0]
        if lstnrs.size>2:
          lstnrs=np.r_[lstnrs[0]-1,lstnrs]
          p.plot(strains[lstnrs],L300[lstnrs,5],'k-d',label=r'$\frac{1}{Lg}_{X_5-X_6}$',ms=4.)
p.xlabel('Strain')
p.ylabel('Recrystallised Volume Fraction Average Lattice Incompatibility')
p.legend(loc='upper left')
if savefigs:
  p.savefig(figdir+'f00_T300_R0001_Incomp.png')
  fnum -= 1
fnum+=1


p.close(fnum)
p.figure(fnum)
p.plot(strains,D300,'k-s',label=r'$X_1$',ms=4.)
p.xlabel('Strain')
p.ylabel(r'Average grain size $\bar{D}$')
if savefigs:
  p.savefig(figdir+'f00_D300_R0001_RecXsize.png')
  fnum -= 1
fnum+=1

p.close(fnum)
p.figure(fnum)
p.plot(strains,P300[:,0],'k-s',label=r'$\varepsilon^\mathrm{p}_{1-X_1}$',ms=4.)
lstnrs =np.where(Y300[:,1]>1.)[0]
if lstnrs.size>2:
  lstnrs=np.r_[lstnrs[0]-1,lstnrs]
  p.plot(strains[lstnrs],P300[lstnrs,1],'k-^',label=r'$\varepsilon^\mathrm{p}_{X_1-X_2}$',ms=4.)
  lstnrs =np.where(Y300[:,2]>1.)[0]
  if lstnrs.size>2:
    lstnrs=np.r_[lstnrs[0]-1,lstnrs]
    p.plot(strains[lstnrs],P300[lstnrs,2],'k-v',label=r'$\varepsilon^\mathrm{p}_{X_2-X_3}$',ms=4.)
    lstnrs =np.where(Y300[:,3]>1.)[0]
    if lstnrs.size>2:
      lstnrs=np.r_[lstnrs[0]-1,lstnrs]
      p.plot(strains[lstnrs],P300[lstnrs,3],'k-<',label=r'$\varepsilon^\mathrm{p}_{X_3-X_4}$',ms=4.)
      lstnrs =np.where(Y300[:,4]>1.)[0]
      if lstnrs.size>2:
        lstnrs=np.r_[lstnrs[0]-1,lstnrs]
        p.plot(strains[lstnrs],P300[lstnrs,4],'k->',label=r'$\varepsilon^\mathrm{p}_{X_4-X_5}$',ms=4.)
        lstnrs =np.where(Y300[:,5]>1.)[0]
        if lstnrs.size>2:
          lstnrs=np.r_[lstnrs[0]-1,lstnrs]
          p.plot(strains[lstnrs],P300[lstnrs,5],'k-d',label=r'$\varepsilon^\mathrm{p}_{X_5-X_6}$',ms=4.)
p.xlabel('Strain')
p.ylabel('Plastic strain')
p.legend(loc='upper left')
if savefigs:
  p.savefig(figdir+'f00_D300_R0001_Plast.png')
  fnum -= 1
fnum+=1

p.close(fnum)
p.figure(fnum)
p.plot(strains,eqp300,'k-s',label=r'$\varepsilon^\mathrm{p}_{1-X_1}$',ms=4.)
p.xlabel('Strain')
p.ylabel('Equivalent plastic strain')
if savefigs:
  p.savefig(figdir+'f00_D300_R0001_Eqpl.png')
  fnum -= 1
fnum+=1




# #
# # ~~~~ ABOVE ONLY SHOWS VOLUME FRACTION GROWTH ETC FOR RECRYSTALLISATION AT 300DEG 0.0001/S
# ##~~~~~ ACTUAL figures 1-.. 
# V['ADDinfo'] = False
# stressT0 = mm.stresscalc(strains,strains[1]/1000.,25.,**V)
# stressT1 = mm.stresscalc(strains,strains[1]/1000.,25.,**V)
# stressT2 = mm.stresscalc(strains,strains[1]/1000.,100.,**V)
# stressT3 = mm.stresscalc(strains,strains[1]/1000.,200.,**V)
# stressT4 = mm.stresscalc(strains,strains[1]/1000.,300.,**V)
# stressT5 = mm.stresscalc(strains,strains[1]/1000.,400.,**V)
# stressT6 = mm.stresscalc(strains,strains[1]/1000.,500.,**V)
#
# p.close(fnum)
# p.figure(fnum)
# p.plot(strains,stressT1,'k-s',label=r'25$^o$C',ms=4.)
# p.plot(strains,stressT2,'k-^',label=r'100$^o$C',ms=4.)
# p.plot(strains,stressT3,'k-v',label=r'200$^o$C',ms=4.)
# p.plot(strains,stressT4,'k-<',label=r'300$^o$C',ms=4.)
# p.plot(strains,stressT5,'k->',label=r'400$^o$C',ms=4.)
# p.plot(strains,stressT6,'k-d',label=r'500$^o$C',ms=4.)
# p.xticks(np.array(range(11))/10.)
# ymax = int(p.ylim()[1]/100)+1
# p.yticks(np.array(range(ymax))*100.)
# p.grid(True)
# p.legend(loc='lower right')
# p.xlabel('Strain')
# p.ylabel('Stress [MPa]')
# p.title(r'Temperature dependence at $\dot{\varepsilon}=1000s^{-1}$')
# if savefigs:
#   p.savefig(figdir+'f01_Trates1000.png')
#   fnum -= 1
# fnum+=1
#   
# stressT0 = mm.stresscalc(strains,strains[1]/0.0001,25.,**V)
# stressT1 = mm.stresscalc(strains,strains[1]/0.0001,25.,**V)
# stressT2 = mm.stresscalc(strains,strains[1]/0.0001,100.,**V)
# stressT3 = mm.stresscalc(strains,strains[1]/0.0001,200.,**V)
# stressT4 = mm.stresscalc(strains,strains[1]/0.0001,300.,**V)
# stressT5 = mm.stresscalc(strains,strains[1]/0.0001,400.,**V)
# stressT6 = mm.stresscalc(strains,strains[1]/0.0001,500.,**V)
#
# p.close(fnum)
# p.figure(fnum)
# p.plot(strains,stressT1,'k-s',label=r'25$^o$C',ms=4.)
# p.plot(strains,stressT2,'k-^',label=r'100$^o$C',ms=4.)
# p.plot(strains,stressT3,'k-v',label=r'200$^o$C',ms=4.)
# p.plot(strains,stressT4,'k-<',label=r'300$^o$C',ms=4.)
# p.plot(strains,stressT5,'k->',label=r'400$^o$C',ms=4.)
# p.plot(strains,stressT6,'k-d',label=r'500$^o$C',ms=4.)
# p.xticks(np.array(range(11))/10.)
# ymax = int(p.ylim()[1]/100)+1
# p.yticks(np.array(range(ymax))*100.)
# p.grid(True)
# p.legend(loc='lower right')
# p.xlabel('Strain')
# p.ylabel('Stress [MPa]')
# p.title(r'Temperature dependence at $\dot{\varepsilon}=0.0001s^{-1}$')
# if savefigs:
#   p.savefig(figdir+'f01_Trates0001.png')
#   fnum -= 1
# fnum+=1
#   
#   
#
# stressT0 = mm.stresscalc(strains,strains[1]/1,25.,**V)
# stressT1 = mm.stresscalc(strains,strains[1]/1,25.,**V)
# stressT2 = mm.stresscalc(strains,strains[1]/1,100.,**V)
# stressT3 = mm.stresscalc(strains,strains[1]/1,200.,**V)
# stressT4 = mm.stresscalc(strains,strains[1]/1,300.,**V)
# stressT5 = mm.stresscalc(strains,strains[1]/1,400.,**V)
# stressT6 = mm.stresscalc(strains,strains[1]/1,500.,**V)
#
# p.close(fnum)
# p.figure(fnum)
# p.plot(strains,stressT1,'k-s',label=r'25$^o$C',ms=4.)
# p.plot(strains,stressT2,'k-^',label=r'100$^o$C',ms=4.)
# p.plot(strains,stressT3,'k-v',label=r'200$^o$C',ms=4.)
# p.plot(strains,stressT4,'k-<',label=r'300$^o$C',ms=4.)
# p.plot(strains,stressT5,'k->',label=r'400$^o$C',ms=4.)
# p.plot(strains,stressT6,'k-d',label=r'500$^o$C',ms=4.)
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
#   fnum -= 1
# fnum+=1
#   
#   
#
#
#
# stressR0 = mm.stresscalc(strains,strains[1]/0.0001,25.,**V)
# stressR1 = mm.stresscalc(strains,strains[1]/10000.,25.,**V)
# stressR2 = mm.stresscalc(strains,strains[1]/10,25.,**V)
# stressR3 = mm.stresscalc(strains,strains[1]/0.01,25.,**V)
# stressR4 = mm.stresscalc(strains,strains[1]/0.001,25.,**V)
# stressR5 = mm.stresscalc(strains,strains[1]/0.0001,25.,**V)
# stressR6 = mm.stresscalc(strains,strains[1]/0.00001,25.,**V)
#
# p.close(fnum)
# p.figure(fnum)
# p.plot(strains,stressR1,'k-s',label=r'10000$s^{-1}$',ms=4.)
# p.plot(strains,stressR2,'k-^',label=r'10$s^{-1}$',ms=4.)
# p.plot(strains,stressR3,'k-v',label=r'0.01$s^{-1}$',ms=4.)
# #p.plot(strains,stressR4,'k<',label=r'0.001$s^{-1}$',ms=4.)
# #p.plot(strains,stressR5,'k>',label=r'0.0001$s^{-1}$',ms=4.)
# p.plot(strains,stressR6,'k-d',label=r'0.00001$s^{-1}$',ms=4.)
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
#   fnum -= 1
# fnum+=1
#
# #
# stressR00 = mm.stresscalc(strains,strains[1]/0.0001,250.,**V)
# stressR01 = mm.stresscalc(strains,strains[1]/10000.,250.,**V)
# stressR02 = mm.stresscalc(strains,strains[1]/10,250.,**V)
# stressR03 = mm.stresscalc(strains,strains[1]/0.01,250.,**V)
# stressR04 = mm.stresscalc(strains,strains[1]/0.001,250.,**V)
# stressR05 = mm.stresscalc(strains,strains[1]/0.0001,250.,**V)
# stressR06 = mm.stresscalc(strains,strains[1]/0.00001,250.,**V)
#
# p.close(fnum)
# p.figure(fnum)
# p.plot(strains,stressR01,'k-s',label=r'10000$s^{-1}$',ms=4.)
# p.plot(strains,stressR02,'k-^',label=r'10$s^{-1}$',ms=4.)
# p.plot(strains,stressR03,'k-v',label=r'0.01$s^{-1}$',ms=4.)
# p.plot(strains,stressR04,'k-<',label=r'0.001$s^{-1}$',ms=4.)
# p.plot(strains,stressR05,'k->',label=r'0.0001$s^{-1}$',ms=4.)
# p.plot(strains,stressR06,'k-d',label=r'0.00001$s^{-1}$',ms=4.)
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
#   fnum -= 1
# fnum+=1
#
#
#
#
#
#
#
#
# #change in temperature:
# temps=np.ones(strains.size)*25.
# temps[np.where(strains<0.3)[0]]=200.
# temps[np.where(strains>0.6)[0]]=500.
# stressTC = mm.stresscalc(strains,strains[1]/1,temps,**V)
#
# p.close(fnum)
# p.figure(fnum)
# p.plot(strains,stressT1,'k-',alpha=0.3)
# p.plot(strains,stressT6,'k-',alpha=0.3)
# p.plot(strains,stressT3,'k-',alpha=0.3)
# p.plot(strains,stressTC,'ks',label=r'150$^o$C,250.,50',ms=4.)
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
#   fnum -= 1
# fnum+=1
#
# dtimes = np.ones(strains.size)*strains[1]/0.00001
# dtimes[np.where(strains<0.3)[0]]=strains[1]/10000.
# dtimes[np.where(strains>0.6)[0]]=strains[1]/10.
# stressRC = mm.stresscalc(strains,dtimes,25.,**V)
#
# p.close(fnum)
# p.figure(fnum)
# p.plot(strains,stressR1,'k-',label=r'10000$s^{-1}$',alpha=0.3)
# p.plot(strains,stressR2,'k-',label=r'10$s^{-1}$',alpha=0.3)
# p.plot(strains,stressR6,'k-',label=r'0.01$s^{-1}$',alpha=0.3)
# p.plot(strains,stressRC,'ks',label=r'10,0.01,10000.$s^{-1}$,',ms=4.)
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
#   fnum -= 1
# fnum+=1
#
#
#
# ltp5 = np.where(strains<=0.5)[0]
# mtp5 = np.where(strains>0.5)[0]
# lst1 = np.r_[ltp5,mtp5+1]
# lst5 = np.r_[ltp5,mtp5+5]
# lst20 = np.r_[ltp5,mtp5+20]
# strainsh1 = np.r_[strains[ltp5],0.5,strains[mtp5]]
# strainsh5 = np.r_[strains[ltp5],0.5*np.ones(5),strains[mtp5]]
# strainsh20 = np.r_[strains[ltp5],0.5*np.ones(20),strains[mtp5]]
#
# T1,T2 = 25,300
#
# rate = 0.0001
# dtime = strains[1]/rate
# dtimes0 = np.r_[np.ones(ltp5.size)*dtime,np.ones(mtp5.size)*dtime]
# dtimes1 = np.r_[np.ones(ltp5.size)*dtime,60.,np.ones(mtp5.size)*dtime]
# dtimes2 = np.r_[np.ones(ltp5.size)*dtime,4*60.*np.ones(5),np.ones(mtp5.size)*dtime]
# dtimes3 = np.r_[np.ones(ltp5.size)*dtime,3*60.*24*np.ones(20),np.ones(mtp5.size)*dtime]
# temps0 = np.r_[np.ones(ltp5.size)*T1,np.ones(mtp5.size)*T2]
# temps1 = np.r_[np.ones(ltp5.size)*T1,T2,np.ones(mtp5.size)*T2]
# temps2 = np.r_[np.ones(ltp5.size)*T1,T2*np.ones(5),np.ones(mtp5.size)*T2]
# temps3 = np.r_[np.ones(ltp5.size)*T1,T2*np.ones(20),np.ones(mtp5.size)*T2]
# stress0min = mm.stresscalc(strains,dtimes0,temps0,**V)
# stress1min = mm.stresscalc(strainsh1,dtimes1,temps1,**V)
# stress20min = mm.stresscalc(strainsh5,dtimes2,temps2,**V)
# stress60min = mm.stresscalc(strainsh20,dtimes3,temps3,**V)
# stress0u = mm.stresscalc(strains,strains[1]/0.0001,T1,**V)
# stress0l = mm.stresscalc(strains,strains[1]/0.0001,T2,**V)
#
# p.close(fnum)
# p.figure(fnum)
# p.plot(strains,stress0u,'k-',alpha=0.3)
# p.plot(strains,stress0l,'k-',alpha=0.3)
# #p.plot(strains,stress0min,'k',label='no hold')
# p.plot(strainsh1[lst1],stress1min[lst1],'ks',label='1 min hold',ms=4.)
# p.plot(strainsh5[lst5],stress20min[lst5],'k^',label='20 min hold',ms=4)
# p.plot(strainsh20[lst20],stress60min[lst20],'kv',label='24 hrs hold',ms=4.)
# p.title(r'%s$^o$C then %s$^o$C at $\dot{\varepsilon}=0.0001s^{-1}$'%(T1,T2))
# p.xticks(np.array(range(11))/10.)
# ymax = int(p.ylim()[1]/100)+1
# p.yticks(np.array(range(ymax))*100.)
# p.grid(True)
# p.xlabel('Strain')
# p.ylabel('Stress [MPa]')
# p.legend(loc='lower right')
# if savefigs:
#   p.savefig(figdir+'f06_holdtime.png')
#   fnum -= 1
# fnum+=1
#
#
#
#
#
#
#
# # strains0 = np.linspace(0,0.5,26)
# # straincyc = np.r_[strains0,strains0[-2::-1],-strains0[1:],-strains0[-2::-1]]#,strains0[1:]]
# # strains = np.r_[straincyc,straincyc[1:],strains0[1:]]#,straincyc[1:]]
# # strains = np.exp(strains)-1
# #
# # stresscyc1,straincyc1 = mm.stresscalc(strains,strains[1]/1,25.,**V)
# # stresscyc2,straincyc2 = mm.stresscalc(strains,strains[1]/0.0001,250.,**V)
# # stresscyc3,straincyc3 = mm.stresscalc(strains,strains[1]/1000,25.,**V)
# #
# # p.close(fnum)
# # p.figure(fnum)
# # p.plot(straincyc1,stresscyc1,'k-',alpha=0.3)
# # p.plot(straincyc1,stresscyc1,'ks',ms=4.)
# # p.xlabel('Strain')
# # p.ylabel('Stress [MPa]')
# # p.title(r'Cyclic behaviour at $\dot{\varepsilon}=1s^{-1}$ 25$^o$C')
# # p.xticks(np.linspace(-0.6,0.6,13))
# # ymax = int(p.ylim()[1]/100)
# # #p.yticks((np.array(range(2*ymax+1))-ymax)*100.)
# # p.grid(True)
# # if savefigs:
# #  p.savefig(figdir+'f07_cycle1.png')
# #   fnum -= 1
# # fnum+=1
# #
# # p.close(fnum)
# # p.figure(fnum)
# # p.plot(straincyc2,stresscyc2,'k-',alpha=0.3)
# # p.plot(straincyc2,stresscyc2,'ks',ms=4.)
# # p.xlabel('Strain')
# # p.ylabel('Stress [MPa]')
# # p.title(r'Cyclic behaviour at $\dot{\varepsilon}=0.0001s^{-1}$ 200$^o$C')
# # p.xticks(np.linspace(-0.6,0.6,13))
# # ymax = int(p.ylim()[1]/100)
# # #p.yticks((np.array(range(2*ymax+1))-ymax)*100.)
# # p.grid(True)
# # if savefigs:
# #  p.savefig(figdir+'f08_cycle2.png')
# #
# # p.close(fnum)
# # p.figure(fnum)
# # p.plot(straincyc3,stresscyc3,'k-',alpha=0.3)
# # p.plot(straincyc3,stresscyc3,'ks',ms=4.)
# # p.xlabel('Strain')
# # p.ylabel('Stress [MPa]')
# # p.title(r'Cyclic behaviour at $\dot{\varepsilon}=1000s^{-1}$ 25$^o$C')
# # p.xticks(np.linspace(-0.6,0.6,13))
# # ymax = int(p.ylim()[1]/100)
# # #p.yticks((np.array(range(2*ymax+1))-ymax)*100.)
# # p.grid(True)
# # if savefigs:
# #  p.savefig(figdir+'f09_cycle3.png')