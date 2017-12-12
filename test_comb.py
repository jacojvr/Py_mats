# -*- coding: utf-8 -*-
import numpy as np
#datadir = './TannerData/'
import matplotlib.pyplot as p
from scipy.interpolate import interp1d
import scipy.optimize as opt
import material_comb as mm
p.ion()

reload(mm)

V = { 'emod':1000., # Young's modulus
        'sa':50., # Initial yield stress
        'dsi':0., # Linear sotropic hardening modulus
	'dsk':250., # Linear kinematic hardening modulkus
	}

strains0 = np.linspace(0,1.,51)
straincyc = np.r_[strains0,strains0[-2::-1],-strains0[1:],-strains0[-2::-1]]
strains = np.r_[straincyc,straincyc[1:],straincyc[1:]]#,straincyc[1:]]
displ = np.exp(strains)-1

stressi,straini = mm.stresscalc(displ,50.,25.,**{'dsi':250.,'dsk':0.})
stressk,straink = mm.stresscalc(displ,50.,25.,**{'dsi':0.,'dsk':250.})
stressc,strainc = mm.stresscalc(displ,50.,25.,**{'dsi':125.,'dsk':125.})

p.close(1)
p.figure(1)
p.plot(straini,stressi,'k--',label='Isotripic')
p.plot(straink,stressk,'k-.',label='Kinematic')
p.plot(strainc,stressc,'k-',label='Combined')
p.legend(loc='lower right')
p.title('Combined Isotropic-Kinematic Hardening')
