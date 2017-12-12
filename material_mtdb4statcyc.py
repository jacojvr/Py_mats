import numpy as np

def stresscalc(displ=np.linspace(0,0.5,101),dtime=1.,temp=25.,**kwargs):
  """ Mechanical Threshold Stress material:
  
  The function takes displ, time step sizes and temperature values into account
  to calculate the relevant stress state history
  
  INPUT:
     displ : a one dimensional array or list containing strain histories
             : default - numpy linspace array from 0. to 1. in 0.01 increments
     dtime   : a one dimensional array or list containing the time step increments
               associated with the strain history.
             or
	       a single value used over the entire history
	     : default - 1.
     temp    : a one dimensional array or list containing the temperature values
               associated with the strain history.
             or
	       a single value used over the entire history
	     : default - 25^oC
     ** key word arguments and default values:
     
       'mu0' = 85000.
        'd0' = 8300.
        't0' = 200.
        'nu' = 0.33
        'sa' = 35.
	's0' = 100.
        'Y0' = 1.
	'C1' = 100.
       'C20' = 10.
       'lam0' = 0.
        'cD' = 50.
        'C0' = 100.
         'r' = 1.
        'si' = 0.
       'a0e' = 2.
       'a02' = 0.55
       'a0i' = 0.9
     'rate0' = 1e7
    'rate02' = 1e10
        'qe' = 1.2
        'pe' = 0.5
        'qi' = 1.5
        'pi' = 2./3
'tempformat' = 'Celcius'

  # Dictionary containing variables and default values        
#   V = {'mu0':100000.,
#        'd0':10000.,
#        't0':200.,
#        'nu':0.3,
#        'sa':10.,
#        's0':100.,
#        'Y0':1.,
#        'C1':100.,
#        'C20':10.,
#        'lam0':0.,
#        'cD':10.,
#        'C0':100.,#100
#        'r':1.,
#        'si':0.,
#        'a0e':2.,
#        'a02':.5,
#        'a0i':0.9,
#        'rate0':1e7,
#        'rate02':1e10,
#        'qe':1.2,
#        'pe':0.5,
#        'qi':1.5,
#        'pi':2./3,
#        'r0':10., #0.1
#        'beta':1.,
#        'U0':5000.,
#        'Aback':100.,
#        'Qback':30.,
#        'mback':0.5,
#        'Cback1':100.,
#        'Cback2':10.,
#        'C5':100.,
#        'C6':0.0001,       
#        'tempformat':'Celcius'}
       
  """
  
  V={'mu0':100000.,
     'd0':10000.,
     't0':200.,
     'nu':0.3,
     'sa':10.,
     'sig0':100.,
     'a0e':2.,
     'C0':0,#10.,
     'cD':10.,
     'C0m':1.,
     'C1':100.,
     'C20':10.,
     'a02':2.,
     'rate0':1e7,
     'rate02':1e10,
     'qe':1.,
     'pe':0.5,
     'C30':0.,#0.5,
     'C3m':0.1,
     'C3U':5000.,
     'C4':100.,
     'C50':0.5,
     'C5U':5000.,
     'qC1':50.,
     'qC2':0.,
     'tempformat':'Celcius',
     'Y0':1.,
     'lam0':0.}
	  
  toler = 1e-12
  # Check key word arguments defined during function call and set values in variable dictionary
  givedef = 0.
#   if kwargs.__len__()>0:
#     for vardef in kwargs.keys():
#       if vardef in V.keys():
# 	V[vardef] = kwargs[vardef]
#       else:
#         givedef = 1
#         print "\n\t*** WARNING *** keyword - "+vardef+" - unknown."
#   if givedef:
#     print "\n ACCEPTED LIST OF KEYWORDS:\n"
#     for vardef in V.keys():
#       print "\t"+vardef+" with default value: ",V[vardef]
	
  # Initialise equivalent plastic strain and kinematic plasticity
  eqpl = 0.
  sykp = 0. 
  Yp = V['Y0']
  Lp = V['lam0']
  Zp = 0.
  # Make sure inputs are consistent:
  displ = np.array(displ,float)
  dtime = np.array(dtime,float)
  temp = np.array(temp,float)
  nrs = displ.size
  if dtime.size<nrs:
    if dtime.size>1:
      print "\n\t* Time step array size inconsistent with displacements, using only the first value: ",dtime[0]
      dtime = dtime[0]
    dtime = np.ones(nrs)*dtime
  if temp.size<nrs:
    if temp.size>1:
      print "\n\t* Temperature array size inconsistent with displacement, using only the first value : ",temp[0]
      temp = temp[0]
    temp = np.ones(nrs)*temp
        
  # Remove possible zero strain at the start:
  stress = []
  plastic = []
  elastic = []
  strain = []
  eqplhist = []
  logstrainprev = 0.
  elasticprev=0.
  plasticprev=0.
  eqpl = 0.
  backstressprev=0.
  stressprev = 0.
  alpha = 0.
  flowprev = 1.
  if displ[0]==0:
    nrs -=1
    stress = [0.]
    strain = [0.]
    plastic = [0.]
    elastic =[0.]
    eqplhist =[0.]
    displ = displ[1::]
    dtime = dtime[1::]
    temp = temp[1::]
    if V['tempformat']=='Celcius':
      temp+=273.
    
  # Constants
  c1,c2 = 2./3, np.sqrt(2./3)
  # Begin solution:
  # Get variables:
  [mu0,d0,t0,nu] = [V[kn] for kn in ['mu0','d0','t0','nu']]
  
  for k in range(nrs):
    #print "*** INCREMENT ",k+1
    logstrain = np.log(1.+displ[k])
    dstrain = logstrain-logstrainprev
    logstrainprev = logstrain
    # Elastic guess for strain and stress increment:
    elasticcur = elasticprev+dstrain
    deqpl = 0.
    plasticcur = plasticprev
    #
    # current temperature dependent shear modulus:
    if temp[k]>t0:
      mu = mu0 - d0/(np.exp(t0/temp[k])-1.)
    else:
      mu = mu0
    # elastic modulus:
    emod = 2.*mu*(1.+nu)
    
    stresscur = emod*elasticcur
    
     
    smises = np.abs(stresscur-alpha)
    flowdir = (stresscur-alpha)/smises
    #flowdir = stresscur/np.abs(stresscur)
    # Absolute value
    
    # stress zero or stress reversal: wall dislocations back into forest
    if np.abs(smises)<toler:
      Yp = np.max([Yp-Zp,1])
      Zp = 0.
    elif (flowdir==flowprev)==False:
      Yp = np.max([Yp-Zp,1])
      Zp = 0.
    flowprev = flowdir
    
    
    
        
    syield0,dsi,sk,dsk,Yc,Lc,Zc=mtdb(Yp,Lp,Zp,backstressprev,toler,dtime[k],temp[k-1],temp[k]-temp[k-1],V)
    syield0fd,dsifd,skfd,dskfd,Ycfd,Lcfd,Zcfd=mtdb(Yp,Lp,Zp,backstressprev,toler+toler,dtime[k],temp[k-1],temp[k]-temp[k-1],V)
    dsi = (syield0fd-syield0)/toler
    dsk = (skfd-sk)/toler
    Yc = Yp
    Lc = Lp
    Zc = Zp 
    backstresscur=backstressprev
        
    f = smises - syield0
    
    if f>toler:
      syield = syield0
      counter = 0.
      while (np.abs(f)>toler)&(counter<10):
        counter+=1
	f = smises-emod*deqpl-syield-backstresscur+backstressprev
        deqpl=deqpl+f/(emod+dsi+dsk)
#         f = smises-emod*deqpl-syield
#         deqpl=deqpl+f/(emod+dsi)
        deqpl = np.abs(deqpl)
	syield,dsiA,sk,dskA,Yc,Lc,Zc=mtdb(Yp,Lp,Zp,backstressprev,deqpl,dtime[k],temp[k-1],temp[k]-temp[k-1],V)
	syieldfd,dsifd,sfd,dskfd,Ycfd,Lcfd,Zcfd=mtdb(Yp,Lp,Zp,backstressprev,deqpl+toler,dtime[k],temp[k-1],temp[k]-temp[k-1],V)	
	dsi = (syieldfd-syield)/toler
        dsk = (skfd-sk)/toler
	#dsi,dsk = dsiA,dskA
	
	backstresscur=sk
	print " "
 	print "FD vs analytical: "
	print "Isotropic: ",dsi,' ',dsiA
 	print "Kinematic: ",dsk,' ',dskA
      if np.abs(f)>toler:
        print " ***  WARNING:: newton on plasticity did not converge"
        print "    increment:: ",k+1
	print "        deqpl:: ",deqpl
	print " yield stress:: ",syield
	print "  back stress:: ",sk
	print "        alpha:: ",alpha
      
      eqpl += deqpl
      alpha += (backstresscur-backstressprev)*flowdir
      #alpha = backstresscur#-backstressprev)*flowdir
      plasticcur += flowdir*deqpl
      elasticcur -= flowdir*deqpl
      stresscur = flowdir*syield+alpha
    
    elasticprev=elasticcur
    plasticprev=plasticcur
    Yp = Yc
    Lp = Lc
    Zp = Zc
    backstressprev = backstresscur
    
    stressprev = stresscur
  
    stress += [stresscur]
    strain += [logstrain]
    plastic += [plasticcur]
    elastic += [elasticcur]
    eqplhist += [eqpl]
   
  return np.array(stress),np.array(strain)
  





# ----------------------------------------------------------------
# ----------------------------------------------------------------
#       Mechanical Threshold Density Based Hardening
# ---------------------------------------------------------------- 
# ----------------------------------------------------------------
def mtdb(Yp,Cp,Zp,sykp,depl,dtime,tempp,dtemp,props):
#
  keynames = ['mu0','d0','t0','sa','sig0','a0e','C0','cD','C0m','C1','C20','a02','rate0','rate02','qe','pe','C30','C3m','C3U','C4','C50','C5U','qC1','qC2']
  [emu0,d0,t0,sa,sig0,a0e,C0,cD,C0m,C1,C20,a02,rate0,rate02,qe,pe,C30,C3m,C3U,C4,C50,C5U,qC1,qC2] = [props[kn] for kn in keynames]
  toler = 1e-12
  half=0.5
  one=1.
  two=2.
#     check 
  depl = np.abs(depl)
  if(depl<1e-10):
    depl = 1e-10
  rate = depl/dtime
  
  tempc = tempp+dtemp 
#
  if tempc>t0:
    emu = emu0 - d0/(np.exp(t0/tempc)-1.)
    emup = emu0 - d0/(np.exp(t0/tempp)-1.)
  else:
    emu = emu0
    emup = emu0
  emusf = emu/emu0
  
  sfe0 = tempc/(a0e*emu)
  dsfe = np.log(rate0/rate)*sfe0
  sfe = (one-dsfe**(one/qe))**(one/pe)
  c2expc = tempc/(emu*a02)
  c2expp = tempp/(emup*a02)
  C2c = C20*np.exp(c2expc*np.log(rate/rate02))
  C2p = C20*np.exp(c2expp*np.log(rate/rate02))
  dC2cdrt = C2c*c2expc
  dC2pdrt = C2p*c2expp

#   Yp = Xp(1)
#   Cp = Xp(2)
#   Zp = Xp(3)
#c
#c updated lattice incompatibility:
  Cc = Cp + depl*cD
#c initial guess for statistical dislocations Yc=Yp
  Yc = Yp
  sqYc = np.sqrt(Yc)
  sqYp = np.sqrt(Yp)
#c
#c static recrystalisation guess:
  C3c0 = C30*np.exp(-C3U/tempc)/rate
  C3p0 = C30*np.exp(-C3U/tempp)/rate      
  C3c = C3c0*Yc**C3m      
  C3p = C3c0*Yp**C3m
      
  dYdep = C0*Cp**C0m + C1*sqYp - C2p*Yp - C3p
  dYdec = C0*Cc**C0m + C1*sqYc - C2c*Yc - C3c
  dC3cdYc = C3c0*C3m*Yc**(C3m-one)
    
  Yinc = depl*(dYdep+dYdec)/two
  dYinc = depl*(half*C1/sqYc-C2c-dC3cdYc)/two
#c
  f = -Yinc
  dFdY = one-dYinc
      
      #write(*,*) 'F ',f
# do knewton=1,15
  counter = 0.
  while (np.abs(f)>toler)&(counter<10):
    counter+=1
    Yc = Yc - f/dFdY
    Yc = np.abs(Yc)
    sqYc = np.sqrt(Yc)
    C3c = C3c0*Yc**C3m
    dYdec = C0*Cc**C0m + C1*sqYc - C2c*Yc - C3c
    Yinc = depl*(dYdep+dYdec)/two
#c new function value
    f = Yc-Yp-Yinc
	#write(*,*) "|F|: ",np.abs(f)
	#if(np.abs(f).lt.toler) goto 20
    dC3cdYc = C3c0*C3m*Yc**(C3m-one)
    dYinc = depl*(half*C1/sqYc-C2c-dC3cdYc)/two
    dFdY = one-dYinc
  if np.abs(f)>toler:
    print "*** WARNING:: microstructure evolution did not converge"
    print " Y : ", Yc," previous value: ",Yp
#c     If Newton did not converge write output to the screen
#       write(*,*) '*** WARNING ***'
#       write(*,*) '*** THERMAL STRESS EVOLUTION DID NOT CONVERGE'
# 20    continue
#c update wall dislocation ISV
      #write(*,*) 'Y_end: ',Yc
  C5c = C50*np.exp(-C5U/tempc)/rate
  C5p = C50*np.exp(-C5U/tempp)/rate
      #Zc = (Zp+depl*C4*sqYc)/(one+depl*(C2c+C5c))
  Zc = (Zp + depl*(C4*(sqYp+sqYc)-C2p*Zp-C5p*Zp/rate)/two)/(one+ depl*(C2c + C5c/rate)/two)
  #Xc = (/Yc,Cc,Zc/) # ISV vector returned
#c
#c calculate equivalent threshold stress and gradient
  sec = sig0*sqYc
#c d(Yc)/d(epl)
      #dYdrt = -dC2cdrt*Yc + C3c/rate -dC2pdrt*Yp + C3p/rate
  dYdrt = half*(-dC2cdrt*Yc-dC2pdrt*Yp+C3c+C3p)*dtime
  dYdCc = half*depl*C0*C0m*Cc**(C0m-one)
  dYdepart = Yinc/depl
  dYcnum = dYdepart + dYdrt/dtime + dYdCc*cD
  dYc = dYcnum/dFdY    
#       dYc = (dYde-depl*dC2de*Yc+depl*C0*C0m*Cc**(C0m-one))/dFdY
#c therefore d(sec)/d(epl)
  dsecdepl = half*sig0*dYc/sqYc
#c d(sfe)/d(epl)
  dsfedepl = (sfe0*(one-dsfe**(one/qe))**(one/pe-one)*dsfe**(one/qe-one)/(pe*qe*rate))/dtime
#c Isotropic stress:
  syi = sa + emusf*sfe*sec
#c d(syi)/d(epl)
  dsyi =  emusf*(sfe*dsecdepl+dsfedepl*sec)
#c
#c Kinematic stress
  skden = one+qC2*depl
  sknum = depl*qC1*sqYc+sykp
  syk = sknum/skden
#c Gradient
  dsyk = (qC1*sqYc + half*depl*qC1*dYc/sqYc)/skden-qC2*sknum/(skden**two)
  
  
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                          FINITE DIFFERENCE PARTIAL GRADIENTS
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  eplFD = [depl+1e-10,depl,depl]
  rateFD = [rate,rate+1e-10,rate]
  CcFD = [Cc,Cc,Cc+1e-10]
  
  YcFD = []
  secFD = []
  syiFD = []
  sykFD = []
   
  for ifd in range(3):
          depl,rate,Cc = eplFD[ifd],rateFD[ifd],CcFD[ifd]
	  
# 	  sfe0 = tempc/(a0e*emu)
# 	  dsfe = np.log(rate0/rate)*sfe0
# 	  sfe = (one-dsfe**(one/qe))**(one/pe)
	  c2expc = tempc/(emu*a02)
	  c2expp = tempp/(emup*a02)
	  C2c = C20*np.exp(c2expc*np.log(rate/rate02))
	  C2p = C20*np.exp(c2expp*np.log(rate/rate02))
	  dC2cdrt = C2c*c2expc
	  dC2pdrt = C2p*c2expp

	#   Yp = Xp(1)
	#   Cp = Xp(2)
	#   Zp = Xp(3)
	#c
	#c updated lattice incompatibility:
	  #Cc = Cp + depl*cD
	#c initial guess for statistical dislocations Yc=Yp
	  YcFD += [Yp]
	  sqYc = np.sqrt(YcFD[ifd])
	  sqYp = np.sqrt(Yp)
	#c
	#c static recrystalisation guess:
	  C3c0 = C30*np.exp(-C3U/tempc)/rate
	  C3p0 = C30*np.exp(-C3U/tempp)/rate      
	  C3c = C3c0*YcFD[ifd]**C3m      
	  C3p = C3c0*Yp**C3m
      
	  dYdep = C0*Cp**C0m + C1*sqYp - C2p*Yp - C3p
	  dYdec = C0*Cc**C0m + C1*sqYc - C2c*YcFD[ifd] - C3c
	  dC3cdYc = C3c0*C3m*YcFD[ifd]**(C3m-one)
    
	  Yinc = depl*(dYdep+dYdec)/two
	  dYinc = depl*(half*C1/sqYc-C2c-dC3cdYc)/two
	#c
	  f = -Yinc
	  dFdY = one-dYinc
      
	      #write(*,*) 'F ',f
	# do knewton=1,15
	  counter = 0.
	  while (np.abs(f)>toler)&(counter<10):
	    counter+=1
	    YcFD[ifd] = YcFD[ifd] - f/dFdY
	    YcFD[ifd] = np.abs(YcFD[ifd])
	    sqYc = np.sqrt(YcFD[ifd])
	    C3c = C3c0*YcFD[ifd]**C3m
	    dYdec = C0*Cc**C0m + C1*sqYc - C2c*YcFD[ifd] - C3c
	    Yinc = depl*(dYdep+dYdec)/two
	#c new function value
	    f = YcFD[ifd]-Yp-Yinc
		#write(*,*) "|F|: ",np.abs(f)
		#if(np.abs(f).lt.toler) goto 20
	    dC3cdYc = C3c0*C3m*YcFD[ifd]**(C3m-one)
	    dYinc = depl*(half*C1/sqYc-C2c-dC3cdYc)/two
	    dFdY = one-dYinc
	  if np.abs(f)>toler:
	    print "*** WARNING:: microstructure evolution did not converge"
	    print " Y : ", YcFD[ifd]," previous value: ",Yp
	#c     If Newton did not converge write output to the screen
	#       write(*,*) '*** WARNING ***'
	#       write(*,*) '*** THERMAL STRESS EVOLUTION DID NOT CONVERGE'
	# 20    continue
	#c update wall dislocation ISV
	      #write(*,*) 'Y_end: ',Yc
	  C5c = C50*np.exp(-C5U/tempc)/rate
	  C5p = C50*np.exp(-C5U/tempp)/rate
	      #Zc = (Zp+depl*C4*sqYc)/(one+depl*(C2c+C5c))
	  Zc = (Zp + depl*(C4*(sqYp+sqYc)-C2p*Zp-C5p*Zp/rate)/two)/(one+ depl*(C2c + C5c/rate)/two)
	  #Xc = (/Yc,Cc,Zc/) # ISV vector returned
	#c
	#c calculate equivalent threshold stress and gradient
	  secFD += [sig0*sqYc]
	#c d(Yc)/d(epl)
	      #dYdrt = -dC2cdrt*Yc + C3c/rate -dC2pdrt*Yp + C3p/rate
	#c Isotropic stress:
	  syiFD += [sa + emusf*sfe*secFD[ifd]]
	#c d(syi)/d(epl)
	#c
	#c Kinematic stress
	  skden = one+qC2*depl
	  sknum = depl*qC1*sqYc+sykp
	  sykFD += [sknum/skden]
	#c Gradient
  print " "	
  print "FD partial gradients vs Analytical of Y : "
  print "d(epl):  ", (YcFD[0]-Yc)/1e-10, " vs ",dYdepart/dFdY
  print "d(rate): ", (YcFD[1]-Yc)/1e-10, " vs ",dYdrt/dFdY
  print "d(Cc):   ", (YcFD[2]-Yc)/1e-10, " vs ",dYdCc/dFdY
  #print "d(sec)/d(epl) : "
#   secFD = [] 
#   syiFD = []
#   sykFD = []
  sqYc = np.sqrt(Yc)
  dYdepartFD = (YcFD[0]-Yc)/1e-10
  dYdrtFD = (YcFD[1]-Yc)/1e-10
  dYdCcFD = (YcFD[2]-Yc)/1e-10
  dYcnumFD = dYdepartFD + dYdrtFD/dtime + dYdCcFD*cD
  dYcFD = dYcnum/dFdY
  
  dsecdeplFDpart = (secFD[0]-sec)/1e-10 + (secFD[1]-sec)/(1e-10*dtime) + cD*(secFD[2]-sec)/1e-10
  
  dsecdeplFD = half*sig0*dYcFD/sqYc
  print "d(sec)d(epl) FD, FDpartial vs Analytical: "
  print dsecdeplFD,', ',dsecdeplFDpart,', ',dsecdepl
#c d(syi)/d(epl)
  dsyiFD =  emusf*(sfe*dsecdeplFD+dsfedepl*sec)
  
  dsyiFDpart = (syiFD[0]-syi)/1e-10 + (syiFD[1]-syi)/(1e-10*dtime) + cD*(syiFD[1]-syi)/1e-10
  print "d(syi)d(epl) FD, FDpartial vs Analytical: "
  print dsyiFD,', ',dsyiFDpart,', ',dsyi
  
  return syi,dsyi,syk,dsyk,Yc,Cc,Zc














# def mtdb0(Yp,Lp,Zp,skp,depl,dtime,temp,props):
# #
#   keynames = ['mu0','d0','t0','sa','s0','si','a0e','C0','cD','r','C1','C20','a02','a0i','rate0','rate02','qe','pe','qi','pi','r0','beta','U0','Aback','Qback','mback','Cback1','Cback2','C5','C6']
#   [mu0,d0,t0,sa,s0,si,a0e,C0,cD,r,C1,C20,a02,a0i,rate0,rate02,qe,pe,qi,pi,r0,beta,U0,Aback,Qback,mback,Cback1,Cback2,C5,C6] = [props[kn] for kn in keynames]
#   toler = 1e-8
# #     check 
#   depl = np.abs(depl)
#   rate = depl/dtime
#   if(rate<1e-10):
#     rate = 1e-10
# #
#   if temp>t0:
#     mu = mu0 - d0/(np.exp(t0/temp)-1.)
#   else:
#     mu = mu0
#   musf = mu/mu0
# #  
#   sfi0 = temp/(a0i*mu)
#   lsfi = np.log(rate0/rate)*sfi0
#   sfi = (1.-lsfi**(1./qi))**(1./pi)
# #
#   sfe0 = temp/(a0e*mu)
#   lsfe = np.log(rate0/rate)*sfe0
#   sfe = (1.-lsfe**(1./qe))**(1./pe)
# #
#   c0 = -temp/(mu*a02)
#   C2 = C20*(rate/rate02)**c0
#   dC2de = C20*c0*((rate/rate02)**(c0-1.))/(dtime*rate02)
# #
#   Yc = Yp
#   Lc = Lp+depl*cD
# #
#   sqY = np.sqrt(Yc)
# #
# # static recrystalisation:
#   R0 = r0*np.exp(-U0/temp)/rate
#   Rstat = R0*Yc**beta
#   dRdY = R0*beta*Yc**(beta-1.)
# # ALTERNATIVE
# #   Rstat = r0*np.exp(-U0/temp)*np.sinh(beta*sqY/temp)/rate
# #   dRdY = r0*beta*np.exp(-U0/temp)*np.cosh(beta*sqY/temp)/(2*temp*sqY*rate)
#   #print Rstat,dRdY*depl
#   C6 = C6*np.exp(-U0/temp)/rate
# #
# # Beta_backstress
#   Betaback = Aback*np.exp(Qback/temp)**mback
# #
#   dYde = C0*Lc**r+C1*sqY-C2*Yc - Rstat
#   dFdY = 1.-depl*(0.5*C1/sqY-C2-dRdY)
#   f = Yc - Yp - depl*dYde
# #
#   counter = 0.
#   while (np.abs(f)>toler)&(counter<10):
#     counter+=1
#     f = Yc - Yp - depl*dYde
#     Yc = Yc - f/dFdY
#     sqY = np.sqrt(Yc)
# #     Rstat = r0*np.exp(-U0/temp)*np.sinh(beta*sqY/temp)/rate
# #     dRdY = r0*beta*np.exp(-U0/temp)*np.cosh(beta*sqY/temp)/(2*temp*sqY*rate)
#     Rstat = R0*Yc**beta
#     dRdY = R0*beta*Yc**(beta-1.)
#     dYde = C0*Lc**r+C1*sqY-C2*Yc-Rstat
#     dFdY = 1.-depl*(0.5*C1/sqY-C2-dRdY)
#   if np.abs(f)>toler:
#     print "*** WARNING:: microstructure evolution did not converge"
#     print " Y : ", Yc," previous value: ",Yp
# # Yield stress
#   if (sqY<0)==(sqY>0):
#     Yc=Yp
#     sqY = np.sqrt(Yc)
#   if Yc<1:
#     Yc=1.
#     sqY = 1.
#     
#   Zc = (Zp+depl*C5*sqY)/(1.+depl*C2+C6*depl)
#
#   sec = s0*sqY
#   sy = sa + musf*(sfi*si+sfe*sec)
# #
# # Partial gradient componenets
# # d(Yc)/d(epl)
#   dYc = (dYde-depl*dC2de*Yc+depl*c0*r*Lc**(r-1.))/dFdY
# # therefore d(sec)/d(epl)
#   dsecdepl = 0.5*s0*dYc/sqY
# # d(sfi)/d(epl)
#   dsfidepl = (sfi0*(1.-lsfi**(1./qi))**(1./pi-1.)*lsfi**(1./qi-1.)/(pi*qi*rate))/dtime
# # d(sfe)/d(epl)
#   dsfedepl = (sfe0*(1.-lsfe**(1./qe))**(1./pe-1.)*lsfe**(1./qe-1.)/(pe*qe*rate))/dtime
# # Total
# #d(yield)/d(epl)
#   dsy = musf*(sfe*dsecdepl+dsfedepl*sec+dsfidepl*si)
#   
# #
# # kinematic
#   if np.abs(Betaback)>toler:
#     qA = depl*Betaback
#     qB = 1.+Cback2*depl
#     qC = -(depl*Cback1*sqY+skp)
#     sk = (-qB+np.sqrt(qB**2.-4.*qA*qC))/(2.*qA) #take positive root
#   else:
#     sk = (depl*Cback1*sqY+skp)/(1.+Cback2*depl)
#     
#   #sk = skp+depl*(Cback1)#*sqY - Cback2*skp)
#   skden = (1.+Cback2*depl)
#   sknum = (depl*Cback1*sqY+skp)
#   sk = sknum/skden
#   dsk = (Cback1*sqY + 0.5*depl*Cback1*dYc/sqY)/skden - Cback2*sknum/(skden**2.)
#   #sk = (depl*Cback1*sqY+skp)/(1.+Cback2*depl)
#   #sk = skp+depl*(Cback1*sqY - Cback2*skp - Betaback*skp*np.abs(skp))
#     
#   #dsk = 0.
#   
#   return sy,dsy,sk,dsk,Yc,Lc,Zc
#
