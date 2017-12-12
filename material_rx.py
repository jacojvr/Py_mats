import numpy as np

def stresscalc(strain=np.linspace(0,0.5,101),dtime=1.,temp=25.,**kwargs):
  """ Mechanical Threshold Stress material:
  
  The function takes displ, time step sizes and temperature values into account
  to calculate the relevant stress state history
  
  INPUT:
     strain : a one dimensional array or list containing strain histories
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
     
     # Elastic Properties:
       'E_mu0'       : 100000.,
       'E_d0'        : 1000.,
       'E_t0'        : 200.,
       'E_nu'        : 0.33,
     
     # Reference stress values
       'sigma_a'     : 1.,
       'sigma_0'     : 100.,
     
     # ISV's at time=0 
       '1/Lg_0'      : 1.e-10,
       'rho_0'       : 1.,
     
     # Evolution of 1/Lg:      
       'c/grainsize' : 1000.,
       'Lm'          : 0.7, # Exponent
       'grainsize_0' : 60., # Initial grain size
     
     # Stage IV:   
       'C00'         : 100., # Constant
     
     # Storage
       'C1'          : 100.,
     
     # Dynamic recovery:
       'C20'         : 10., # Constant
       'a02'         : 0.5, # k_B / (g02 b^3) related to the activation energy
       'rate02'      : 1e10,
     
     # Thermal recovery
       'C30'         : 2., # Constant
       'C3M'         : 1., # Exponent
       'C3U'         : 5000., # Q / R related to the activation energy
     
     # Recrystallisation
       'C40'         : 1.e2, # Constant
       'C4T'         : 1.5e4, # Temperature dependence
       'C4L'         : 10, # Misorientation Constant
       'C4D'         : 0.1, # Misorientation Scaling Theta/Theta_max = min([C4D/L,1])
       'C4M'         : 4., # Misorientation exponent
       'C4A'         : 0.22, # Interfacial Area Exponent "a"
       'C4B'         : 1.3, # Interfacial Area Exponent "b"
       'C4C'         : 25., # Interfacial Area Constant "c"
       'RX_nr'       : 20, # Number of recrystalised volume fractions to track
       'RX_sizef'    : 0.8, # Recrystallised grain size fraction
     
     # Kinetic equation / scaling function
       'a0e'         : 2., # k_B / (g0e b^3) related to the activation energy
       'rate0'       : 1e7,
       'qe'          : 1.2,
       'pe'          : 0.5,
     
     # Input temperature (specify Celcius or Kelvin is assumed)
       'tempformat'  : 'Celcius',
     
     # Return additional information (Evolution of Recrystallised volume fractions, 1/Lg, rho and Approximate average grain size )
       'ADDinfo'     : False,
       
     # move recrystallised volume fractions if Xv>0.999 (Xvf2-->Xvf1-->Xvf0) 
       'moveXv'      :False,
       
     # use Finite Difference Gradients?
       'useFD':False  
       
     
  """
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##    CHECK INPUT KEYWORD ARGUMENTS
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##    Original / Default value keyword dictionary
##
  V = {
       # Elastic Properties:
       'E_mu0':100000.,
       'E_d0':1000.,
       'E_t0':200.,
       'E_nu':0.33,
       # Reference stress values
       'sigma_a':1.,
       'sigma_0':100.,
       # ISV's at time=0 
       '1/Lg_0':1.e-10,
       'rho_0':1.,
       # Evolution of 1/Lg:      
       'c/grainsize':1000.,
       'Lm':0.7, # Exponent
       'grainsize_0':60., # Initial grain size
       # Stage IV:   
       'C00':100., # Constant
       # Storage
       'C1':100.,
       # Dynamic recovery:
       'C20':10., # Constant
       'a02':0.5, # \frac{k_B}{g02 b^3} related to the activation energy
       'rate02':1e10,
       # Thermal recovery
       'C30':2., # Constant
       'C3M':1., # Exponent
       'C3U':5000., # \frac{Q/R} related to the activation energy
       # Recrystallisation
       'C40':1.e2, # Constant
       'C4T':1.5e4, # Temperature dependence
       'C4L':10., # Misorientation Constant
       'C4D': 0.1, # Misorientation Scaling Theta/Theta_max = min([C4D/L,1])
       'C4M':4., # Misorientation exponent
       'C4A':0.22, # Interfacial Area Exponent "a"
       'C4B':1.3, # Interfacial Area Exponent "b"
       'C4C':25., # Interfacial Area Constant "c"
       'RX_nr':20, # Number of recrystalised volume fractions to track
       'RX_sizef':0.8, # Recrystallised grain size fraction
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
       'moveXv':True,
       # use Finite Difference Gradients?
       'useFD':False
       }
##
##    Check input keywords
##
  givedef = 0.
  if kwargs.__len__()>0:
    for vardef in kwargs.keys():
      if vardef in V.keys():
	V[vardef] = kwargs[vardef]
      else:
        givedef = 1
        print "\n\t*** WARNING *** keyword - "+vardef+" - unknown."
  if givedef:
    print "\n ACCEPTED LIST OF KEYWORDS:\n"
    for vardef in V.keys():
      print "\t"+vardef+"  -  with default value: ",V[vardef]
##
##    Make sure inputs are consistent:
##
  strain = np.array(strain,float)
  dtime = np.array(dtime,float)
  temp = np.array(temp,float)
  nrs = strain.size
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
##
##    Remove possible zero strain at the start:
##
  stress = []
  plastic = []
  elastic = []
  eqplhist = []
  strainprev=0.
  elasticprev=0.
  plasticprev=0.
  eqpl = 0.
  stressprev = 0.
  if strain[0]==0:
    nrs -=1
    stress = [0.]
    plastic = [0.]
    elastic =[0.]
    eqplhist =[0.]
    strain = strain[1::]
    dtime = dtime[1::]
    temp = temp[1::]
    if V['tempformat']=='Celcius':
      temp+=273.
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##    INITIALISE SOLUTION AND ARRAYS CONTAINING HISTORIES
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##    Constants
##
  c1,c2 = 2./3, np.sqrt(2./3)
  toler = 1e-8
##    Number of recrystallised volumes to track
##
  RX_nr = V['RX_nr']
  moveRX = V['moveXv']
##
##    ISV's:
##
  Yp = V['rho_0'] 
  Y_plst = np.ones(RX_nr)*Yp
  Y_hist = [Y_plst]
  Lp = V['1/Lg_0']
  L_plst = np.ones(RX_nr)*Lp
  L_hist = [L_plst]
##
##    Estimate final crystal size in each volume fraction
##
  D_plst = V['grainsize_0']*(V['RX_sizef']**np.array(range(RX_nr)))
  D_eq = [V['grainsize_0']]
##
##    Volume fractions
##
  Xv_plst = (np.ones(RX_nr+1)*toler)
  Xv_plst[0] = 1.
  Xv_hist = [Xv_plst]
##
##    Equivalent plastic strain per Volume fraction
##
  P_plst = np.zeros(RX_nr)
  P_hist = [P_plst]
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     START THE SOLUTION
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##    Elastic Properties
##
  [mu0,d0,t0,nu] = [V[kn] for kn in ['E_mu0','E_d0','E_t0','E_nu']]
##
##    Solution Info
##
  totalnewton = 0
  totalplastic = 0.
##
##    Loop over input strain array
##
  for k in range(nrs):
    dstrain = strain[k]-strainprev
    strainprev = strain[k]
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
    smises = np.abs(stresscur)    
    flowdir = (stresscur)/smises
        
    syield0,dsi,X2lst,Y_clst,L_clst,P_clst,eqD,Xpl = mtdbrec(Xv_plst,Y_plst,L_plst,P_plst,D_plst,toler,dtime[k],temp[k],V)    
    X2lst,Y_clst,L_clst,P_clst = Xv_plst.copy(),Y_plst.copy(),L_plst.copy(),P_plst.copy()
    if V['useFD']:
      syield0FD,dsi,X2lstFD,Y_clstFD,L_clstFD,P_clstFD,eqDFD,XplFD = mtdbrec(Xv_plst,Y_plst,L_plst,P_plst,D_plst,toler*2,dtime[k],temp[k],V) 
      dsi = (syield0FD-syield0)/toler
      
    f = smises - syield0
    
    if f>toler:
      totalplastic += 1
      syield = syield0
      counter = 0.
      while (np.abs(f)>toler)&(counter<10):
        counter+=1
	totalnewton+=1
        f = smises-emod*deqpl-syield
        deqpl=deqpl+f/(emod+dsi)
	deqpl = np.abs(deqpl)
	syield,dsi,X2lst,Y_clst,L_clst,P_clst,eqD,Xpl=mtdbrec(Xv_plst,Y_plst,L_plst,P_plst,D_plst,deqpl,dtime[k],temp[k],V)
        
	if V['useFD']:
          syieldFD,dsiFD,X2lstFD,Y_clstFD,L_clstFD,P_clstFD,eqDFD,XplFD = mtdbrec(Xv_plst,Y_plst,L_plst,P_plst,D_plst,deqpl+toler,dtime[k],temp[k],V) 
	  dsiFD = (syieldFD-syield)/toler
	  print "FD: ",dsiFD," vs A: ",dsi
	  dsi = dsiFD 
	  
      if np.abs(f)>toler:
        print " ***  WARNING:: newton on plasticity did not converge"
        print "    increment:: ",k+1
	print "        deqpl:: ",deqpl
	print " yield stress:: ",syield
	print " Volume Fract:: ",X2lst[np.where(X2lst>toler)[0]]
      eqpl += deqpl
      plasticcur += flowdir*deqpl
      elasticcur -= flowdir*deqpl
      stresscur = flowdir*syield
    
    elasticprev=elasticcur
    plasticprev=plasticcur
    Xv_plst,Y_plst,L_plst,P_plst = X2lst.copy(),Y_clst.copy(),L_clst.copy(),P_clst.copy()
    
#check and move volumes:
    if (Xv_plst[1]>0.999)&moveRX:
      Xv_plst[0:RX_nr] = Xv_plst[1:]
      Y_plst[0:RX_nr-1] = Y_plst[1:]
      L_plst[0:RX_nr-1] = L_plst[1:]
      D_plst[0:RX_nr-1] = D_plst[1:]
      P_plst[0:RX_nr-1] = P_plst[1:]
      Xv_plst[0] = 1.
      Xv_plst[-1] = toler#*Xv_plst[-2]
      Y_plst[-1] = 1.
      L_plst[-1] = 0.
      D_plst[-1] = V['RX_sizef']*D_plst[-2]
      P_plst[-1] = 0.
    
    
    
    
    Xv_hist += [Xv_plst]
    Y_hist += [Y_plst]
    L_hist += [L_plst]
    P_hist += [P_plst]
    
    stress += [stresscur]
    plastic += [plasticcur]
    elastic += [elasticcur]
    eqplhist += [Xpl]
    D_eq += [eqD]
   
  print "Total Newton Loops: ",totalnewton
  print "Average per plastic increment: ",totalnewton/totalplastic
  
  if V['ADDinfo']:
    return np.array(stress),np.array(Xv_hist),np.array(Y_hist),np.array(L_hist),np.array(P_hist),np.array(D_eq),np.array(eqplhist)
  
  return np.array(stress)
  


#
# ----------------------------------------------------------------
# ----------------------------------------------------------------
#    Density Based Hardening, Recovery and Recrystallisation
# ---------------------------------------------------------------- 
# ----------------------------------------------------------------
def getF(Xvec,Xv_p,Yp,Lp,X1,X1R,dX1dE,D1,depl,dtime,temp,props):
#
  keynames = ['E_mu0','E_d0','E_t0','sigma_a','sigma_0','a0e','C00','c/grainsize','Lm','C1','C20','a02','rate0','rate02','qe','pe','C30','C3M','C3U','C40','C4T','C4L','C4D','C4M','C4C','C4A','C4B']
  [mu0,d0,t0,sa,s0,a0e,C0,cD0,r,C1,C20,a02,rate0,rate02,qe,pe,r0,beta,U0,Rx0,Rx1,Qx0,QxD,Qx1,cx,ax,bx] = [props[kn] for kn in keynames]
  toler = 1e-4
#     check 
  depl = np.abs(depl)
  rate = depl/dtime
  if(rate<1e-10):
    rate = 1e-10
#
  if temp>t0:
    mu = mu0 - d0/(np.exp(t0/temp)-1.)
    Rrx = Rx0*np.exp(-Rx1/temp)*mu
  else:
    mu = mu0
    Rrx = Rx0*mu
  musf = mu/mu0
#
  cD = cD0/D1
  Y1,L1 = np.max([np.abs(Xvec[0]),toler]),np.max([np.abs(Xvec[1]),toler])
  sqY = np.sqrt(Y1)
  # determine X2 associated with Y1 and L1:
  misalign = np.min([QxD*L1,1.])
  Qrx = (1.-np.exp(-Qx0*(misalign)**Qx1))
  #
  QRYx = Rrx*Qrx*Y1
  ##
  # guess X1 = previous value
  X2 = Xv_p
  GRX = X1*((X2/X1)**ax)*((1.-X2/X1)**bx)*(1+cx*(1.-X1))
  dGRXdX2 = ax*X1*((X2/X1)**(ax-1.))*((1.-X2/X1)**bx)*(1+cx*(1.-X1)) - bx*X1*((X2/X1)**ax)*((1.-X2/X1)**(bx-1.))*(1+cx*(1.-X1))
  
  #GRX = X1*(X1-X2)
  #dGRXdX2 = -X1
  X2R = np.abs(QRYx*GRX)
  fxc = X2 - Xv_p - dtime*X2R
  counter=0.
  while (np.abs(fxc)>toler)&(counter<15):
    counter+=1
    dfxc = 1.- dtime*QRYx*dGRXdX2
    #print dfxc
    X2 = X2 - fxc/dfxc
    X2 = np.abs(X2)
    GRX = X1*((X2/X1)**ax)*((1.-X2/X1)**bx)*(1+cx*(1.-X1))
    dGRXdX2 = ax*((X2/X1)**(ax-1.))*((1.-X2/X1)**bx)*(1+cx*(1.-X1)) - bx*((X2/X1)**ax)*((1.-X2/X1)**(bx-1.))*(1+cx*(1.-X1))
    #GRX = X1*(X1-X2)
    #dGRXdX2 = -X1
    if X2>X1:
      print "fully recrystallise current volume fraction"
      F = np.array([Y1-1.,L1])
      dFdX = np.eye(2)
      X2vec = np.array([toler,0,0,0])
      GF = np.array([0,0])
      return F,dFdX,X2vec,GF
    X2R = np.abs(QRYx*GRX)
    fxc = X2 - Xv_p - dtime*X2R
  if np.abs(fxc)>toler:
    print "*** WARNING:: Additional Recrystallised volume fraction did not converge"
    print "              attempting to use a volume fraction of = ",X2
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # PARTIAL : CHANGE OF X2 WITH RESPECT TO Y1 AND L1:
  # partial gradients d(X2)/d(Y1)
  
  dX2dY1 = (dtime*Rrx*Qrx*GRX)/(1.-dtime*Rrx*Qrx*Y1*dGRXdX2)
  # partial gradients d(X2)/d(L1)
  dMdL1 = 0.
  if misalign<1:
    dMdL1 = QxD
  dQRXdM1 = Qx0*Qx1*np.exp(-Qx0*(misalign)**Qx1)*(misalign)**(Qx1-1.)
  dX2dM1 = (dtime*Rrx*dQRXdM1*Y1*GRX)/(1.-dtime*Rrx*Qrx*Y1*dGRXdX2)
  # partial gradients d(X2)/d(L1)
  dX2dL1 = dX2dM1*dMdL1
  
  
#   print "d(X2)/d(Y1) FD= ",dX2dYFD," vs analytical ",dX2dY1
#   print "d(X2)/d(L1) FD= ",dX2dLFD," vs analytical ",dX2dL1
#   print "X1: ",X1
#   print "X2: ",X1
#   print "dGRXdX2: ",dGRXdX2
#   print "dX2dY1: ",dX2dY1
#   print "dX2dL1: ",dX2dL1
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  c0 = -temp/(mu*a02)
  C2 = C20*(rate/rate02)**c0
  dC2de = C20*c0*((rate/rate0)**(c0-1.))/(dtime*rate0)
  R0 = r0*np.exp(-U0/temp)/rate
  Rstat = R0*Y1**beta
  Hy = C0*L1+C1*sqY-C2*Y1 - Rstat
  # Function values:   
  Xrec0 = 1./(X1 - X2)
  Xrec = dtime*X1R*Xrec0
  dXrecdX = Xrec*Xrec0
  #
  R1 = Y1 - Yp - Hy*depl + Y1*Xrec
  R2 = L1-Lp-depl*cD*(L1)**(1.-1./r) + L1*Xrec
  F = np.array([R1,R2])
  
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # PARTIAL : CHANGE IN RESIDUAL WITH RESPECT TO Y1 AND L1:
  dHydY1 = 0.5*C1/sqY - C2 - beta*R0*Y1**(beta-1.)
  # partial gradients d(F1)/d(Y1)
  dR1dY0 = 1.-dHydY1*depl + Xrec
  dR1dY = dR1dY0 + Y1*dXrecdX*dX2dY1
  # partial gradients d(F1)/d(L1)
  dR1dL0 = -C0*depl
  dR1dL = dR1dL0 + Y1*dXrecdX*dX2dL1
  # partial gradients d(F2)/d(Y1)
  dR2dY = +L1*dXrecdX*dX2dY1
  # partial gradients d(F2)/d(L1)
  dR2dL0 =  1.- (1.-1./r)*depl*cD*(L1)**(-1./r) + Xrec
  dR2dL = dR2dL0 + L1*dXrecdX*dX2dL1
  #
  dFdX = np.array([[dR1dY,dR1dL],[dR2dY,dR2dL]])
  #
  dR1dE = Hy - depl*dC2de*Y1 + Rstat
  dR2dE = cD*(L1)**(1.-1./r)
  
  GF = np.array([dR1dE,dR2dE])
  #[dY1dE,dL1dE] = np.linalg.solve(dFdX,GF)
  #dX2dE = dX2dY1*dY1dE + dX2dL1*dL1dE
  
  #gradvec = np.array([dY1dE,dX2dE,dX2dY1,dX2dL1,dL1dE])
  
  X2vec = np.array([X2,X2R,dX2dY1,dX2dL1])
  
  return F,dFdX,X2vec,GF
 




def mtdbrec(Xv_plst,Y_plst,L_plst,P_plst,D_plst,depl,dtime,temp,props):
#
  keynames = ['E_mu0','E_d0','E_t0','sigma_a','sigma_0','a0e','C00','c/grainsize','Lm','C1','C20','a02','rate0','rate02','qe','pe','C30','C3M','C3U','C40','C4T','C4L','C4D','C4M','C4C','C4A','C4B']
  [mu0,d0,t0,sa,s0,a0e,C0,cD,r,C1,C20,a02,rate0,rate02,qe,pe,r0,beta,U0,Rx0,Rx1,Qx0,QxD,Qx1,cx,ax,bx] = [props[kn] for kn in keynames]
  toler = 1e-4
#     check 
  depl = np.abs(depl)
  rate = depl/dtime
  if(rate<1e-10):
    rate = 1e-10
#
  if temp>t0:
    mu = mu0 - d0/(np.exp(t0/temp)-1.)
    Rrx = Rx0*np.exp(-Rx1/temp)*mu
  else:
    mu = mu0
    Rrx = Rx0*mu
  musf = mu/mu0
#
  sfe0 = temp/(a0e*mu)
  lsfe = np.log(rate0/rate)*sfe0
  sfe = (1.-lsfe**(1./qe))**(1./pe)
#
  c0 = -temp/(mu*a02)
  C2 = C20*(rate/rate02)**c0
  dC2de = C20*c0*((rate/rate0)**(c0-1.))/(dtime*rate0)
#
  RX_nr = Xv_plst.size-1   
#
  X2lst,Y_clst,L_clst,P_clst = Xv_plst.copy(),Y_plst.copy(),L_plst.copy(),P_plst.copy()
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ISV EVOLUTION AND RECRYSTALLISING VOLUME FRACTION:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  eqD = 0.
  Xpl = 0.
  YEQ = 0.
  dYEQdE = 0
  dX1dE = 0
  X1 = 1.
  X1R = 0.
  
  iXvf = 0
  checknextXv = True
  while (iXvf<RX_nr)&(checknextXv):
    P_clst[iXvf] += depl
    D1 = D_plst[iXvf]
    X2_p = Xv_plst[iXvf+1]      
    Y1_p = Y_plst[iXvf]
    L1_p = L_plst[iXvf]
    Xvec = np.array([Y1_p,L1_p])
    F,dFdX,X2vec,GF = getF(Xvec,X2_p,Y1_p,L1_p,X1,X1R,dX1dE,D1,depl,dtime,temp,props)
    #print "XVEC 0: ", Xvec
    fv = np.sqrt(np.sum(F*F))
    detFmat = (dFdX[0,0]*dFdX[1,1]-dFdX[1,0]*dFdX[0,1])
    counter = 0
    while (np.abs(fv)>toler)&(counter<20):
      counter+=1
      if np.abs(detFmat)>0:
        Xvec = Xvec - np.linalg.solve(dFdX,F)
        if Xvec[0]==1.:
	  # ISV's have been reset completely -- no further checks on recrystallisation
	  checknextXv = False
	else:
          Xvec = np.abs(Xvec)
          F,dFdX,X2vec,GF = getF(Xvec,X2_p,Y1_p,L1_p,X1,X1R,dX1dE,D1,depl,dtime,temp,props)
          fv = np.sqrt(np.sum(F*F))
          detFmat = (dFdX[0,0]*dFdX[1,1]-dFdX[1,0]*dFdX[0,1])
      else:
        print "*** WARNING:: SINGULAR MATRIX"
	print "              Values returned unchanged"
        Xvec = Xvec # unchanged ISV's
        fv=0.
    
    if np.abs(fv)>toler:
      print "*** WARNING:: microstructure evolution did not converge"
      print "Residuals: R(Y) = ",F[0]," R(L) = ",F[1]
      
    [Y1,L1] = Xvec#[0],Xvec[1]
    [X2,X2R,dX2dY1,dX2dL1] = X2vec
    
    Y_clst[iXvf]=np.max([Y1,toler])
    L_clst[iXvf]=L1
    X2lst[iXvf+1]=X2 #np.max([X2,toler]) 
    
    if X2<0.001:
      checknextXv = False
      X2 = 0. # "next recrystallised volume does not contribute to average dislaocation density ratio 
  
    YEQ += Y1*(X1-X2)
    # equivalent grain size
    eqD += D1*(X1-X2)
    # equivalent plastic strain
    Xpl += P_clst[iXvf]*(X1-X2)
    
    if np.abs(detFmat)>toler:
      [dY1dE,dL1dE] = np.linalg.solve(dFdX,GF)
      dX2dE = dX2dY1*dY1dE + dX2dL1*dL1dE
      dYEQdE += dY1dE*(X1-X2) + Y1*(dX1dE-dX2dE)
# update X2 --> X1 and redo for the next volume fraction:
      dX1dE = dX2dE
      X1 = X2
      if X1>0.999:
        X1R = 0.
      else:
        X1R = X2R
    else:
      checknextXv = False

    iXvf+=1
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              END
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sqYEQ = np.sqrt(YEQ)
#
  sec = s0*sqYEQ
  sy = sa + musf*(sfe*sec)
#
# Partial gradient componenets
# therefore d(sec)/d(epl)
  dsecdepl = 0.5*s0*dYEQdE/sqYEQ
# d(sfe)/d(epl)
  dsfedepl = (sfe0*(1.-lsfe**(1./qe))**(1./pe-1.)*lsfe**(1./qe-1.)/(pe*qe*rate))/dtime
# Total
  dsy = musf*(sfe*dsecdepl+dsfedepl*sec)
  return sy,dsy,X2lst,Y_clst,L_clst,P_clst,eqD,Xpl

