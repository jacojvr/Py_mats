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
     
      'emod' = 120000.
        'sa' = 90.
         'C' = 0.
        'C1' = 131.1
       'C20' = 52.35
        'C3' = 0.
        'C4' = 0.
        'C5' = 0.
        'C6' = 0.
        'C7' = 0.
        'C8' = 7.39e4
        'C9' = 1.14e3
         'm' = 5.97
         'n' = 13.25
         'q' = 0.
	'xi' = 1.
     'rate0' = 1.
      'sig0' = 325.3
 'integrate' = 'explicit' alternatively 'implicit'
'tempformat' = 'Kelvin' alternatively 'Celcius'
'timeformat' = 'Total' alternatively 'Incremental'
       
     
  """
  # Dictionary containing variables and default values        
  V = {'emod':120000.,
       'sa':90.,
       'C':0.,
       'C1':131.1,
       'C20':52.35,
       'C3':0.,
       'C4':0.,
       'C5':0.,
       'C6':0.,
       'C7':0.,
       'C8':7.39e4,
       'C9':1.14e3,
       'm':5.97,
       'n':13.25,
       'q':0.,
       'xi':1.,
       'rate0':1.,
       'sig0':325.3,
       'integrate':'implicit',
       'tempformat':'Kelvin',
       'timeformat':'Total'}
	  
  toler = 1e-8
  # Check key word arguments defined during function call and set values in variable dictionary
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
      print "\t"+vardef+" with default value: ",V[vardef]
	
  # Initialise equivalent plastic strain and kinematic plasticity
  eqpl = 0.
  sykp = 0. 
  V['C7star']=V['C7']
  timeformat = V['timeformat']
  Xep,Yep,Zep = 1., 1., 0. # Evolving ISV's of dislocation density ratios
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
      temp +=273.
    
  # Constants
  c1,c2 = 2./3, np.sqrt(2./3)
  # Begin solution:
  # Get variables:
  emod = V['emod']
  sig0 = V['sig0']
  C7 = V['C7']
  for k in range(nrs):
    #print "*** INCREMENT ",k+1
    logstrain = np.log(1.+displ[k])
    dstrain = logstrain-logstrainprev
    logstrainprev = logstrain
    # Elastic guess for strain and stress increment:
    elasticcur = elasticprev+dstrain
    deqpl = 0.
    plasticcur = plasticprev
    stresscur = emod*elasticcur
    #flowdir = stresscur/np.abs(stresscur)
    # Absolute value
    smises = np.abs(stresscur-alpha)
    flowdir = (stresscur-alpha)/smises
    syield0,dsi,backstresscur,dsk,Xec,Yec,Zec = estrin(Xep,Yep,Zep,backstressprev,toler,dtime[k],temp[k],V)
    backstresscur = backstressprev
    Xec,Yec,Zec = Xep,Yep,Zep
    
    f = smises - syield0
    
    if (syield0<sig0*np.sqrt(Yec)):
      V['C7star']=C7
    else:
      V['C7star']=0.
    
    if f>toler:
      deqpl=toler
      syield = syield0
      counter = 0.
      while (np.abs(f)>toler)&(counter<10):
        counter+=1
        f = smises-emod*deqpl-syield-backstresscur+backstressprev
        deqpl=deqpl+f/(emod+dsi+dsk)
	syield,dsi,backstresscur,dsk,Xec,Yec,Zec = estrin(Xep,Yep,Zep,backstressprev,deqpl*flowdir,dtime[k],temp[k],V)
# 	syieldfd,dsifd,backstresscurfd,dskfd,Xecfd,Yecfd,Zecfd = estrin(Xep,Yep,Zep,backstressprev,(deqpl+toler),dtime[k],temp[k],V)
# 	dsi = (syieldfd-syield)/toler
# 	dsk = (backstresscurfd-backstresscur)/toler
	if (syield<sig0*np.sqrt(Yec)):
          V['C7star']=C7
        else:
          V['C7star']=0.
        #print "ABS(F):: ", np.abs(f)
      if np.abs(f)>toler:
        print "*** WARNING:: newton on plasticity did not converge"
      
      eqpl += deqpl
      alpha += (backstresscur-backstressprev)*flowdir
      #alpha = backstresscur*flowdir
      plasticcur += flowdir*deqpl
      elasticcur -= flowdir*deqpl
      stresscur = flowdir*syield + alpha
    
    if stresscur<toler:
      Yec = Yec - Zec
      Zec = 0.
    elif stresscur/np.abs(stresscur)==-(stressprev)/np.abs(stressprev): #stress reversal
      # inject wall dislocations back into forest 
      Yec = Yec - Zec
      Zec = 0.
    
    print "dislocation evolution: ", Xec,Yec,Zec
    #print "strain stress : ", logstrain, stresscur
      
    Xep,Yep,Zep = Xec,Yec,Zec
    elasticprev=elasticcur
    plasticprev=plasticcur
    alphaprev = alpha
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
# dislocation density evolution form 
# "UNIFIED CONSTITUTIVE LAWS OF PLASTIC DEFORMATION"
# by  Yuri Estrin
# ---------------------------------------------------------------- 
# ----------------------------------------------------------------
def estrin(Xep,Yep,Zep,sykp,depl,dtime,temp,props):
#
  keynames = ['sa','C','C1','C20','C3','C4','C5','C6','C7star','C8','C9','m','n','q','xi','rate0','sig0']
  [sa,C,C1,C20,C3,C4,C5,C6,C7star,C8,C9,m,n,q,xi,rate0,sig0] = [props[kn] for kn in keynames]
  toler = 1e-8
#     plastic strain rate
  depl=np.abs(depl)
  rate = np.abs(depl)/dtime
  if(rate<1e-10):
    rate = 1e-10
#
  c00 = -1./n
  C2 = C20*(rate/rate0)**c00
  dC2dep = C20*c00*((rate/rate0)**(c00-1.))/(dtime*rate0)
#
# ##   FORWARD INTEGRATION
  if props['integrate']=='explicit':
    dXdep = q*(-C - C1*  np.sqrt(Yep) - C3*Xep+ C4*Yep/Xep)
    dYdep = C + C1*  np.sqrt(Yep) - C2*Yep -depl*dC2dep*Yep + C3*Xep
    Xec = Xep + depl*dXdep
    Yec = Yep + depl*dYdep
    Zec = Zep + depl*(C5*  np.sqrt(Yep) - C2*Zep - C6*Zep/rate) 
  else:
    X = np.array([Xep,Yep,Zep])
#
    F = np.array([- depl*q*(-C - C1*  np.sqrt(X[1]) - C3*X[0] + C4*X[1]/X[0]+C7star*X[2]),
                - depl*(C + C1*  np.sqrt(X[1]) - C2*X[1] + C3*X[0]),
                - depl*(C5*  np.sqrt(X[1]) - C2*X[2] - C6*X[2]/rate) ])
#
    FF =   np.sqrt(F[0]*F[0]+ F[1]*F[1] + F[2]*F[2])
# 
    counter=0
    while (FF>toler)&(counter<10):
      counter+=1    
      A = np.eye(3)
      A[0,0]-=depl*q*(-C3-C4*X[1]/(X[0]*X[0]))
      A[0,1]=-depl*q*(- 0.5*C1/  np.sqrt(X[1]) + C4/X[0])
      A[0,2]=-depl*q*C7star
      A[1,0]=-depl*C3
      A[1,1]-=depl*(0.5*C1/np.sqrt(X[1]) - C2) 
      A[2,1]=-depl*0.5*C5/  np.sqrt(X[1])
      A[2,2]+=depl*(C2 + C6/rate)
      #print X,A
      X = X - np.linalg.solve(A,F)
      F = np.array([ X[0] - Xep - depl*q*(-C - C1*np.sqrt(X[1]) - C3*X[0] + C4*X[1]/X[0]+C7star*X[2]),
                   X[1] - Yep - depl*(C + C1*np.sqrt(X[1]) - C2*X[1] + C3*X[0]),
                   X[2] - Zep - depl*(C5*np.sqrt(X[1]) - C2*X[2] - C6*X[2]/rate) ])
      FF =   np.sqrt(F[0]*F[0]+ F[1]*F[1] + F[2]*F[2])
    if FF>toler:
      print "*** WARNING:: microstructure evolution did not converge"
    Xec = X[0]
    Yec = X[1]
    Zec = X[2]
#
    dXdep = q*(-C - C1*  np.sqrt(Yec) - C3*Xec+ C4*Yec/Xec)
    dYdep = C + C1*  np.sqrt(Yec) - C2*Yec -depl*dC2dep*Yec + C3*Xec
#    dZdep
#
#     Isotropic hardening
  syi = sa + sig0*((rate/xi)**(1./m))*(Xec**(-1./m))*(Yec**(0.5))
  #print 'Yield stress: ',syi
#
#     isotropic yield gradient components
  dsyA = ((rate/xi)**(1./m - 1.))/(m*xi*dtime)
  dsyB = (-1./m)*(Xec**(-1./m-1.))*dXdep
  dsyC = 0.5*(Yec**(-0.5))*dYdep
#     and gradient
  dsyi = sig0*(dsyA*(Xec**(-1./m))*(Yec**(0.5))+
         ((rate/xi)**(1./m))*dsyB*(Yec**(0.5))+
         ((rate/xi)**(1./m))*(Xec**(-1./m))*dsyC)
#
#     Kinematic Hardening
  syk = (1.-C9*depl)*sykp + C8*np.sqrt(Yec)*depl
#     Gradient:
  dsyk = C8*(np.sqrt(Yec)+depl*0.5*dYdep/np.sqrt(Yec)) - C9*sykp
  
  #print 'Back stress: ',syk

  return syi,dsyi,syk,dsyk,Xec,Yec,Zec
