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
       's0e' = 0.
      's0es' = 450.
        'si' = 0.
       'a0e' = 2.
      'a0es' = 0.55
       'a0i' = 0.9
     'hard0' = 2500.
     'rate0' = 1e7
   'rate0es' = 1e10
        'qe' = 1.2
        'pe' = 0.5
        'qi' = 1.5
        'pi' = 2./3
     'alpha' = 2.
     'model' = 'voce'
'tempformat' = 'Celcius'
       
     
  """
  # Dictionary containing variables and default values        
  V = {'mu0':85000.,
       'd0':8300.,
       't0':200.,
       'nu':0.33,
       'sa':35.,
       's0e':0.,
       's0es':450.,
       'si':0.,
       'a0e':2.,
       'a0es':0.55,
       'a0i':0.9,
       'hard0':2500.,
       'rate0':1e7,
       'rate0es':1e10,
       'qe':1.2,
       'pe':0.5,
       'qi':1.5,
       'pi':2./3,
       'alpha':2.,
       'model':'tanh',
       'tempformat':'Celcius'}
	  
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
  sep = V['s0e']
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
    #flowdir = stresscur/np.abs(stresscur)
    # Absolute value
    smises = np.abs(stresscur)
    flowdir = (stresscur)/smises
        
    syield0,dsi,sec=mtshard(sep,toler,dtime[k],temp[k],V)
    sec = sep
    
    
    f = smises - syield0
    
    if f>toler:
      syield = syield0
      counter = 0.
      while (np.abs(f)>toler)&(counter<10):
        f = smises-emod*deqpl-syield
        deqpl=deqpl+f/(emod+dsi)
	syield,dsi,sec=mtshard(sep,deqpl,dtime[k],temp[k],V)
      if np.abs(f)>toler:
        print "*** WARNING:: newton on plasticity did not converge"
      
      eqpl += deqpl
      plasticcur += flowdir*deqpl
      elasticcur -= flowdir*deqpl
      stresscur = flowdir*syield
    
    elasticprev=elasticcur
    plasticprev=plasticcur
    sep = sec
  
    stress += [stresscur]
    strain += [logstrain]
    plastic += [plasticcur]
    elastic += [elasticcur]
    eqplhist += [eqpl]
   
  return np.array(stress),np.array(strain)
  





# ----------------------------------------------------------------
# ----------------------------------------------------------------
#       Mechanical Threshold Stress Hardening
# ---------------------------------------------------------------- 
# ----------------------------------------------------------------
def mtshard(sep,depl,dtime,temp,props):
#
  keynames = ['mu0','d0','t0','sa','s0es','si','a0e','a0es','a0i','hard0','rate0','rate0es','qe','pe','qi','pi','alpha']
  [mu0,d0,t0,sa,s0es,si,a0e,a0es,a0i,hard0,rate0,rate0es,qe,pe,qi,pi,alpha] = [props[kn] for kn in keynames]
  toler = 1e-8
  cmname = props['model']
#     check 
  depl = np.abs(depl)
  rate = depl/dtime
  if(rate<1e-10):
    rate = 1e-10
#
  if temp>t0:
    mu = mu0 - d0/(np.exp(t0/temp)-1.)
  else:
    mu = mu0
  musf = mu/mu0
#  
  sfi0 = temp/(a0i*mu)
  lsfi = np.log(rate0/rate)*sfi0
  sfi = (1.-lsfi**(1./qi))**(1./pi)
#
  sfe0 = temp/(a0e*mu)
  lsfe = np.log(rate0/rate)*sfe0
  sfe = (1.-lsfe**(1./qe))**(1./pe)
#
  c0 = temp/(mu*a0es)
  sat = s0es*(rate/rate0es)**c0
#
  sec = sep
#
  if cmname.lower() == 'tanh':
    hard = hard0*(1.-(np.tanh(alpha*sec/sat)/np.tanh(alpha)))
    dhard = depl*hard0*alpha/(np.tanh(alpha)*sat*np.cosh(alpha*sec/sat)**2.)
  else:
    hard = hard0*(1.-(sec/sat))**alpha
    dhard = (depl*hard0*alpha*(1.-(sec/sat))**(alpha-1.))/sat
  
  f = sec - sep - depl*hard
#
  if np.abs(f)>toler:
    counter = 0.
    while (np.abs(f)>toler)&(counter<10):
       f = sec - sep - depl*hard
       sec = sec - f/(1.-dhard)
       if cmname.lower() == 'tanh':
         hard = hard0*(1.-(np.tanh(alpha*sec/sat)/np.tanh(alpha)))
	 dhard = depl*hard0*alpha/(np.tanh(alpha)*sat*np.cosh(alpha*sec/sat)**2.)
       else:
         hard = hard0*(1.-(sec/sat))**alpha
	 dhard = (depl*hard0*alpha*(1.-(sec/sat))**(alpha-1.))/sat
  if np.abs(f)>toler:
    print "*** WARNING:: microstructure evolution did not converge"
# Yield stress
  sy = sa + musf*(sfi*si+sfe*sec)
#
# Partial gradient componenets
# d(sec)/d(epl)
  dsecdepl = hard/(1.+dhard)  
# d(sfi)/d(epl)
  dsfidepl = sfi0*(1.-lsfi**(1./qi))**(1./pi-1.)*lsfi**(1./qi-1.)/(pi*qi*rate)
# d(sfe)/d(epl)
  dsfedepl = sfe0*(1.-lsfe**(1./qe))**(1./pe-1.)*lsfe**(1./qe-1.)/(pe*qe*rate)
# d(sec)/d(rate)
  dsecdrt = c0*dhard*s0es*sec*(rate/rate0es)**(c0-1.)/((dhard+1.)*rate0es*sat)
# Total
#d(yield)/d(epl)
  dsy =  musf*(sfe*(dsecdepl+dsecdrt/dtime)+dsfedepl*sec/dtime+dsfidepl*si/dtime)
  return sy,dsy,sec

