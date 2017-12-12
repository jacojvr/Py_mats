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
     ** key word arguments:
       
     
  """
  # Dictionary containing variables and default values  
  V = { 'emod':1000., # Young's modulus
        'sa':50., # Initial yield stress
        'dsi':0., # Linear sotropic hardening modulus
	'dsk':250., # Linear kinematic hardening modulkus
	}
	  
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
    temp = temp[1::]+273.
    
  # Constants
  c1,c2 = 2./3, np.sqrt(2./3)
  # Begin solution:
  # Get variables:
  emod = V['emod']
  sa = V['sa']
  dsi = V['dsi']
  dsk = V['dsk']
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
        
    syield0,backstresscurr = combhard(eqpl+deqpl,sa,dsi,dsk)
    backstresscur = backstressprev
    
    f = smises - syield0
    
    if f>toler:
      syield = syield0
      counter = 0.
      while (np.abs(f)>toler)&(counter<10):
        f = smises-emod*deqpl-syield-backstresscur+backstressprev
        deqpl=deqpl+f/(emod+dsi+dsk)
        syield,backstresscur = combhard(eqpl+deqpl,sa,dsi,dsk)
      
      if np.abs(f)>toler:
        print "*** WARNING:: newton on plasticity did not converge"
      
      eqpl += deqpl
      alpha += (backstresscur-backstressprev)*flowdir
      plasticcur += flowdir*deqpl
      elasticcur -= flowdir*deqpl
      stresscur = flowdir*syield + alpha
    
    elasticprev=elasticcur
    plasticprev=plasticcur
    alphaprev = alpha
    backstressprev = backstresscur
  
    stress += [stresscur]
    strain += [logstrain]
    plastic += [plasticcur]
    elastic += [elasticcur]
    eqplhist += [eqpl]
   
  return np.array(stress),np.array(strain)
  

def combhard(eqpl,sa,dsi,dsk):
  si = sa+dsi*eqpl
  sk = dsk*eqpl
  return si,sk

