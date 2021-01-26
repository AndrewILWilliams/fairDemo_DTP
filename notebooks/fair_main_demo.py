import numpy as np
from scipy.optimize import root

def fair_scm(emissions=False,
             other_rf=0.0,
             q=np.array([0.33,0.41]),
             tcrecs=np.array([1.6,2.75]),
             d=np.array([239.0,4.1]),
             a=np.array([0.2173,0.2240,0.2824,0.2763]),
             tau=np.array([1000000,394.4,36.54,4.304]),
             r0=32.40,
             rc=0.019,
             rt=4.165,
             F_2x=3.74,
             C_0=278.0,
             ppm_gtc=2.123,
             iirf_max=97.0,
             restart_in=False,
             restart_out=False):

  # If TCR and ECS are supplied, calculate the q1 and q2 model coefficients 
  # (overwriting any other q array that might have been supplied)
  # ref eq. (4) and (5) of Millar et al ACP (2017)
  k = 1.0 - (d/70.0)*(1.0 - np.exp(-70.0/d))
  if type(tcrecs) in [np.ndarray,list]:
    q =  (1.0 / F_2x) * (1.0/(k[0]-k[1])) * np.array([tcrecs[0]-tcrecs[1]*k[1],tcrecs[1]*k[0]-tcrecs[0]])

  #Set up the output timeseries variables
  # emissions must be a numpy array for this to work
  if type(emissions) in [np.ndarray,list]:
    carbon_boxes_shape = tuple(list(emissions.shape) + [4])
    thermal_boxes_shape = tuple(list(emissions.shape) + [2])
    integ_len = len(emissions)
  elif type(other_rf) in [np.ndarray,list]:
    carbon_boxes_shape = tuple(list(other_rf.shape) + [4])
    thermal_boxes_shape = tuple(list(other_rf.shape) + [2])
    integ_len = len(other_rf)
    emissions = np.zeros(integ_len)
  else:
    raise ValueError("Neither emissions or other_rf is defined as a timeseries")

  RF = np.zeros(integ_len)
  C_acc = np.zeros(integ_len)
  iirf = np.zeros(integ_len)
  R_i = np.zeros(carbon_boxes_shape)
  T_j = np.zeros(thermal_boxes_shape)

  C = np.zeros(integ_len)
  T = np.zeros(integ_len)

  if restart_in:
    R_i[0]=restart_in[0]
    T_j[0]=restart_in[1]
    C_acc[0] = restart_in[2]

  else:
    #Initialise the carbon pools to be correct for first timestep in numerical method
    R_i[0,:] = a * emissions[0,np.newaxis] / ppm_gtc

  C[0] = np.sum(R_i[0,:],axis=-1)

  if type(other_rf) == float:
    RF[0] = (F_2x/np.log(2.)) * np.log((C[0] + C_0) /C_0) + other_rf
  else:
    RF[0] = (F_2x/np.log(2.)) * np.log((C[0] + C_0) /C_0) + other_rf[0]

  if restart_in == False:
    #Update the thermal response boxes
    T_j[0,:] = (q/d)*(RF[0,np.newaxis])

  #Sum the thermal response boxes to get the total temperature anomlay
  T[0]=np.sum(T_j[0,:],axis=-1)

  for x in range(1,integ_len):
      
    #Calculate the parametrised iIRF and check if it is over the maximum allowed value
    iirf[x] = rc * C_acc[x-1]  + rt*T[x-1]  + r0
    if iirf[x] >= iirf_max:
      iirf[x] = iirf_max
      
    #Linearly interpolate a solution for alpha
    if x == 1:
      time_scale_sf = (root(iirf_interp_funct,0.16,args=(a,tau,iirf[x])))['x']
    else:
      time_scale_sf = (root(iirf_interp_funct,time_scale_sf,args=(a,tau,iirf[x])))['x']

    #Multiply default timescales by scale factor
    tau_new = tau * time_scale_sf

    #Compute the updated concentrations box anomalies from the decay of the pervious year and the additional emisisons
    R_i[x,:] = R_i[x-1,:]*np.exp(-1.0/tau_new) + a*(emissions[x,np.newaxis]) / ppm_gtc
    #Summ the boxes to get the total concentration anomaly
    C[x] = np.sum(R_i[...,x,:],axis=-1)
    #Calculate the additional carbon uptake
    C_acc[x] =  C_acc[x-1] + 0.5*(emissions[x] + emissions[x-1]) - (C[x] - C[x-1])*ppm_gtc

    #Calculate the total radiative forcing
    if type(other_rf) == float:
      RF[x] = (F_2x/np.log(2.)) * np.log((C[x] + C_0) /C_0) + other_rf
    else:
      RF[x] = (F_2x/np.log(2.)) * np.log((C[x] + C_0) /C_0) + other_rf[x]

    #Update the thermal response boxes
    T_j[x,:] = T_j[x-1,:]*np.exp(-1.0/d) + q*(1-np.exp((-1.0)/d))*RF[x,np.newaxis]
    #Sum the thermal response boxes to get the total temperature anomaly
    T[x]=np.sum(T_j[x,:],axis=-1)

  if restart_out:
    restart_out_val=(R_i[-1],T_j[-1],C_acc[-1])
    return C + C_0, T, restart_out_val
  else:
    return C + C_0, T

def iirf_interp_funct(alp_b,a,tau,targ_iirf):
	# ref eq. (7) of Millar et al ACP (2017)
    iirf_arr = alp_b*(np.sum(a*tau*(1.0 - np.exp(-100.0/(tau*alp_b)))))
    return iirf_arr   -  targ_iirf

def load_rcp_emissions(emms_file,Y_max=2200):

  #Loads a CO2 emissions timeseries
  emm_data = np.loadtxt(emms_file,delimiter=',')
  #Extract time, FF+LUC CO2 emissions
  i_out=emm_data[:,0]<=Y_max
  return emm_data[i_out,0], emm_data[i_out,1]+emm_data[i_out,2]

def load_rcp_forcing(forc_file,Y_max=2200):

  #Loads a radiative forcing timeseries
  forc_data = np.loadtxt(forc_file,delimiter=',')
  #Extract time, total, anthropogenic, GHG and CO2 radiative forcing
  i_out=forc_data[:,0]<=Y_max
  return forc_data[i_out,0], forc_data[i_out,1] , forc_data[i_out,4], forc_data[i_out,5], forc_data[i_out,8]

def load_hadcruttemps(d_file,Y_max=2016):

  hc_data = np.genfromtxt(d_file)

  i_out=hc_data[:,0]<=Y_max
  return hc_data[i_out,0],hc_data[i_out,1]

def damage_function(D0, gamma, T):
    T1=np.zeros_like(T)
    T1[np.where(T > 0)]=T[np.where(T > 0)]
    return ((D0*T1**gamma)/(1+D0*T1**gamma))
        




