import numpy as np

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

def load_rcp_emissions(emms_file,Y_max=2200):
    #Loads a CO2 emissions timeseries
    emm_data = np.loadtxt(emms_file,delimiter=',')
    #Extract time, FF+LUC CO2 emissions
    i_out=emm_data[:,0]<=Y_max
    return emm_data[i_out,0], emm_data[i_out,1]+emm_data[i_out,2]