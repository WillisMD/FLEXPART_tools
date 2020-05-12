# -*- coding: utf-8 -*-
"""
Reads h5 file with summed FLEXPART output and NSIDC netCDF file to compare
    PES with sea ice concentration and calculate residence time over open water
A lot of this code came from Calculate_PES_NETCARE2014.py and WRFseaice_OpenWater_ResidenceTime.py
Created on October 13, 2016
@author: meganwillis
"""

################################
import numpy as np                         
import matplotlib.pyplot as plt            
from matplotlib import cm                  
import matplotlib.colors as mcolors         
from mpl_toolkits.basemap import Basemap  
import h5py    
from netCDF4 import Dataset   
import sys    
import pandas as pd 
import datetime    
import pyresample                   
#################################

#########################
####-EXTERNAL INPUT:#####
#########################
if len(sys.argv) == 4:
    date  = sys.argv[1] # DATE in YYYYMMDD
    flight_ID = str(sys.argv[2]) # Jennie's flight name for FLEXPART file structure
    Nreleases = int(sys.argv[3])
else:
    print('You forgot the input arguements!')

#########################
####-INPUT PATHS:########
#########################
# base path
path0="C:\\Users\\JA3\\Documents\\FLEXPART_output_2014\\"
#path to a specific flight
path1 = path0+flight_ID+"-backward\\"
#set header file name
n_header = "header_d01.nc" 

####################################################
######-------- Function Definitions--------#########
####################################################

def get_var(pfilex,variable): #pfilex is the path to the file, variable is a string (the name of the variable)
    # function to get variables (mostly from the header.nc file), except for the tracer
    ff = h5py.File(pfilex, mode='r') #ff is an instance of a HDF5 file object
    # if you want to check which variables are in the file, uncomment the following:
    #print ff.visititems(print)
    data=np.array(ff[variable]) #make a numpy array so we can do something with it
    ff.close()
    return data #returns whatever variable as an np array
    
def getvar_nc(pfilex,variable): #pfilex is the path to the file, variable is a string (the name of the variable)
    # function to get variables from the WRF output file (NETCDF3_64BIT_OFFSET) and other NetCDF3 files (e.g., NSIDC sea ice)
    ff = Dataset(pfilex,'r') #ff is an instance of a netCDF file object
    # if you want to check which variables are in the file, uncomment the following:
    #print ff.visititems(print)
    data = np.array(ff[variable])
    ff.close()
    return data #returns whatever variable as an np array
    
    
####################################################
######----------------MAIN:-----------------########
####################################################
#initialize two numpy arrays and then put them into a dataframe
dates  = pd.date_range(date, periods=Nreleases)
Tot_res = np.zeros(shape=(Nreleases))
Rel_res = np.zeros(shape=(Nreleases))
dataout = pd.DataFrame({'t_series':dates, 'residence_time':Tot_res, 'relative_residencetime':Rel_res})

for n in range(0, Nreleases):
    dfcount = 1 #initialize this as 1
    dfcount = dfcount + n*10 #ten minutes between each output 
    dfolder = path1+"output_"+str(dfcount).zfill(5)+"\\" 
    p_header = dfolder + n_header
    
    print('*** release number'+str(n)+'***')
    
    ###############################################################################
    #############----look at the FLEXPART file header:-------######################
    ###############################################################################
    f_head = h5py.File(p_header, mode='r')
    #f_head.visititems(print) #print out items in the HDF file
    f_head_attr = list(f_head.attrs.items())#read attributes to a list
    #get the start time that should be in the data file name
    start_time = f_head_attr[9][1] #array of size (1,)
    temp_str = str(start_time).strip('[]')
    if temp_str[2:4] == '00':
        start_str = str(start_time-4100).strip('[]')
    else:
        start_str = str(start_time-100).strip('[]')
    
    print('***start time = ' +start_str)
    
    ################################################################################
    ###############--------build the data file names:---------------################
    ################################################################################
    #summed PES HDF5 file created from Calculate_PES_NETCARE2014.py
    n_Fdata = "summedPES_10days_"+str(date)+"_"+start_str+ ".h5"
    p_Fdata = dfolder + n_Fdata
    #NSIDC daily data for July 2014
    n_dailydata = "NSIDC-0051_91602_daily.sub.nc"
    p_dailydata = "C:\\Users\\JA3\\Documents\\NSIDC_seaice_dailyJuly2014\\" + n_dailydata
    
    ###############################################################################
    #########################--- get variables:---#################################
    ###############################################################################
    #get the grid information
    xlat = get_var(p_header, 'XLAT')
    Nlat = len(xlat)  
    xlon = get_var(p_header, 'XLONG')
    Nlon = len(xlon)   
    ztop = get_var(p_header, 'ZTOP')
    Ntop = len(ztop)
    
    xlatice = getvar_nc(p_dailydata, 'latitude')
    xlonice = getvar_nc(p_dailydata, 'longitude')
    times = getvar_nc(p_dailydata, 'time')# this is a time in days from 1601-01-01
    #Ntimes = len(times)
    
    #get the previously summed PES information - check Calculate_PES_NETCARE2014.py for how I named the variables in this file
    Conc_pCol = get_var(p_Fdata, 'Conc_pCol300') 
    #Conc_tCol = get_var(n_Fdata, 'Conc_tCol') #total column
    #Conc_tsum = get_var(n_Fdata, 'Conc_tsum') ##D array summed over time
    
    #get NSIDC data for the appropriate day
    dayindex = int(date[6:8]) - 1
    ff = Dataset(p_dailydata, mode='r') 
    temp = np.array(ff['Sea_Ice_Concentration_with_Final_Version'][dayindex,:,:])
    seaice = temp.astype(dtype='float64')
    ff.close()
    
    ###############################################################################
    #########################--- fix  the sea ice:---##############################
    ###############################################################################
    #make land areas and missing data have sea ice concentration of NaN
    #From http://nsidc.org/data/docs/daac/nsidc0051_gsfc_seaice.gd.html
    # 0-250 = sea ice concentrations, 251 = mask to cover data gap at the pole, 252 = unused, 253 = coastlines, 254 = superimposed land mask, 255 = missing data
    it = np.nditer(seaice, flags = ['multi_index'], op_flags = ['readwrite']) #iterate over the elements of the np.array seaice
    
    for i in it: #for each element in the iterator (it's an array, I think)
        if seaice[it.multi_index] == 251. or seaice[it.multi_index]==252. or seaice[it.multi_index]==253. or seaice[it.multi_index]==255.: #remove various flags
            seaice[it.multi_index] = np.nan 
        elif seaice[it.multi_index] == 254.:
            seaice[it.multi_index] = np.nan 
            
    seaice_norm = seaice/250  #sea ice fraction goes from 0 to 1
    
    #re-grid NSIDC sea ice 
    orig_def = pyresample.geometry.SwathDefinition(lons=xlonice, lats=xlatice)
    targ_def = pyresample.geometry.SwathDefinition(lons=xlon, lats=xlat)
    seaice_regrid = pyresample.kd_tree.resample_nearest(orig_def, seaice_norm, targ_def, radius_of_influence=500000, fill_value=None)
    
    ###############################################################################
    #########################--- sum the residence times:---#######################
    ###############################################################################
    
    it2 = np.nditer(seaice_regrid, flags = ['multi_index'], op_flags = ['readwrite'])
    
    restime = np.float64(0) #initialize residence time as zero
    for j in it2:
        if seaice_regrid[it2.multi_index] < 0.2:
            #print(seaice_regrid[it2.multi_index])
            restime += Conc_pCol[it2.multi_index]
    
    print("{0:.15f}".format(restime))
    dataout.loc[n, 'residence_time'] = restime #put the final number for this release into an array for later
    
    rel_restime = restime/(np.sum(Conc_pCol))
    dataout.loc[n, 'relative_residencetime']=rel_restime
    
    ###############################################################################
    #####################---create a datetime object:---###########################
    ###############################################################################
    strdate = str(date)
    time = datetime.datetime.strptime(strdate+start_str,'%Y%m%d%H%M%S')
    
    dataout.loc[n,'t_series'] = time 
    
    
 ##############################################################
 #####################---save output:---#######################
 ##############################################################

outfile = path1 + "NSIDC_10days_openwater_"+flight_ID+".csv"
print('*** Saving to :'+ outfile)
dataout.to_csv(outfile)
    
