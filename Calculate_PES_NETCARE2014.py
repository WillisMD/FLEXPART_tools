# -*- coding: utf-8 -*-
"""
Calculates vertically resolved PES as a residence time, plus total & partial column PES, from LATMOS FLEXPART-WRF output (run for NETCARE by J.L. Thomas)
Created on Thu Aug 25, 2016
@author: Megan Willis

To run this script:
Calculate_PES_NETCARE2014.py YYYYMMDD = flight date FLIGHT_ID = e.g., flight-05July NUMRELEASES = integer number of releases for the flight
"""
################################
import numpy as np                         
import matplotlib.pyplot as plt            
from matplotlib import cm                  
import matplotlib.colors as mcolors         
from mpl_toolkits.basemap import Basemap   
import h5py                                 
import sys                                  
import os 
#################################

#########################
####-EXTERNAL INPUT:#####
#########################
if len(sys.argv) == 4:
    date  = sys.argv[1] # DATE in YYYYMMDD
    flight_ID = str(sys.argv[2]) # Jennie's flight name for FLEXPART file structure
    Nreleases = int(sys.argv[3])
    #IntTime = int(sys.argv[4])# Use this input arguement if time integration less than Ntimes is required --> change line 153
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
n_header = "header_d01.nc" #names

####################################################
######-------- Function Definitions--------#########
####################################################

def get_var(pfilex,variable): #pfilex is the path to the file, variable is a string (the name of the variable)
    # function to get variables (mostly from the header.nc file), except for the tracer
    print ('*** '+pfilex+' ***')
    print ('*** variable: '+variable+' ***')
    ff = h5py.File(pfilex, mode='r') #ff is an instance of a HDF5 file object
    # if you want to check which variables are in the file, uncomment the following:
    #print ff.visititems(print)
    dset = ff[variable] #dset is still an hdf5 file object
    data=np.array(dset) #make dset a numpy array so we can do something with it
    del dset
    ff.close()
    return data #returns whatever variable as an np array


def get_tracer(pfiley, variable, dimlist, Ntimes):
    print ("***  getting ", variable, " from file "+pfiley+" with dimensions ", dimlist, " ***")
    ff = h5py.File(pfiley, mode='r')
    dset = np.zeros((dimlist[3], dimlist[4], dimlist[5]))
    
    for i in range(1,Ntimes):#start with the second time step, the first one contains zero concentration
        
        if i%1000 == 0: #check progress
            print(i)
        
        ds = ff[variable][i,0,0,:,:,:] #The variable CONC has shape (time, ageclass, releases, height(Z  - top/bottom), latitude (X - north/south), longitude (Y - west/east) )
        dset = ds + dset
        
    data = np.array(dset)
    del(dset)
    ff.close()
    return data
        
####################################################
######----------------MAIN-----------------#########
####################################################
for i in range(0, Nreleases):
    dfcount = 1 #initialize this as 1
    dfcount = dfcount + i*10 #ten minutes between each output 
    dfolder = path1+"output_"+str(dfcount).zfill(5)+"\\" 
    p_header = dfolder + n_header
    
    ###############################################################################
    #############----look at the file header-------################################
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

    ################################################################################
    ############---other ways to interogate what is in the file:--------############
    ################################################################################
    #f_head.visit(print) #this is the same as using a for loop
    #for name in ff:
        #print(name) #see the names of all the variables in ff, seems like there is just one group in the header file
    #f_head_attr = list(ff.attrs.items()) 
    #print(f_head_attr) #print the file attributes    
    #f_head_keys = list(ff.keys()) #some things for exploring what is in the file
    #print(f_head_keys)
    #ff_values = list(ff.values())  #includes shape and data type, also works: items()
    #print(ff_values)
    ################################################################################
    ###############--------build the data file name:---------------#################
    ################################################################################
    n_data = "flxout_d01_"+str(date)+"_"+start_str+ ".nc" 
    p_data = dfolder + n_data
    print("***Getting data from "+ p_data+"***")
    ###############################################################################
    #############----look at the main Flexpart output .nc file-------##############
    ###############################################################################
    f_tracer = h5py.File(p_data, mode = 'r')
    #f_tracer.visititems(print)
    ###############################################################################
    ###########--- get grid and dimension information---###########################
    ###############################################################################
    xlat = get_var(p_header, 'XLAT')
    Nlat = len(xlat)
     
    xlon = get_var(p_header, 'XLONG')
    Nlon = len(xlon)
    
    ztop = get_var(p_header, 'ZTOP') #Ztop = top of vertical levels in meters = [10, 100, 200, 300, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 20,000]
    Ntop = len(ztop)
    
    times = get_var(p_data, 'Times')
    Ntimes = times.shape[0] #same as len(times)
    
    ageclass = get_var(p_data, 'ageclass')
    Nage = len(ageclass)
    
    releases = get_var(p_data, 'releases')
    Nrel = len(releases)
    
    dimlist = [Ntimes, Nage, Nrel, Ntop, Nlat, Nlon]
    print("Dimensions of CONC variable = " + str(dimlist))
    
    ###############################################################################
    ###########--- get PES(residence time) 3D and 2D arrays---#####################
    ###############################################################################
    
    Conc = get_tracer(p_data, 'CONC', dimlist, Ntimes) #put Ntimes here for the whole amount of time, or (e.g.,) 5700 for 4 days back
    #Conc is the vertically resolved PES, a 3D array
    Conc_tCol = np.sum(Conc, axis=0)
    
    temp1 = Conc[0:4,:,:]
    Conc_pCol300 = np.sum(temp1, axis=0)
    del(temp1)
    
    temp2 = Conc[0:1,:,:]
    Conc_pCol10 = np.sum(temp2, axis=0)
    del(temp2)
    
    temp3 = Conc[0:2,:,:]
    Conc_pCol100 = np.sum(temp3, axis=0)
    del(temp3)
    
    temp4 = Conc[0:3,:,:]
    Conc_pCol200 = np.sum(temp4, axis=0)
    del(temp4)
    
    temp5 = Conc[0:5,:,:]
    Conc_pCol500 = np.sum(temp5, axis=0)
    del(temp5)

    ###############################################################################
    ###########--- Save all three outputs to an HDF5 file---#######################
    ###############################################################################
    outfile = dfolder + "summedPES_" +str(date)+ "_" +start_str+ ".h5" #change the output file name when changing total summing time

    outHDF5 = h5py.File(outfile, mode = 'w') 
    outHDF5.create_dataset('Conc_tsum', data = Conc)
    outHDF5.create_dataset('Conc_tCol', data = Conc_tCol)
    outHDF5.create_dataset('Conc_pCol300', data = Conc_pCol300)
    outHDF5.create_dataset('Conc_pCol10', data = Conc_pCol10)
    outHDF5.create_dataset('Conc_pCol100', data = Conc_pCol100)
    outHDF5.create_dataset('Conc_pCol200', data = Conc_pCol200)
    outHDF5.create_dataset('Conc_pCol500', data = Conc_pCol500)
    
    
    outHDF5.close()
 
###############################################################################
###########------------------- END------------------###########################
###############################################################################

    

    

 
 