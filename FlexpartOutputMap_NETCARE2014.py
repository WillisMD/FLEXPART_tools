
# -*- coding: utf-8 -*-
"""
Creates a stereographic projection of FLEXPART-WRF output from HDF5 files created from Calculate_PES_NETCARE2014.py (put these files in your working directory)
Created on Wed Aug 31 13:30:52 2016
@author: meganwillis
"""

################################
import numpy as np                         
import matplotlib.pyplot as plt            
from matplotlib import cm                  
import matplotlib.colors as mcolors         
from mpl_toolkits.basemap import Basemap  
import h5py 
import pandas as pd                                
#################################

#########################
####-INPUT PATHS:########
#########################

n_header = 'header_d01.nc' #put a full path here, or change your working directory
n_data = 'summedPES_20140712_203925.h5' #this is output from Calculate_PES_2014.py

p_Flighttracks = '/Users/meganwillis/Documents/Data/FlightTrack_Maps/FlightTracks2014.csv'
####################################################
######-------- Function Definitions--------#########
####################################################

def get_var(pfilex,variable): #pfilex is the path to the file, variable is a string (the name of the variable)
    # function to get variables (mostly from the header.nc file), except for the tracer
    #print ('*** '+pfilex+' ***')
    #print ('*** variable: '+variable+' ***')
    ff = h5py.File(pfilex, mode='r') #ff is an instance of a HDF5 file object
    # if you want to check which variables are in the file, uncomment the following:
    #print ff.visititems(print)
    dset = ff[variable] #dset is still an hdf5 file object
    data=np.array(dset) #make dset a numpy array so we can do something with it
    del dset
    ff.close()
    return data #returns whatever variable as an np array

####################################################
######----------------MAIN-----------------#########
####################################################  
#get flight tracks
FlightTracks = pd.DataFrame.from_csv(p_Flighttracks)

                   
#define the figure
plt.clf()
fig = plt.plot(figsize=(10,10))

#get the data
#get the grid information
xlat = get_var(n_header, 'XLAT')
Nlat = len(xlat)  
xlon = get_var(n_header, 'XLONG')
Nlon = len(xlon)   
ztop = get_var(n_header, 'ZTOP')
Ntop = len(ztop)
#get the previously summed PES information
#Conc_pCol = get_var(n_data, 'Conc_pCol300') #change this or comment it out
Conc_tCol = get_var(n_data, 'Conc_tCol') #total column
                       
#create a basemap instance
m = Basemap(width=4000000,height=4000000, resolution='l',projection='stere',lat_0=74,lon_0=-90.)
#m = basemap(projection='npstere',resolution='c',boundinglat=40.,lon_0=-133.4)

#re-project the FLEXPART grid
x,y=m(xlon,xlat)

#draw the map details
m.drawcoastlines(linewidth=0.5)
m.drawmapboundary(fill_color = 'lightcyan', zorder = 0)
m.fillcontinents(color = 'oldlace', lake_color = 'lightcyan', zorder=0)
parallels = np.arange(0.,81,10.) # labels = [left,right,top,bottom]
m.drawparallels(parallels,labels=[False,True,True,False])
meridians = np.arange(10.,351.,20.)
m.drawmeridians(meridians,labels=[True,False,False,True])

#make the PES contour plot 
levels = ([0.01,0.05,0.25,0.5,2.5,5.0,25.,50.,250.,500.,2500.])#levels=([0.01,0.02,0.05,0.1,0.2,0.5,1.0,2.0,5.0,10.,20.])
cols=cm.jet
norm = mcolors.BoundaryNorm(levels,ncolors=cols.N, clip=False)
cs = m.contourf(x,y,Conc_tCol, levels=levels, cmap=cols, extend='max', norm=norm) # this code is currently only plotting the total column PES

#add a scale bar
#cbar_ax = fig.add_axes([0.15, 0.04, 0.7, 0.015])
cbar=plt.colorbar(cs, format="%.2f", orientation="horizontal",fraction=.06, pad=0.08)
cbar.set_label('FLEXPART-WRF total column residence time (s)')

#add one flight track
long_f6 = np.array(FlightTracks.loc[::10, 'Longitude_deg6'])
lat_f6 = np.array(FlightTracks.loc[::10, 'Latitude_deg6'])
xf6,yf6 = m(long_f6, lat_f6)
m.plot(xf6,yf6, color='grey', linewidth=0, marker='.', markersize=0.5 )

plt.savefig("PES_20140712_203925.png", format="png", dpi=300)
plt.show()