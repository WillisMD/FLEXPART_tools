# -*- coding: utf-8 -*-
# ######################################### #
# Time-stamp: <2017-02-22 mdwillis> #
# ######################################### #

'''Reads output of NETCARE 2015 FLEXPART-ECMWF (run for NETCARE by D. Kunkel at UMainz) simulations for tracing plumes along flight paths. Contact Megan or Daniel to get the FLEXPART data.
Returns total residence times for each plume that is initialized, and saves an HDF5 file for later use. Adjust the saved variables for your own purposes.

To run this script:
plot_tracer_netcare_plumes_backward.py YYYYMMDD HHMMSS nflight

Adapted from D. Kunkel's code (Main_plot_tracer_netcare_backward.py)
'''

# ######################################### # runs with:
import numpy as np                          # Version: 1.11.1
import matplotlib.pyplot as plt             # Version: 1.5.1
from matplotlib import cm                   # Version: 1.5.1
import matplotlib.colors as mcolors         # Version: 1.5.1
from mpl_toolkits.basemap import Basemap    # Version: 1.0.7
import h5py                                 # Version: 2.6.0
import sys                                  # 
import os                                   # 
# ######################################### #

########################################
######----External Input:----###########
########################################
if len(sys.argv) == 4:
    date  = sys.argv[1] # DATE in YYYYMMDD
    itime = sys.argv[2] # TIME in HHMMSS, this is the time in each file name
    nflight = str(sys.argv[3]) # Flight name e.g., Flight1NC2015
else:
    print ("                         ")
    print (" START HELP PAGE         ")
    print ("                         ")
    print (" ----------------------- ")
    print (" HOW TO USE THIS SCRIPT: ")
    print (" ----------------------- ")
    print (" PLEASE TYPE:            ")
    print ("                         ")
    print (" python ",sys.argv[0], "date itime nflight")
    print ("                         ")
    print (" date  = DATE in YYYYMMDD")
    print (" itime = TIME in HHMMSS")
    print (" nflight= flight name")
    print (" date and itime have to be as given in the netcdf file name")
    quit(" END HELP PAGE ")
# -------------------------------------------------#
print ("                                           ")
print (" ***************************************** ")
print ("    START OF FILE                          ")
print ("     ", sys.argv[0])
print (" ***************************************** ")

########################################
######----INPUT PATHS:----##############
########################################
# base path
path0="C:\\Users\\JA3\\Documents\\FLEXPART_output_2015\\"

# path to data of specific flight
path1=path0+nflight+"b\\"
# filename of flexpart output file
fname="grid_time_"+str(date)+str(itime)+".nc"
# full path to flexpart output file
pfile=path1+fname
# 
print ('*** file: '+pfile)

# path where the flight tracks are saved:
ftpath=path0+"flight_tracks\\"

# outpath for h5 files
opath=path0+"PES\\"

# check if data paths exist, create it if not
xtest=os.path.exists(opath)
if xtest == False:
    os.mkdir(opath)

# trafile
trafile=path1+"trajectories_new.txt"
xtest=os.path.exists(trafile)
if xtest == False:
    # run a bash script with a sed command
    # to replace all \- with \s\-
    print ("!!! please run sed \"s/\-/\ \-/g\" trajectories.txt >> trajectories_new.txt")
    print ("!!! to replace all minus with a space+minus in the trajectories.txt file")
    quit( "!!! trajectories_new.txt is missing -> quit" )

# general marker of the campaign
marker="netcare"

# tracer name
var="spec001_mr"

####################################################
######----- Start Function Definitions------########
####################################################
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx #array[idx]
    
####################################################   
def get_flight_track(dpath,ident):
    '''Get the right flight track for each flight. Flight tracks should be in ftpath. Called as flightname=get_flight_track(ftpath,nflight)'''
    files=os.listdir(dpath) #return a list containing the names of the files in a directory
    for n in range(0, len(files)): #loop over all the files in this folder
        if ident==files[n][9:]: #checking that we have the right flight track file
            print ("*** flight track file: "+str(files[n]))
            fn=dpath+files[n] #A string with the full path to a particular flight track file (str(date)+"_"+nflight)
            break
    return fn #full path to flight track
    
####################################################
def read_flight_track(pflight):
    '''Reads flight track data. Called as flight_data = read_flight_track(flightname). flightname is the output of get_flight_track. Output: time, lon, lat, alt'''
    date,time,lon,lat,alt=np.loadtxt(pflight,unpack=True)#unpack transposes the output so that we can assign it to variables
    return [date,time,lon,lat,alt] #we are returning a list here, so flight_data is a list of arrays (floats by default)

####################################################
def read_trajectories(ptra):
    '''Reads trajectories file from flexpart. Called as read_trajectories(trafile) i.e., path to trajectories_new.txt. Outputs data in trajectory file'''
    # open file and read general header
    with open(ptra,'r') as ft: #this syntax ensures that the files closes at the end of the code chunk without problems
    # read first line into the variable ft (general info about sim start and flexpart version)
        line1=ft.readline() #not used
    # skip second line (parameter: method, lsubgrid, lconvection)
        line2=ft.readline() #not used
    # third line gives the number of release points == number of center of mass trajectories
        line3=ft.readline()
    ntra = int(line3) # == number of center of mass trajectories as an integer
    # now read start info of release points == start of center of mass trajectorys
    with open(ptra, 'r') as ft: #open the file again...
        relinfo = [] # initialize empty lists
        relname = []
        tradata = []
        for i in range(0,3): #for i = 0, 1, 2 skip to the next line
            ft.readline() #skip the first three lines
        ii = 0 
        for i in range(0,2*ntra): #for i = 0,..., 2*ntra-1
            help1=ft.readline() #store the value in help1
            if i%2: #does this say if i is odd? (i.e., i%2 = true or non-zero?)
                relname.append(help1)
                ii=ii+1
            else:
                relinfo.append([x for x in help1.split()]) #help1 is a string and we are splitting it here
        for line in ft:
            tradata.append(line) #append each line to tradata
    # tradata: only number of release points, time, lon, lat, alt ---> these are the release points along the flight track
    tracenter = np.zeros((5, len(tradata))) #returns a 5 by len(tradata) array filled with zeros
    for i in range(0,len(tradata)):
        tracenter[:,i] = tradata[i].split()[0:5] 
    
    return [relname, relinfo, tracenter] #tracenter is an array of the release points(time, lat lon, alt), put into the variable centertra at call

####################################################
def get_index0(ndim, idx):
    '''ndim: length of dimension, idx: current index. Output: sidx= adjusted index as string, adds leading zeros to index if necessary'''
    if ndim < 100 and ndim >=10:
        if int(idx) < 10:
            sidx = "0"+str(int(idx))
        else:
            sidx = str(int(idx))
    elif ndim < 1000 and ndim >=100:
        if int(idx) < 10:
            sidx = "00"+str(int(idx))
        elif int(idx) > 9 and int(idx) < 100:
            sidx = "0"+str(int(idx))
        else:
            sidx = str(int(idx))
    elif ndim < 10 and ndim >= 1:
        sidx = str(int(idx))
    else:
        quit("!!! no time available" )
    return sidx

####################################################
def get_var(pfile1,var0): #var0 is a string
    '''function to get variables from the NetCDF file, except the species given in var (i.e., the tracer)'''
    print ('*** variable: '+var0+' ***')
    ff = h5py.File(pfile1, mode='r') #ff is an instance of a HDF5 file object
    # if you want to check which variables are in the file, uncomment the following:
    #ff.visititems(print)
    dset = ff[var0] #dset is still an hdf5 file object
    data=np.array(dset) #make dset a numpy array so we can do something with it
    del dset
    ff.close()
    return data #returns whatever variable as an np.array

####################################################
def get_spec(pfile1,var0,dimlist,idx_sum):
    '''Funtion to get 'var' = spec001_mr. The issue is that var is too large to get at once, so the retrieval is split up. dimlist = [nage,npoint,ntime,nlev,nlat,nlon]
    Called as: trac = get_spec(pfile, var, [nage,npoint,ntime,nlev,nlat,nlon],2)'''
    print ("***  get ", var0, " from file "+pfile1[50:]+" with dimensions ", dimlist, " ***")
    ff = h5py.File(pfile1, mode='r') # ff is an HDF5 file object, we only 'open' it once
    if dimlist[0]==1: #if nage = 1(only one ageclass)
        if int(idx_sum) != 2: #the only thing we can sum over is the time dimension
            print ("!!! not yet implemented")
            quit("!!! only mean over time allowed" ) # this code is actually taking the mean residence time
        # set idx_sum to an integer value
        idx = int(idx_sum)
        nlen = int(dimlist[idx]) #in current code dimlist[idx] = dimlist[2] = ntime
        # initialize a 4D array of the correct size, containing zeros
        dset = np.zeros((dimlist[1],dimlist[3],dimlist[4],dimlist[5])) #create a 4D array of zeros npoint x nlev x nlat x nlon
        
        #Loop to sum over the times for all the releases during one flight
        for i in range(0,nlen): # nlen is about equal to 10 days (ntime) -- looping over all the times backwards from the release point
            ds = ff[var0][0,:,i,:,:,:] # from the tracer (var) in ff, grab the whole array at time stamp i
            dset = ds + dset #add the existing dset to ds --> this is summing over all the times backwards from the release point
            print ("*** Num and shape: ", i, dset.shape, " with max values (ds, dset/nlen): ", i, ds.max(), dset.max()/float(nlen)) #print info about each array being added
        #dset = dset/float(nlen) #finally divide dset by nlen (total time) --- this gives the mean residence time over the past 10 days in each grid cell
        
    else:
        print( "!!! not yet implemented")
        quit( "!!! more than one age class" )
    data=np.array(dset) #make dset an np.array (npoint x nlev x nlat x nlon) ---> 4D array that is the residence time summed over time for each release along the flight track
    del dset
    ff.close()
    return data 
    
####################################################
######----- Endt Function Definitions------#########
####################################################


####################################################
######----------------MAIN-----------------#########
####################################################
def main():
    ### MAIN:

    ###########################
    ######-get grid info-######
    ###########################
    # latitude
    lat= get_var(pfile,"latitude") #1D array --> shape = (260,)
    nlat = len(lat)
    # longitue
    lon= get_var(pfile,"longitude") #1D array --> shape = (1440,)
    nlon = len(lon)
    # levels
    lev= get_var(pfile,"height") #10 levels/heights in meters are [50, 100, 200, 500, 1000, 2000, 4000, 6000, 8000, 10000]
    nlev = len(lev)
    # times
    time= get_var(pfile,"time") #tracking particles every 6 hours back in time from release, up to 10.5 days --> shape = (42,)
    ntime = len(time)
    # numspec
    nspec=get_var(pfile, "numspec") #shape = (1,)
    numpoint=get_var(pfile, "numpoint") #this seems to be identical to pointspec
    ageclass=get_var(pfile, "nageclass")
    nage = len(ageclass)
    pointspec=get_var(pfile, "pointspec") #this is a 1D array of zeros of size (nreleases,)
    npoint = len(pointspec)
    # release heights --> +/- 25m in z from aircraft altitude
    # bottom
    relbot=get_var(pfile,"RELZZ1")
    nrelbot=len(relbot)
    # top
    reltop=get_var(pfile,"RELZZ2")
    nreltop=len(reltop)
    assert nreltop == nrelbot #number of releases is constant
    #calculate aircraft altitude at release
    P6alt = relbot + ((reltop-relbot)/2) #aircraft altitude at release
    
    relstart=get_var(pfile,"RELSTART") #related to the time of release, but unclear what this is relative to, maybe itime?
    relend  =get_var(pfile,"RELEND")

    #dimlist = [nage, npoint, ntime, nlev, nlat, nlon]
    dimlist = [nage, npoint, 20, nlev, nlat, nlon] #to integrate over 1 day ntime=4, to integrate over 5 days ntime=20
    #dimlist = [nage, npoint, 32, nlev, nlat, nlon] #currently set to sum only to 15.5 days back for consistency
    
    print("Dimensions of tracer variable = " + str(dimlist))
   
   
    ###########################
    ########-get tracer-#######
    ###########################
    #4 dimensions of the tracer, with shape [npoint,nlev,nlat,nlon]
    tracer = get_spec(pfile, var, dimlist, 2)
   
    tracer_tcol = np.sum(tracer, axis=1) #sum over nlev, this is now a 3D [npoint,nlat,nlon] array holding total column PES for all the releases along flight track
    
    temp0 = tracer[:,0:2,:,:]
    tracer_pcol100= np.sum(temp0, axis=1)
    del(temp0)
    
    temp1 = tracer[:,0:3,:,:]
    tracer_pcol200= np.sum(temp1, axis=1)
    del(temp1)
    
    temp2 = tracer[:,0:5,:,:]
    tracer_pcol1000=np.sum(temp2, axis=1)
    del(temp2)
    
    temp3 = tracer[:,0:6,:,:]
    tracer_pcol2000=np.sum(temp3, axis=1)
    del(temp3)
    
    ###########################
    ########-get traj.-########
    ###########################
    #read flight track
#    flightname=get_flight_track(ftpath,nflight)
#    #flight data contains date,time,lon,lat,alt
#    flight_data = read_flight_track(flightname)
    #read trajectories
    tra_data = read_trajectories(trafile)
    center_traj=tra_data[2]
  
    ###################################################################
    ###########--- Save all outputs to an HDF5 file---#################
    ###################################################################
    outfile = opath + "summedPES_"  +nflight+ "_5d.h5" 

    outHDF5 = h5py.File(outfile, mode = 'w')
    outHDF5.create_dataset('PES', data = tracer)
    outHDF5.create_dataset('PES_tCol', data = tracer_tcol)
    outHDF5.create_dataset('PES_pCol100', data = tracer_pcol100)
    outHDF5.create_dataset('PES_pCol200', data = tracer_pcol200)
    outHDF5.create_dataset('PES_pCol1000', data = tracer_pcol1000)
    outHDF5.create_dataset('PES_pCol2000', data = tracer_pcol2000)
    outHDF5.create_dataset('latitude', data = lat)
    outHDF5.create_dataset('longitude', data = lon)
    outHDF5.create_dataset('zlevels', data = lev)
    outHDF5.create_dataset('aircraft_alt', data = P6alt)
    outHDF5.create_dataset('plume_center', data = center_traj)
    outHDF5.create_dataset('Nreleases', data = npoint)
    outHDF5.create_dataset('time', data = time)
    outHDF5.create_dataset('RelEnd', data = relend)
    outHDF5.create_dataset('RelStart', data = relstart)
    
    outHDF5.close()  
    
####################################################
######----------Execute Main-------------###########
####################################################   
main() 
####################################################
#########----------THE END-------------#############
####################################################