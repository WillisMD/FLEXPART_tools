#!/usr/bin/python
# ######################################### #
# Time-stamp: <2016-07-07 12:21:54 dkunkel> #
# ######################################### #
'''
DESCRIPTION:
Read FLEXPART output of simulations for tracing
plumes along flight paths. 
Returns mean residence times for each plume that
is initialized. 

PLEASE GO THROUGH THE PART INPUT PATHS and adjust
the variables for your needs!

USAGE OF THE SCRIPT:
python plot_tracer_netcare_plumes_backward.py YYYYMMDD HHMMSS FLIGHT_IDENT

Edited by Megan Willis to be compatible with Python3.5 (2016/08/25)
'''
# ######################################### # runs with:

import numpy as np                          # Version: 1.8.2
import matplotlib.pyplot as plt             # Version: 1.4.2
from matplotlib import cm                   # Version: 1.4.2
import matplotlib.colors as mcolors         # Version: 1.4.2
from mpl_toolkits.basemap import Basemap    # Version: 1.0.7
import h5py                                 # Version: 2.2.1
import sys                                  # Version:
import os                                   # Version:
# ######################################### #

# External input:
if len(sys.argv) == 4:
    date  = sys.argv[1] # DATE in YYYYMMDD
    itime = sys.argv[2] # TIME in HHMMSS --> this might need to be a string for when itime = 000000, follow up on this later
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
    print (" date and itime have to be as given in the netcdf file name    ")
    print ()
    print (" ANY QUESTIONS? -> dkunkel@uni-mainz.de                        ")
    print (" Last revised: 07.07.2016                                      ")
    print ("                                                               ")
    quit(" END HELP PAGE ")
# -------------------------------------------------
print ("                                           ")
print (" ***************************************** ")
print ("    START OF FILE                          ")
print ("     ", sys.argv[0])
print (" ***************************************** ")
# ######################################### #
#
# INPUT PATHS:
#
# ######################################
# base path
path0="C:\\Users\\JA3\\Documents\\FLEXPART_output_2015\\"
# ######################################
# path to data of specific flight
path1=path0+nflight+"b\\"
# filename of flexpart output file
fname="grid_time_"+str(date)+str(itime)+".nc"
# full path to flexpart output file
pfile=path1+fname
# 
print ('*** file: '+pfile)
# ######################################
# path where the flight tracks are saved:
ftpath=path0+"flight_tracks\\"
# ######################################
# outpath where to put the pdf/png
opath=path1+"my_plots\\"
# check if data paths exist
xtest=os.path.exists(opath)
if xtest == False:
    os.mkdir(opath)
# ######################################
# trafile
trafile=path1+"trajectories_new.txt"
xtest=os.path.exists(trafile)
if xtest == False:
    # run a bash script with a sed command
    # to replace all \- with \s\-
    print ("!!! please run sed \"s/\-/\ \-/g\" trajectories.txt >> trajectories_new.txt")
    print ("!!! to replace all minus with a space+minus in the trajectories.txt file")
    quit( "!!! trajectories_new.txt is missing -> quit" )
# ######################################### #
# general marker of the campaign
marker="netcare"
# ######################################### #
# tracer name
var="spec001_mr"

# ######################################### #
# END DEFINITION PART 
# ######################################### #

# ######################################### #
#### ROUTINES FOR GETTING DATA:

def get_flight_track(dpath,ident):
    # get the right flight track
    # for each flight
    # flight tracks should be in the dir given below
    #
    files=os.listdir(dpath) #return a list containing the names of the files in a directory
    for n in range(0, len(files)):
        if ident==files[n][9:]:
            print ("*** flight track file: "+str(files[n]))
            fn=dpath+files[n] #I think this is a string?
            break
    return fn #full path to flight track

def read_flight_track(pflight):
    # read flight track data
    # input: path to file name file (fn, also flightname)
    # output: time, lon, lat, alt
    date,time,lon,lat,alt=np.loadtxt(pflight,unpack=True)#unpack transposes the output so that we can assign it to variables
    return [date,time,lon,lat,alt] #we are returning a list here, so flight_data(L335) must be a list containing arrays (floats by default)
#
def read_trajectories(ptra):
    # read trajectories file from flexpart
    # input : path to tra file
    # output: data in tra file
    # open file
    # read general header
    with open(ptra,'r') as ft: #this syntax ensures that the files closes at the end of the code chunk without problems
    # read first line into the variable ft (general info about sim start and flexpart version)
        line1=ft.readline()
    # skip second line (parameter: method, lsubgrid, lconvection)
        line2=ft.readline()
    # third line gives the number of release points
    # == number of center of mass trajectories
        line3=ft.readline()
    #
    ntra = int(line3) # == number of center of mass trajectories as an integer
    #
    # now read start info of release points == start of center of mass trajs
    with open(ptra, 'r') as ft: #open the file again...
        relinfo = [] # initialize empty lists
        relname = []
        tradata = []
        for i in range(0,3): #for i = 0, 1, 2 skip to the next line
            ft.next() #skip to the next line...not sure what this is doing?
        ii = 0 #what is ii for??
        for i in range(0,2*ntra): #for i = 0,..., 2*ntra-1
            help1=ft.next() #store the value in help1
            if i%2: #does this say if i is odd? (i.e., i%2 = true or non-zero?)
                relname.append(help1)
                ii=ii+1
            else:
                relinfo.append([x for x in help1.split()]) #help1 is a string and we are splitting it here
        for line in ft:
            tradata.append(line) #append each line to tradata
            
    # tradata: only number of release point, time, lon, lat, alt ---> these are the release points along the flight track
    tracenter = np.zeros((5, len(tradata))) #returns a 5 by len(tradata) array filled with zeros
    for i in range(0,len(tradata)):
        tracenter[:,i] = tradata[i].split()[0:5] #use square braces to index an array
    
    return [relname,relinfo, tracenter] #tracenter is an array of the release points(time, lat lon, alt), put into the variable centertra at call

#
def get_index0(ndim, idx):
    # input:
    # ndim: length of dimension
    # idx: current index
    #
    # output:
    # sidx: adjusted index as string
    # adds leading zeros to index,
    # if necessary
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


def get_var(pfile1,var0): #var0 is a string
    # function to get variables
    # except the species given in var (ie, the tracer)
    #print '*** '+pfile1+' ***'
    print ('*** variable: '+var0+' ***')
    ff = h5py.File(pfile1, mode='r') #ff is an instance of a HDF5 file object
    # if you want to check which variables are
    # in the file, uncomment the following:
    #print ff.keys()
    dset = ff[var0] #dset is still an hdf5 file object
    data=np.array(dset) #make dset a numpy array so we can do something with it
    del dset
    ff.close()
    return data #returns whatever variable as an np array

def get_spec(pfile1,var0,dimlist,idx_sum):
    # funtion to get var --- i.e., the tracer 'var' = spec001_mr
    # the issue is that var has too large to get everything of
    # it at once, so the retrieval is splitted
    #dimlist = [nage,npoint,ntime,nlev,nlat,nlon]
    print ("***  get ", var0, " from file "+pfile1[50:]+\
        " with dimensions ", dimlist, " ***")
    ff = h5py.File(pfile1, mode='r') # ff is an HDF5 file object, we only 'open' it once
    if dimlist[0]==1: #if nage = 1(only one ageclass)
        # if len(idx_sum) > 1:
        #     print "not yet implemented"
        #     quit( "more than one dimension to sum" )
        if int(idx_sum) != 2: #the only thing we can sum over is the time dimension
            print ("!!! not yet implemented")
            quit("!!! only mean over time allowed" ) # this code is actually taking the mean residence time
            
        # set idx_sum to an integer value
        idx = int(idx_sum)
        nlen = int(dimlist[idx]) #in current code dimlist[idx] = dimlist[2] = ntime
        
        # needed: dynamical array allocation
        dset = np.zeros((dimlist[1],dimlist[3],dimlist[4],dimlist[5])) #create a 4D array of zeros npoint x nlev x nlat x nlon
        
        ############---- Big for loop to sum over the times for one release along the flight track -----######
        for i in range(0,nlen): # nlen is about equal to 10 days (ntime) -- looping over all the times backwards from the release point
            # needed dynamical allocation
            ds = ff[var0][0,:,i,:,:,:] # from the tracer (var) in ff, grab the whole array at time stamp i
            dset = ds + dset #add the existing dset to ds --> this is summing over all the times backwards from the release point
            print ("*** Num and shape: ", i, dset.shape, \
                " with max values (ds, dset/nlen): ",\
                i, ds.max(), dset.max()/float(nlen)) #print info about each array being added
        dset = dset/float(nlen) #finally divide dset by nlen (total time) --- this gives the mean residence time over the past 10 days in each grid cell
        
    else:
        print( "!!! not yet implemented")
        quit( "!!! more than one age class" )
    data=np.array(dset) #make dset an np.array (npoint x nlev x nlat x nlon) ---> 4D array that is the residence time summed over time for each release along the flight track
    del dset
    ff.close()
    return data #assigned to trac later

# #################################### #
# END FUNCTION SECTION
# #################################### #

# #################################### #
# MAIN PART OF THE PROGRAM:
# #################################### #

def main():
    ### MAIN:

    ###### get grid info 
    # latitude
    lat= get_var(pfile,"latitude") #1D arrays
    nlat = len(lat)
    # longitue
    lon= get_var(pfile,"longitude")
    nlon = len(lon)
    # level
    lev= get_var(pfile,"height")
    nlev = len(lev)
    # time
    time= get_var(pfile,"time")
    ntime = len(time)
    # numspec
    nspec=get_var(pfile, "numspec")
    numpoint=get_var(pfile, "numpoint") #this seems to be identical to pointspec
    ageclass=get_var(pfile, "nageclass")
    nage = len(ageclass)
    pointspec=get_var(pfile, "pointspec") #this is a 1D array of zeros of size (nreleases,)
    npoint = len(pointspec)
    # release heights
    # bottom
    relbot=get_var(pfile,"RELZZ1")
    nrelbot=len(relbot)
    # top
    reltop=get_var(pfile,"RELZZ2")
    nreltop=len(reltop)
    assert nreltop == nrelbot #number of heights relative to bottom and top should be equal - otherwise break??
    
    P6alt = relbot + ((reltop-relbot)/2) #aircraft altitude at release
    
    relstart=get_var(pfile,"RELSTART")
    relend  =get_var(pfile,"RELEND")

    ###### get full tracer
    trac = get_spec(pfile, var, [nage,npoint,ntime,nlev,nlat,nlon],2)
    # new tracer size
    # 4 dimensions, currently: [npoint,nlev,nlat,nlon] 
    # so now plot only first level (0-200m) and total column
    # for the total column --> sum over nlev
    # --> ptrac 4D array
    ptrac = np.sum(trac,axis=1) #sum over nlev, this is now a 3D [npoint,nlat,nlon] array holding total column PES for all the release times
    # this is for the total column, sums array elements over a given axis where 1=rows, 0=columns ---> i.e., sum trac over its rows <---
    ptracl= trac[:,0,:,:]       # this is for the lowest level, just level zero --- should be just up to 50m?
    '''
    EXAMPLE to obtain data for specific alitutde ranges
    ---
    If you want to plot the level up to a specific altitude, e.g., 200m 
    you have to sum up the lowest levels
    ! --------------------------
    level of alitudes in the netcdf files (upper borders):
    height = 50, 100, 200, 500, 1000, 2000, 4000, 6000, 8000, 10000 ;
    ! --------------------------
    --> so for 0-200 m you need the first 3 levels:
    ! --------------------------
    phelp = trac[:,0:3,:,:]
    ptracl= np.sum(phelp,axis=1)
    ! --------------------------
    then ptracl holds the data for 0-200 m
    '''
    del trac
    #print ptrac.shape

    ##### read flight track data
    # file for flighttrack
    #flightname="/home/dkunkel/output/2016_julia_netcare15_flexpart/pyscripts/"+str(date)+"_"+nflight
    flightname=get_flight_track(ftpath,nflight)

    ##### flight data contains date,time,lon,lat,alt
    flight_data = read_flight_track(flightname)
    ###### read trajectories
    tra_data = read_trajectories(trafile)
    centertra=tra_data[2]

    ######## START THE PLOT ROUTINE FOR EACH RELEASE POINT:
    for jd in range(0,npoint): #pull each release out of ptrac in a for loop ptrac [npoint, nlat, nlon]
        # p0 should be a 3D array
        # at one time step jd
        # with p0[numpoints, nlat, nlon]
        p0 = ptrac[jd,:,:]
        p1 = ptracl[jd,:,:]
        #print "*** Some info: ", numpoint[jd],p0.shape, \
        #    p1.shape, p0.max(),p0.min()

        ############ GET CENTER OF MASS TRAJECTORY
        # get lon / lat of center tra
        clon = []
        clat = []
        calt = []
        for i in range(0, len(centertra[0,:])):
            if int(centertra[0,i]) == int(jd)+1 and centertra[1,i] < 9.e5:

                clon.append(float(centertra[2,i]))
                clat.append(float(centertra[3,i]))
                calt.append(float(centertra[4,i]))

        # get the plume number as string:
        sjd = get_index0(npoint,jd)
        # specify the output file name:
        outfile=opath+marker+"_plume_"+sjd+"_"+var+"_"+str(date)+str(itime)+"_backward.pdf"
                    
        ######################################### #
        ######################################### #
        # plot data
        # ######################################### #
        
        print ("*** PLOT DATA "+sjd+" "+outfile)
        # ######################################### #
        # ######################################### #
        # initialize the plot
        # set up a 2x2 plot matrix
        # so each of the release altitudes gets one panel
        #
        #   title
        #  --------
        # |  TC    |
        # |________|
        # |0-50 m |
        # |________|
        #  color bar           
        #
        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(5.5,11)) #name the figure fig and the subplots ax
        # set the grid
        lons, lats = np.meshgrid(lon,lat) #make coordinate matrices from coordinate vecotors (for 2014, the header contains coordinate matrices already)
        # # set the colorlevels
        levels=([0.01,0.02,0.05,0.1,0.2,0.5,1.0,2.0,5.0,10.,20.])
        # set the colormap --> yellow orange red
        cols=cm.YlOrRd
        ### adjust colors
        norm = mcolors.BoundaryNorm(levels,ncolors=cols.N, clip=False)
        ### set the variable and the title of each panel:
        for jn in range(0,2): #for jn=0,1 i.e., for the two panels, set pp depending on whether we want the total or partial column
            if jn==0:
                pp = p0
                ii=0
                stitle="Total column (0-10 km)"
            elif jn==1:
                pp = p1
                ii=1
                stitle="0-50 m"
            
            # setup polar stereo plot
            m = Basemap(ax=ax[ii],projection='npstere',
                        resolution='c',boundinglat=40.,lon_0=-133.4)
            # add coastlines 
            m.drawcoastlines()
            # add latitudes
            parallels = np.arange(0.,75,15.)
            m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
            # add longitudes
            meridians = np.arange(0.,360.,30.)
            m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)
            # set the boarders
            x, y = m(lons, lats)               
            # set the title:
            ax[ii].set_title(stitle, fontsize=14)
            # add the date
            if jn==0:
                # time
                simstarth=int(itime[0:2])
                simstartm=int(itime[2:4])
                simstarts=int(itime[4:6])
                #
                # print simstarth,simstartm,simstarts
                # release time relative to absolute
                reh = int(relend[jd]/3600.)#convert to hours (release end hours)
                rem = int((relend[jd]/3600. - int(relend[jd]/3600.))*60./100. * 100.)
                teh = simstarth + reh
            
                if teh < 0:
                    teh = 24+teh
                tem = simstartm + rem
                #print reh, rem, teh, tem
                if tem < 0:
                    tem = 60 + tem
                    teh = teh - 1
                    if teh < 0:
                        teh = 24+teh
                #
                # print teh, tem
                steh = get_index0(24, teh)
                stem = get_index0(60, tem)
                st = steh+" "+stem+" h" #create a string that is the releas time
                # add the release time to the figure
                ax[ii].text(0.15, 2., st, color="black", \
                            bbox=dict(facecolor='white', \
                                      pad=2, edgecolor='white'),  fontsize=16)
            #
            # do the contour filled plot
            cs0 = ax[ii].contourf(x,y,pp[:,:], levels=levels,\
                                  cmap=cols,extend="max",norm=norm) #pp is either p0 or p1 for total or partial column
            # white color below last levels contour
            cs0.cmap.set_under('w') 
            # add the flight path
            xf,yf = m(flight_data[2],flight_data[3])
            cf=ax[ii].plot(xf,yf,'-',color='blue' ,linewidth=.5)
            # add the center of mass trajectories
            if jn==0: # only in the total column plot
                xc,yc = m(clon,clat)
                cc=ax[ii].plot(xc,yc,'-',color='cyan' ,linewidth=2.)


            # # axis labels
            if jn==1:
                ax[ii].set_xlabel(r'longitude (deg)', \
                                  fontsize=18, labelpad=15)
                ax[ii].set_ylabel(r'latitude (deg)', \
                                  fontsize=18, labelpad=20)
            # # add a colorbar
            if jn==1:
                # add_axes[left, bottom, width, height]  
                cbar_ax = fig.add_axes([0.15, 0.04, 0.7, 0.015])
                cbar=plt.colorbar(cs0, cax=cbar_ax, \
                                      format="%.2f", orientation="horizontal")
                cbar.set_label('residence time (s)')
            # delete the plot variable
            del pp
        # save figure
        fig.savefig(outfile, dpi=150)
        # print the name of the outfile:
        # delete more variables
        del p0, p1, xc, yc, xf, yf
        # close current figure:
        plt.close(fig)
        plt.clf()
    del ptrac, ptracl
    #
# #########################################
# execute main part of the program:       #
# #########################################
main()                                    #
# #########################################
# FINISH                                  #
# ####################################### #
quit("\n !!! THIS IS IT !!!\n" )          #
# ####################################### #
