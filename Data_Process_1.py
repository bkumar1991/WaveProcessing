# This script is for pre-processing surface wave data
# read in mseed data and stationxml
# remove channel response
# resample the data with sampling_rate=1
# remove trend and mean
# rotate NE to RT (need event information to calculate back_azimuth)
# save the Z and T component data as SAC files by the order of epicentral distance


import numpy as np
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
clientfdsn = Client("IRIS")
from obspy import read, read_inventory
from obspy.clients.iris import Client
clientiris = Client()
from obspy.core.stream import Stream
from obspy.core.event import read_events
from obspy.core.trace import Trace
from obspy.io.sac.sactrace import SACTrace
from obspy.core.inventory import Inventory
import sys,os
import pandas as pd
import csv
import shutil # shutil.rmtree('path') to remove all files in the directory
import glob #(find file names with matching chacracters)
from obspy.geodetics.base import gps2dist_azimuth, kilometer2degrees

yeardata = "./Example_Data/"
catalog = 'catalog-1'
CatFile = pd.read_csv(yeardata + catalog + '.csv')
#print(CatFile.columns)

## loop over events in the catalog
nev = len(CatFile)
#print(nev)
for i in range(nev):
#for i in range(0,2):   # test for individual event
    otime = UTCDateTime(CatFile.DateTime[i])
    evlat = CatFile.Latitude[i]
    evlon = CatFile.Longitude[i]
    evdp = CatFile.Depth[i]
    evmg = CatFile.Magnitude[i]
    print(nev, otime, evlat, evlon, evdp, evmg)
    # generate path name by event            
    evotime = otime.isoformat()
    evotime = evotime[0:10] + '-' + evotime[11:13] + '-' + evotime[14:16]
    evpath = yeardata + evotime
    mseedpath = evpath + '/waveforms/'
    stationxml = evpath + '/stations/*xml'
    #print(mseedpath)
    #print(stationxml)

    ##  make path for LHZ and LHT and erase the existing path and subdirectories
    LHZpath = evpath + '/LHZ/'
    LHTpath = evpath + '/LHT/'
    if not os.path.exists(LHZpath):
        os.makedirs(LHZpath)
    else:
        shutil.rmtree(LHZpath)
        os.makedirs(LHZpath)
    
    if not os.path.exists(LHTpath):
        os.makedirs(LHTpath)
    else:
        shutil.rmtree(LHTpath)
        os.makedirs(LHTpath)
    # Read in mseed files and stationxml files from event directory     
    inv = read_inventory(stationxml)
    inv0 = inv.copy()
    nsta = len(inv)
    #print("the number of stations:", nsta)
    st = read(mseedpath + '*mseed')

    st0 = st.copy()
    # print("the number of traces:", len(st))
    
    # loop over stations to get the Z and T components
    
    stZ = Stream()
    stT = Stream()
    fig1, ax1= plt.subplots(figsize= (6.5,4)) 
    # fig2, ax2= plt.subplots(figsize= (6.5,4))   
    for j in range(nsta):
    #for j in range(2,3):   # test for individual station
        stlat = inv0[j].stations[0].latitude
        stlon = inv0[j].stations[0].longitude
        stel = inv0[j].stations[0].elevation
        stcode = inv0[j].stations[0].code
        stnet = inv0[j].code
        # print(j, stcode, stnet, stlat, stlon, stel)
        
        # calculate back_azimuth for rotating components EN to RT 

        
        dist_m, az, bacaz = gps2dist_azimuth(stlat, stlon, evlat, evlon)
        dist = kilometer2degrees(dist_m / 1000.0)
        # distbaz = clientiris.distaz(stalat=stlat, stalon=stlon, evtlat=evlat, evtlon=evlon)
        # bacaz = distbaz['backazimuth']
        # dist = distbaz['distance']
                
        # form a new stream for one station based on station name
        st1 = st.select(station = stcode)
        #st1.plot(type = 'relative' )


        # Only keep stations with three components

        if len(st1) == 3:                        
            
           # remove response and define a filter band to prevent apmlifying noise during deconvolution
           
            st1.remove_response(inventory=inv, water_level=60, plot=False, zero_mean=True)
            #st1.remove_response(inventory=inv, water_level=60, plot=True, pre_filt = None, taper_fraction= 0, zero_mean=True)
 
           # Find the common time span
            common_start = max(trace.stats.starttime for trace in st1)
            common_end = min(trace.stats.endtime for trace in st1)

           # Trim all traces to the common time span        
            st1.trim(starttime=common_start, endtime=common_end)
           

            #remove trend
            st1.detrend('spline', order=3, dspline=300, plot=False)
            #st1.detrend('spline', order=3, dspline=300, plot=True)
            
            # apply a broadband filters, using zerophase = True to avoid phase shift
            st1.filter('bandpass', freqmin = 0.002, freqmax = 0.1, zerophase = 'True')
            #st1.plot()
            # scale the amplitdues
            for k in range(3):
                st1[k].data = st1[k].data*1e7

            # remove stations that have traces without data (zero amplidues) 
            if (st1[0].max() == 0 or st1[1].max() == 0 or st1[2].max() == 0):
                ampratio1 = 0.
                ampratio2 = 0.
            else:    
                ampratio1 = abs(st1[0].max()/st1[2].max())
                ampratio2 = abs(st1[1].max()/st1[2].max())
            
            # if the data are ok, go ahead to rotate NEZ-RTZ    
            if ampratio1 > 0.1 and ampratio1 < 10 and ampratio2 > 0.1 and ampratio2 < 10:
            
                st2 = st1.rotate(method='NE->RT', back_azimuth=bacaz)
                #st2.plot()

                sacZ = SACTrace.from_obspy_trace(st2[2])
                sacT = SACTrace.from_obspy_trace(st2[0])
                
                # add header information
                sacZ.evla = evlat; sacZ.evlo = evlon; sacZ.evdp = evdp; sacZ.mag = evmg; sacZ.stel = stel 
                sacZ.stla = stlat; sacZ.stlo = stlon; sacZ.gcarc = dist; sacZ.baz = bacaz; sacZ.o = 0
                sacT.evla = evlat; sacT.evlo = evlon; sacT.evdp = evdp; sacT.mag = evmg; sacT.stel = stel 
                sacT.stla = stlat; sacT.stlo = stlon; sacT.gcarc = dist; sacT.baz = bacaz; sacT.o = 0

                # from SAC to obspy trace with header information
                trZ = sacZ.to_obspy_trace()
                trT = sacT.to_obspy_trace()
                
                # add distance in meters for the section plot
                trZ.stats.distance = trZ.stats.sac.gcarc*110*1000
                trT.stats.distance = trT.stats.sac.gcarc*110*1000

                # form new stream for Z and T components
                stZ.append(trZ)
                stT.append(trT)
            
    ## make plots and save in files for all Z components and T copmonents, respectively
    # stZ.plot(type = 'section', time_down = True,shaw = True)
    # stT.plot(type = 'section', time_down = True, shaw = True)
    stZ.plot(type = 'section',fig=fig1,  time_down = True, scale=1.0, linewidth=0.6)
    ax1=fig1.axes[0]
    ax1.set_xlabel('Epicentral Distance (Km)', fontsize=10)
    ax1.set_ylabel('Time (s)', fontsize=10)
    fig1.tight_layout()
    plt.text(0.02, 1.04, '(a)',transform=plt.gca().transAxes, fontsize=10, va='top',ha='left')
    plt.savefig(evpath + '-Z.svg', bbox_inches='tight')
    
    fig2, ax2= plt.subplots(figsize=(6.4,4))
    stT.plot(type = 'section',fig=fig2,  time_down = True, scale=1.0, linewidth=0.6)
    ax2=fig2.axes[0]
    ax2.set_xlabel('Epicentral Distance (Km)', fontsize=10)
    ax2.set_ylabel('Time (s)', fontsize=10)
    fig2.tight_layout()
    plt.text(0.02, 1.04, '(b)',transform=plt.gca().transAxes, fontsize=10, va='top',ha='left')
    plt.savefig(evpath + '-T.svg', bbox_inches='tight')
           
    
    
    # sort traces by distance (for Z and T )
    # wrie out the files in SAC format in the order of distance
    
    #for rl in range (1,2): # Love wave
    #for rl in range(1): # Rayleigh  wave
    for rl in range(2): # both       
        if rl == 0:
           ntr = len(stZ)
           stX = stZ.copy()
           Datapath = LHZpath
           print("Rayleigh wave")

        if rl == 1:
           ntr = len(stT)
           stX = stT.copy()
           Datapath = LHTpath
           print("Love wave")

        # form an array for all epicentral distances
        adist = []
        for j in range(ntr):
            adist = np.append(adist,np.array([[stX[j].stats.distance]]))
        # return the indices of the sorted distance array 
        ind = np.argsort(adist)
            
        # write out the SAC files for Z components in the order of distance (near to far)
        for j in range(ntr):
            net = stX[ind[j]].stats.network
            sta = stX[ind[j]].stats.station
            
            
            if j < 9:
                fname = Datapath + 'D00' + str(j+1) + '.' + net + '.' + sta
                stX[ind[j]].write(fname + '.SAC', format = 'SAC')
                
            if j >= 9 and j < 99:
                fname = Datapath + 'D0' + str(j+1) + '.' + net + '.' + sta
                stX[ind[j]].write(fname + '.SAC', format = 'SAC')
            else:
                fname = Datapath + 'D'+ str(j+1) + '.' + net + '.' + sta
                stX[ind[j]].write(fname + '.SAC', format = 'SAC')

            if(j<9):
                tfname = Datapath + 'D' + str(j+1) + '.' + net + '.' + sta + '.SAC'
                if os.path.exists(tfname):
                    os.remove(tfname)
