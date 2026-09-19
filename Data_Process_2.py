# This script reads in SAC files produced by Data_process_1.py and applys filters
# filter data (Z and T) at selected frequencies
# window the filtered data based on energy location at the diseried f
# could also try using autopick or triggering for surface wave
# clean the data before and after the window and form a new trace
# save the new trace as sac format 
# the file name should include event, station, and filter information
####################################################################

import numpy as np
import matplotlib.pyplot as plt
from obspy.clients.fdsn import Client
clientfdsn = Client("IRIS")
from obspy import UTCDateTime
from obspy import read
from obspy.clients.iris import Client
clientiris = Client()
from obspy import read, read_inventory
from obspy.core.stream import Stream
from obspy.core.event import read_events
from obspy.core.trace import Trace
from obspy.io.sac.sactrace import SACTrace
from obspy.core.inventory import Inventory
import sys,os
from obspy.signal.filter import envelope
from obspy.signal.trigger import trigger_onset, plot_trigger
import pandas as pd
import csv
import shutil # shutil.rmtree('path') to remove all files in the directory
import glob #(find file names with matching chacracters)

def pickmin(a,npw,npt,itmax):
    # pick the minimum amplitude value for the envelope that is closet to the largest amplitude
    # a: the array
    # itmax: the index of the largest amplitdue 
    # npt: the number of points before and after itmax for searching the minimum value
    abmin = a[itmax]
    aamin = a[itmax]
    ibmin = 0
    iamin = 0
    for i in range(npw):
     ib = itmax - i
     ia = itmax + i
     if ib < 1:
         ib = 1
     if (a[ib] < abmin and ibmin == 0):
         abmin = a[ib]
         if abmin <= a[ib -1]:
             ibmin = ib
     if ia > npt-2:
         ia = npt -2
     if (a[ia] < aamin and iamin == 0):
         aamin = a[ia]
         if aamin <= a[ia + 1]:
             iamin = ia

     if (ibmin > 0 and iamin >0):
         break

    if (ibmin == 0):
         ibmin = itmax - npw
         
    if (iamin == 0):
         iamin = itmax + npw
   
    return ibmin, iamin
         

## define the array containing the interested central frequencies (Hz = (1/T(s))
freq = []
# freq = [0.005, 0.0055, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011, 0.0125, 0.015, 0.0175,
#         0.02, 0.0225, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05]
freq = [ 0.0225]
nfrq = len(freq)
## define the array for Rayleigh and Love wave phase velocities based on the global model
# Rvel = [4.6, 4.5, 4.42, 4.31, 4.23, 4.18, 4.15, 4.12, 4.09, 4.04, 4.01,
#         3.98, 3.95, 3.91, 3.83, 3.75, 3.67, 3.59, 3.52]
# Lvel = [4.9, 4.84, 4.79, 4.72, 4.66, 4.61, 4.57, 4.54, 4.49, 4.43, 4.37,
#         4.31, 4.27, 4.20, 4.11, 4.02, 3.95, 3.89, 3.84]
Rvel = [ 4.09]
Lvel = [ 4.49]

## working direcroty
#fdir = "./python codes"
yeardata = "./Example_Dataset/"
catalog = 'catalog-1'
CatFile = pd.read_csv(yeardata + catalog + '.csv')
#print(CatFile.columns)

fw = open(catalog + '.output', 'a')  #change 'w' to 'a' if want to append to the existing file at the end of the file

#print(CatFile[0:1])
## loop over events in the catalog
nev = len(CatFile)
#print(nev)
for i in range(nev): 
    otime = UTCDateTime(CatFile.DateTime[i])
    evlat = CatFile.Latitude[i]
    evlon = CatFile.Longitude[i]
    evdp = CatFile.Depth[i]
    evmg = CatFile.Magnitude[i]
    print(nev, otime, evlat, evlon, evdp, evmg)
    fw.writelines(str(otime) + '\n')  # output the event info
    # generate path name by event            
    evotime = otime.isoformat()
    evotime = evotime[0:10] + '-' + evotime[11:13] + '-' + evotime[14:16]
    evpath = yeardata + evotime

    # directory LHZ and LHT should have been created using Data_Process_1.py
    LHZpath = evpath + '/LHZ/'
    LHTpath = evpath + '/LHT/'
    
    # Read in SAC files from the event directory     
    
    stZ = read(LHZpath + '/D*SAC', fromat = 'SAC')
    stT = read(LHTpath + '/D*SAC', fromat = 'SAC') 
    nsta = len(stZ) # the number of traces for Z and T should be the same.

    for j in range(nsta):
        stZ[j].stats.distance = stZ[j].stats.sac.gcarc*110*1000
        stT[j].stats.distance = stT[j].stats.sac.gcarc*110*1000

    # stZ.plot(type = 'section', time_down = True,shaw = True)
    # stT.plot(type = 'section', time_down = True, shaw = True)

    distn = stZ[0].stats.distance
    distf = stZ[nsta-1].stats.distance
    
    #loop over Rayleigh and Love wave,

    #for rl in range (1,2): # Love wave
    #for rl in range(1): # Rayleigh  wave
    for rl in range(2): # both       
        if rl == 0:
           stX = stZ.copy()
           Datapath = LHZpath
           vel = Rvel.copy()
           print("Rayleigh wave")
           fw.writelines('  Rayleigh wave    No. traces: ' + str(nsta) +  '\n')

        if rl == 1:
           stX = stT.copy()
           Datapath = LHTpath
           vel = Lvel.copy()
           print("Love wave")
           fw.writelines('  Love wave    No. traces: ' + str(nsta) +  '\n')

        
## loop over frequencies
        
        for k in range(nfrq):  
        #for k in range(7,9): 
            stXflt = Stream()
            stXfltnew = Stream()
            stXfltkeep = Stream()
            Fltpath = Datapath + 'f-' + str(freq[k])
            if not os.path.exists(Fltpath):
               os.makedirs(Fltpath)
            else:   
               shutil.rmtree(Fltpath)
               os.makedirs(Fltpath)
               
            freqmin = freq[k] - 0.005
            freqmax = freq[k] + 0.005
            
            if freqmin == 0:
                freqmin = 0.001
                freqmax = freq[k] + 0.004
         # the long surface wave window
            wnb =  otime + float("{:.0f}".format(distn/(1000*vel[0]))) - 240
            wne =  otime + float("{:.0f}".format(distf/(1000*vel[nfrq-1]))) + 480
      
       
                
            #loop over traces to determine surface wave window more accurately based on phase arrival of 
            # the largest ampitdue (tmax). 
            #divide each trace to three segments: beofre, surface-wave, and after.
            #Multiply zero to the before and after segments and merge the three togeter.
            
            Rttaa = [[],[],[],[]] #array to store [tph, tmax, Amax, Aavg, dist] 
            Lttaa = [[],[],[],[]] #
            # loop over traces
            
            for l in range(nsta):
            #for l in range(5):    
                tr = stX[l]
                net = stX[l].stats.network
                sta = stX[l].stats.station
                npts = stX[l].stats.npts
                samprate = stX[l].stats.sampling_rate
                twnb = (wnb-otime)/samprate 
                tph = tr.stats.distance/(1000*vel[k])  # phase arrival time
                print("network \n", net)
                print("station \n", sta)
                trflt =tr.copy()
                trflt.filter('bandpass', freqmin = freqmin, freqmax = freqmax, zerophase = 'True')
                ## plot figure a

                # figa = plt.figure(figsize=(6.5, 4))
                # tr.plot(type = 'section',fig=figa, starttime= otime, linewidth=1.0, color='k', orientation = 'horizontal')
                # plt.text(0.02, 1.04, '(a)',transform=plt.gca().transAxes, fontsize=10, va='top',  ha='left')
                # plt.gca().yaxis.set_visible(False)
                # plt.gca().spines['top'].set_visible(False)
                # plt.gca().spines['right'].set_visible(False)
                # plt.gca().spines['left'].set_visible(False)
                # plt.xlabel("Time (s)")
                # plt.gca().xaxis.grid(False)
                # plt.show() 

                ## plot figure b

                # figb= plt.figure(figsize=(6.5,4))
                # trflt.plot(type = 'section',fig=figb,  starttime = otime, linewidth=1.0, orientation = 'horizontal')
                # plt.text(0.02, 1.04, '(b)',transform=plt.gca().transAxes, fontsize=10, va='top',  ha='left')
                # plt.gca().yaxis.set_visible(False)
                # plt.gca().spines['top'].set_visible(False)
                # plt.gca().spines['right'].set_visible(False)
                # plt.gca().spines['left'].set_visible(False)
                # plt.xlabel("Time (s)")
                # plt.gca().xaxis.grid(False)
                # plt.show()
                # use trace envelope, sort the ampitdue, find their time
                # the largest amplitude should be close to the phase arrival time
                # the ratio of the 2nd largest amplitude to the largest amplitdue should be wihtin a threshold
                trflt0 = trflt.slice(starttime = wnb, endtime = wne) 
                stXflt.append(trflt0)

                

                evlp = envelope(trflt0.data)
                tmax = np.argmax(evlp)
                Amax = evlp[tmax]
                Aavg = np.average(evlp)
                
                tmin1, tmin2 = pickmin(evlp,240,len(evlp),tmax)
                
                tmax = twnb + tmax
                tmin1 = twnb + tmin1
                tmin2 = twnb + tmin2

#                print(sta,float("{:.0f}".format(tr.stats.distance/1000)), float("{:.4f}".format((freqmin+freqmax)/2.0)), \
#                      float("{:.1f}".format(Amax/Aavg)), float("{:.1f}".format(tmax-tph)))

                if rl == 0:
                    Rttaa[0].append(tph), Rttaa[1].append(tmax), Rttaa[2].append(Amax), Rttaa[3].append(Aavg) 
                
                if rl == 1:
                    Lttaa[0].append(tph), Lttaa[1].append(tmax), Lttaa[2].append(Amax), Lttaa[3].append(Aavg) 

                
                tevlp = np.arange(0, len(trflt0)/samprate, 1/samprate) + twnb


                ## plot figure c

                # fig = plt.figure(figsize=(6.5, 4))
                # plt.plot(tevlp, trflt0.data, "k", linewidth=1)
                # plt.text(0.02, 1.04, '(c)',transform=plt.gca().transAxes, fontsize=10, va='top',
                #   ha='left')
                # plt.plot(tevlp, evlp, 'k:', markersize= 3.0)
                # plt.vlines(tmin1, 0, Amax/2, "r", linewidth=1)
                # plt.vlines(tmin2, 0, Amax/2, "r", linewidth= 1)
                # plt.plot(tmax,Amax,"ro", markersize= 3.0)
                # plt.plot(tph,0,"bo", markersize= 3.0)
                # plt.gca().yaxis.set_visible(False)
                # plt.gca().spines['top'].set_visible(False)
                # plt.gca().spines['right'].set_visible(False)
                # plt.gca().spines['left'].set_visible(False)
                # plt.xlabel("Time (s)")
                # plt.axis('off')
                # plt.show()


                ## find the time of the largest amplittude
                   
                ct1 = otime + tmin1 - 75   ## set the beginning of the signal window
                ct2 = otime + tmin2 + 90   ## set the ending of the signal window
                #ct1 = otime + tmax - 300 
                #ct2 = otime + tmax + 300
                ## take care of anomalous windows
                if ct1 < wnb:
                    ct1 = wnb + 60.
                    ct2 = wne - 60. 
                if ct2 > wne:
                    ct2 = wne - 60.
                    ct1 = wnb + 60.
                   
                trflt1 = trflt0.slice(starttime = wnb, endtime = ct1)
                trflt2 = trflt0.slice(starttime = ct1, endtime = ct2)
                trflt3 = trflt0.slice(starttime = ct2, endtime = wne)
                trflt1.data = trflt1.data * 0.
                trflt3.data = trflt3.data * 0.
                trflt2 = trflt2.taper(max_percentage = 0.05, side = 'both')
                trfltnew = trflt1 + trflt2 + trflt3
                #if(trfltnew.data,np.ma.masked_array):
                #    trfltnew.data = trfltnew.data.filled()
                
                stXfltnew.append(trfltnew)
                # fig = plt.figure(figsize=(6.5, 4))
                # plt.plot(tevlp, trfltnew.data, color='k', linewidth=1.0)
                 
            for stXfltnew_ in stXfltnew:
                network = stXfltnew_.stats.network
                station = stXfltnew_.stats.station
                print(network, station)
                # figd = plt.figure(figsize=(6.5, 4))
                # stXfltnew_.plot(type = 'section',fig=figd, time_down = True, scale = 2, orientation= "horizontal",color= "k", linewidth= 1.0)
                # plt.text(0.02, 1.04, '(d)',transform=plt.gca().transAxes, fontsize=10, va='top',ha='left')
                # plt.gca().yaxis.set_visible(False)
                # plt.gca().spines['top'].set_visible(False)
                # plt.gca().spines['right'].set_visible(False)
                # plt.gca().spines['left'].set_visible(False)
                # plt.gca().xaxis.grid(False)
                # plt.xlabel("Time (s)")
                # plt.show()
                # plt.close(figd)
                                   
        # stXflt.plot(type = 'section', time_down = True, scale = 2)
        # stXfltnew.plot(type = 'section',time_down = True, scale = 2)
            
## below are rlues to remove low SNR traces based on time and amplitdue of the maximum amplitude phase
## these can be changed based on your own tests on events with different magnitudes and at different distances

            

            if rl == 0:  ## Rayleigh wave
                nkp = 0
                dtt = []
             # the maximum (tmax-tph) depends on distance
                if distf > 90*110*1000:
                    dtmax = 800
                if distf > 60*110*1000 and distf <= 90*110*1000:
                    dtmax = 600
                if distf <= 60*110*1000:
                    dtmax = 400
                
                ## remove noisy traces
                for ll in range(nsta):
                    tph = Rttaa[0][ll]
                    tmax = Rttaa[1][ll]
                    Amax = Rttaa[2][ll]
                    Aavg = Rttaa[3][ll]
                    net = stXfltnew[ll].stats.network
                    sta = stXfltnew[ll].stats.station
                    
                    #print(tph, tmax, Amax, Aavg)
                    
                    if (k <= 7 and Amax/Aavg > 3.5 and ((tmax-tph) > -50 and (tmax -tph) < dtmax)) or \
                       (k >= 8 and k < 15 and Amax/Aavg > 3.5 and ((tmax-tph) > -50 and (tmax -tph) < dtmax -100)) or \
                       (k >= 15 and Amax/Aavg > 3.5 and abs(tmax-tph) < (dtmax - 200)):

                        stXfltkeep.append(stXfltnew[ll])
                        dtt.append(tmax - tph)
                        nkp= nkp + 1
                        
                        if ll < 9:
                            fname = Fltpath + '/' + 'D00' + str(ll+1) + '.' + net + '.' + sta + '.f' + str(freq[k])
                            stXfltnew[ll].write(fname + '.SAC', format = 'SAC')
                            #print (j, fname)
                        if ll >= 9 and ll < 99:
                            fname = Fltpath + '/' + 'D0' + str(ll+1) + '.' + net + '.' + sta + '.f' + str(freq[k])
                            stXfltnew[ll].write(fname + '.SAC', format = 'SAC')
                            #print(j, fname)
                        else:
                            fname = Fltpath + '/' + 'D'+ str(ll+1) + '.' + net + '.' + sta + '.f' + str(freq[k])
                            stXfltnew[ll].write(fname + '.SAC', format = 'SAC')
                        if ll < 9:
                            tfname = Fltpath + '/' + 'D' + str(ll+1) + '.' + net + '.' + sta + '.f' + str(freq[k]) + '.SAC'
                            if os.path.exists(tfname):
                                os.remove(tfname)
              
            if rl == 1: ##Love wave
                nkp = 0
                dtt = []
             # the maximum (tmax-tph) depends on distance
                if distf > 90*110*1000:
                    dtmax = 500
                if distf > 60*110*1000 and distf <= 90*110*1000:
                    dtmax = 300
                if distf <= 60*110*1000:
                    dtmax = 250

                ## remove noisy traces
                for ll in range(nsta):    
                    tph = Lttaa[0][ll]
                    tmax = Lttaa[1][ll]
                    Amax = Lttaa[2][ll]
                    Aavg = Lttaa[3][ll]
                    net = stXfltnew[ll].stats.network
                    sta = stXfltnew[ll].stats.station

                    if (k <= 7 and Amax/Aavg > 3.5 and ((tmax-tph) > -100 and (tmax -tph) < dtmax)) or \
                       (k >= 8 and k < 15 and Amax/Aavg > 3.5 and ((tmax-tph) > -100 and (tmax -tph) < dtmax -100)) or \
                       (k >= 15 and Amax/Aavg > 3.5 and abs(tmax-tph) < (dtmax - 200)):
                        stXfltkeep.append(stXfltnew[ll])
                        dtt.append(tmax - tph)
                        nkp = nkp + 1
                        
                        if ll < 9:
                            fname = Fltpath + '/' + 'D00' + str(ll+1) + '.' + net + '.' + sta + '.f' + str(freq[k])
                            stXfltnew[ll].write(fname + '.SAC', format = 'SAC')
                            #print (j, fname)
                        if ll >= 9 and ll < 99:
                            fname = Fltpath + '/' + 'D0' + str(ll+1) + '.' + net + '.' + sta + '.f' + str(freq[k])
                            stXfltnew[ll].write(fname + '.SAC', format = 'SAC')
                            #print(j, fname)
                        else:
                            fname = Fltpath + '/' + 'D'+ str(ll+1) + '.' + net + '.' + sta + '.f' + str(freq[k])
                            stXfltnew[ll].write(fname + '.SAC', format = 'SAC')
                        if ll < 9:
                            tfname = Fltpath + '/' + 'D' + str(ll+1) + '.' + net + '.' + sta + '.f' + str(freq[k]) + '.SAC'
                            if os.path.exists(tfname):
                                os.remove(tfname)    
            if nkp != 0:
                stXsave = stXfltkeep.copy()
                dttmean = np.mean(dtt)
                dttstd = np.std(dtt)
                for p in range(nkp):
                    if abs(dtt[p]- dttmean) > 2.5 * dttstd:
                #      print((dtt[p]-dttmean), dttstd)
                      stXfltkeep.remove(stXsave[p])
                      sta = stXsave[p].stats.station
                      dfile = glob.glob(Fltpath + '/*' + sta + '*')
                      sdfile = ''.join(dfile)
                      os.remove(sdfile)
                      
            if len(stXfltkeep) != 0:
                print(rl, k, freq[k], len(stXfltkeep))
                fw.writelines ('    ' + str(freq[k]) + '  ' + str(len(stXfltkeep)) + '\n')
                # stXfltkeep.plot(type = 'section', time_down = True, scale = 2, offset_min = distn-5.e4, offset_max = distf + 5.e4)
                # fig1, ax1= plt.subplots(figsize=(6.5, 4))
                # fig2, ax2= plt.subplots(figsize=(6.5, 6))  
                # stXfltnew.plot(type = 'section', time_down = True, fig=fig1, scale = 2, outfile = Fltpath + '/f-' + str(freq[k]) + '.svg')
                # ax = fig1.axes[0]
                # ticks = ax.get_xticks()
                # ax.set_xticklabels([f'{x/1000:.0f}' for x in ticks])
                # ax.set_ylabel('Time (s)', fontsize=10)
                # ax.set_xlabel("Epicentral Distance (Km)", fontsize=10)
                # plt.tight_layout()
                # plt.savefig(Fltpath + '/f-' + str(freq[k]) + '.svg')


                # stXfltkeep.plot(type = 'section', time_down = True,fig=fig2, scale = 2, outfile = Fltpath + '/f-' + str(freq[k]) + 'kp.svg')
                # ax2.set_ylabel('Time (sec)', fontsize=10)
                # ax2.set_xlabel("Epicentral Distance (Km)", fontsize=10)
                # plt.tight_layout()
                # plt.savefig(Fltpath + '/f-' + str(freq[k]) + 'kp.svg')
            # discard the directory if the total number less than 20 traces
            if len(stXfltkeep) < 20: #the value is given based on the size of the network
                shutil.rmtree(Fltpath)
fw.close()

    
