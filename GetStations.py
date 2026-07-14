## get stations in a rectangel area ##

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from obspy.clients.fdsn import Client
clientfdsn = Client("IRIS")
from obspy import UTCDateTime
from obspy import read
from obspy.core.event import Catalog, read_events
from obspy.clients.iris import Client
clientiris = Client()
import pandas as pd


## read event file
#cat = read_events("events-2007.txt", format="ZMAP")

## set up time and area for available stations

stt = UTCDateTime("2018-01-01")        # Start time 
edt = UTCDateTime("2018-08-31")        # End time
minlat = 56.0
maxlat = 72.0
minlon = -170.0
maxlon = -132.0


invt = clientfdsn.get_stations(network='TA', channel='LH?', starttime=stt, endtime=edt,
                                           minlatitude=minlat,  maxlatitude=maxlat,
                                        minlongitude=minlon, maxlongitude=maxlon, level="response")

invz = invt.select(channel='LHZ')

#netwk = invt.networks

with open('2018_stations.txt','w') as file:
    nsta=0
    for nwk in invz.networks:
        sta = nwk.stations
        nsta = nsta + len(sta)  
        #print(nwk.code)
        #print('station number : ', len(sta))
        for sta in sta:
         line = "{} {} {} {}\n".format(nwk.code, sta.code, sta.longitude, sta.latitude)
#         print(nwk.code, sta.code, sta.latitude, sta.longitude)
         file.write(line)
print(nsta)

invt.plot(label=False)

#print(invz)

## set up new projection so that the center can be changed to the studt area  ##
projection = ccrs.PlateCarree(central_longitude=-120.0)
fig = plt.figure(dpi=150)
ax = fig.add_subplot(111, projection=projection)
ax.set_extent((-125, -113, 28, 50))
ax.coastlines()
ax.gridlines()

#invt.plot(fig=fig)

invz.write("stations-2018.xml", format= "STATIONXML")

## format writing an array content to a file

#with open('test.txt','w') as file:
    #    for i in range(2):
    #        line = "{0:6} {1:4d}\n".format(arr[i][0],arr[i][1])
    #        file.write(line)
file.close()
