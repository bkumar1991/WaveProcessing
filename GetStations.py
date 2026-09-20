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

stt = UTCDateTime("2017-01-01")
edt = UTCDateTime("2017-01-10")
minlat = 56.0
maxlat = 72.0
minlon = -170.0
maxlon = -132.0


invt = clientfdsn.get_stations(network='*', channel='LH?', starttime=stt, endtime=edt,
                                           minlatitude=minlat,  maxlatitude=maxlat,
                                        minlongitude=minlon, maxlongitude=maxlon, level="response")

invz = invt.select(channel='LH?')

#netwk = invt.networks

with open('Test_Data.txt','w') as file:
    #nsta=0
    total_sta = sum(len(nwk.stations) for nwk in invz.networks)
    file.write(f"{total_sta}\n")   

    for nwk in invz.networks:
        sta = nwk.stations
        #nsta = nsta + len(sta)  
        #print(nwk.code)
        #print('station number : ', len(sta))
        for sta in sta:
            line = "{} {} {} {}\n".format(nwk.code, sta.code, sta.longitude, sta.latitude)
            #print(nwk.code, sta.code, sta.latitude, sta.longitude)
            file.write(line)
#print(nsta)

    
#invt.plot(label=True)

#print(invz)

## set up new projection so that the center can be changed to the studt area  ##
projection = ccrs.PlateCarree(central_longitude=-160.0)
fig = plt.figure(dpi=150)
ax = fig.add_subplot(111, projection=projection)
ax.set_extent((-160, -140, 54, 72))
ax.coastlines()
ax.gridlines()
ax.add_feature(cfeature.OCEAN.with_scale('50m'),
               facecolor='lightblue')
ax.set_title("Station distribution in Alaska")

invt.plot(fig=fig, label= False)


invz.write("Test_Data.xml", format= "STATIONXML")

#with open('test.txt','w') as file:
    #    for i in range(2):
    #        line = "{0:6} {1:4d}\n".format(arr[i][0],arr[i][1])
    #        file.write(line)
file.close()
