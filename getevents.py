##select events for a distance range of 30-120 degree ##
## remove events that are too close (0.5 degree) ##

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from obspy.clients.fdsn import Client
clientfdsn = Client("IRIS")
from obspy import UTCDateTime
from obspy import read
from obspy.core.event import Catalog
from obspy.clients.iris import Client
clientiris = Client()
from obspy.clients.fdsn.mass_downloader import RectangularDomain, \
    Restrictions, MassDownloader
import pandas as pd


stlat = 64.0
stlon = -152.0


## get_events
StEvTime = UTCDateTime("2017-08-01")
EndEvTime = UTCDateTime("2018-12-31")
minlat = 56.0
maxlat = 72.0
minlon = -170.0
maxlon = -132.0
cat = clientfdsn.get_events(starttime=StEvTime, endtime=EndEvTime,
                        minmagnitude=6.0, catalog=None)
cat1 = Catalog()
cat.plot()

#print(cat)
#print(cat.__str__(print_all=True))




for  i in range(0, len(cat)):
    otime = cat[i].origins[0].time
    evlat = cat[i].origins[0].latitude
    evlon = cat[i].origins[0].longitude
    evdp = cat[i].origins[0].depth/1000
    evmag = cat[i].magnitudes[0].mag
    distaz = clientiris.distaz(stlat, stlon, evlat, evlon)
    dist = distaz['distance']
    
    if dist > 30 and dist < 120 and evdp < 150:
        #print(dist)
        ## check the distance between two events in the sequence
        ## and do not add it to the catalog if they are too close
        print(dist, evdp)
       # ncount = 0
       # for j in range(0,len(cat1)):
       #     evlat0 = cat1[j].origins[0].latitude
       #     evlon0 = cat1[j].origins[0].longitude
       #     distaz = clientiris.distaz(evlat0, evlon0, evlat, evlon)
       #     dist0 = distaz['distance']
       #     #print(dist0)
       #     if dist0 < 0.5:
       #         ncount = ncount + 1
       # if ncount == 0:
        cat1.append(cat[i])
        
        
print(cat1.__str__(print_all=True))


cat1.write("AE2017-18.xml", format="QUAKEML")
cat1.write("AE2017-18.txt", format="ZMAP")

## set up new projection so that the center can be changed to the studt area  ##
projection = ccrs.PlateCarree(central_longitude=-120.0)
fig = plt.figure(dpi=150)
ax = fig.add_subplot(111, projection=projection)
ax.set_extent((-180, 180, -90, 90))
ax.coastlines()
ax.gridlines()

cat1.plot(fig=fig)



