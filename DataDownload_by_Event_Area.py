# This script is for mass-downloading data for a given rectangle area at selected events
# The events can be pulled from data centers for the time period with selected magnitude and epicentral distance
# For surface wave studies, the massdownload allows choosing interstation distance to avoid very close stations 
import matplotlib.pyplot as plt
import sys,os
from obspy.clients.fdsn import Client
clientfdsn = Client("IRIS")
from obspy import UTCDateTime
from obspy import read
from obspy.core.event import Catalog, read_events
from obspy.clients.iris import Client
clientiris = Client()
from obspy.clients.fdsn.mass_downloader import RectangularDomain, \
    Restrictions, MassDownloader
import pandas as pd

## get_events
StEvTime = UTCDateTime("2018-01-01")
EndEvTime = UTCDateTime("2018-08-31")
minlat = 56.0
maxlat = 72.0
minlon = -170.0
maxlon = -132.0

stlat = 64.0
stlon = -150.0

## read event file
cat = read_events("events-2018.txt", format="ZMAP")


print(cat)
print(cat[0].origins[0].time, StEvTime, EndEvTime)
print('number of events : ', len(cat))

# root path for waveform downloaded by year

datapath = './' + str(StEvTime.year) + '/'

print(datapath)
if not os.path.exists(datapath):
    os.makedirs(datapath)

# Set up the name of catalog for downloaded events
data_catalog = 'catalog-' + str(StEvTime.month) + '.csv' 
# Create an empty dataframe for downloaded events 
data_cat_file = pd.DataFrame(columns = ['DateTime',  'Latitude', 
                'Longitude', 'Depth', 'Magnitude', 'Distance'])


Nev_downloaded = 0

#for  i in range(0, len(cat)):
    #otime = cat[i].origins[0].time
    #evlat = cat[i].origins[0].latitude
    #evlon = cat[i].origins[0].longitude
for  eq in cat.events:
    otime = eq.origins[0].time
    evlat = eq.origins[0].latitude
    evlon = eq.origins[0].longitude
    stt = UTCDateTime(otime)
    if otime >= StEvTime and otime <= EndEvTime:
        #print (otime)
        distaz = clientiris.distaz(stlat, stlon, evlat, evlon)
        dist = distaz['distance']
        #print(dist)
    #evdp = eq.origins[0].depth/1000
    #evmag = eq.magnitudes[0].mag

    # The window length can be different for events with different distances
    
        if dist <= 60:
            wndlen = 3000
        if dist > 60 and dist <= 80:
            wndlen = 3600
        if dist > 80 and dist <= 100:
            wndlen = 4200
        if dist > 100:
            wndlen = 4800
        #print (wndlen)
        # generate path name by event            
        strotime = otime.isoformat()
        strotime = strotime[0:10] + '-' + strotime[11:13] + '-' + strotime[14:16]
        print(strotime)
        
        #download data for one event for all stations within a rectangle range 
        # Rectangular domain containing NWUS
                            
        # populate dataframe with filtered by distance events for csv catalog file
        data_cat_file.loc[len(data_cat_file.index)] = [
             #eq.origins[0].time.date.strftime("%Y%m%d"),
             eq.origins[0].time,
             #eq.origins[0].time.strftime("%H:%M:%S"),
             eq.origins[0].latitude,
             eq.origins[0].longitude,
             eq.origins[0].depth/1000,
             eq.magnitudes[0].mag,
             dist]  

        domain = RectangularDomain(minlatitude=minlat, maxlatitude=maxlat,
                   minlongitude=minlon, maxlongitude=maxlon)
        

## wndlen is used to calculate the end time of the event
        restrictions = Restrictions(
        starttime = stt, 
        endtime = stt + wndlen,
        reject_channels_with_gaps=True,
        minimum_length=0.95,
        minimum_interstation_distance_in_m=30E3,
        channel_priorities=["LH[ZNE]"],
        location_priorities=["", "00", "10"])

        mseed_storage = datapath + strotime + '/waveforms'
        stationxml_storage = datapath + strotime + '/stations'

        mdl = MassDownloader(providers=["IRIS"])
        mdl.download(domain, restrictions, mseed_storage,
                 stationxml_storage)


    data_cat_file.to_csv(datapath + data_catalog, index=False)
    ## index = false means the default row indexes are not included in the csv file
                 

            



