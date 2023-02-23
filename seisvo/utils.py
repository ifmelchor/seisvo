

import math
import datetime as dt
from obspy import UTCDateTime
from obspy.clients.fdsn.client import Client
from obspy.taup import TauPyModel
from geopy.distance import geodesic
import matplotlib.ticker as mtick
import matplotlib.dates as mdates

def distance_in_degree(coord1, coord2):
    # Convert latitude and longitude to
    # spherical coordinates in radians.
    (lat1, long1) = coord1
    (lat2, long2) = coord2

    degrees_to_radians = math.pi / 180.0

    # phi = 90 - latitude
    phi1 = (90.0 - lat1) * degrees_to_radians
    phi2 = (90.0 - lat2) * degrees_to_radians

    # theta = longitude
    theta1 = long1 * degrees_to_radians
    theta2 = long2 * degrees_to_radians

    # Compute spherical distance from spherical coordinates.
    cos = (math.sin(phi1) * math.sin(phi2) * math.cos(theta1 - theta2) +
           math.cos(phi1) * math.cos(phi2))
    arc = math.acos(cos)

    # Multiply arc by the radius of the earth in km
    return arc


def get_catalog(startime, endtime, receiver):
    
    USGSclient = Client("USGS")
    Model_IASP91 = TauPyModel(model="iasp91")
    
    cat = USGSclient.get_events(starttime=UTCDateTime(startime), endtime=UTCDateTime(endtime), minmagnitude=1, orderby='time')
    
    if len(cat) == 0:
        return None
    
    else:
        quakes = []
        for key in cat:
            lat_s = key.origins[0]['latitude']
            lon_s = key.origins[0]['longitude']
            source = (lat_s, lon_s)
            time = key.origins[0]['time'].datetime
            distance_to_receiver = geodesic(source, receiver).km
            dist_degree = distance_in_degree(source, receiver)*180/np.pi
            depth = float(key.origins[0]['depth'])/1000
            
            if key.magnitudes:
                magnitude = '%.1f %s' % (key.magnitudes[0]['mag'], key.magnitudes[0]['magnitude_type'])
            
            else:
                magnitude = None
            
            try:
                arrivals = Model_IASP91.get_travel_times(
                    source_depth_in_km=depth,
                    phase_list=['ttbasic'],
                    distance_in_degree=dist_degree
                    )
                first_arrival = float(arrivals[0].time)
                time_P = time + dt.timedelta(seconds=first_arrival)
            except:
                time_P = None

            quakes.append([time, source, distance_to_receiver, depth, magnitude, dist_degree, time_P])
    
        return quakes


def get_times_bounds(starttime_bound, endtime_bound, st, et):
    
    in_in = st >= starttime_bound and et <= endtime_bound # start and end in
    st_in = endtime_bound > st > starttime_bound and et > endtime_bound # start in, end out
    et_in = starttime_bound > st and starttime_bound < et < endtime_bound # start out, end in
    out = st < starttime_bound and endtime < et # start and end out

    if in_in:
        return (st, et)

    if st_in:
        return (st, endtime_bound)

    if et_in:
        return (starttime_bound, et)

    if out:
        return(starttime_bound, endtime_bound)


def get_time_format(datetime, day_interval):
    if datetime:
        if day_interval <= 1:
            major_locator = mdates.HourLocator(interval=1)
            major_formatt = mdates.DateFormatter('%d %b\n%H:%M')
            minor_locator = mdates.MinuteLocator(byminute=[15, 30, 45])
            minor_formatt = mtick.NullFormatter()

        elif day_interval <= 10:
            major_locator = mdates.DayLocator(interval=1)
            major_formatt = mdates.DateFormatter('%d %b %H:%M')
            minor_locator = mdates.HourLocator(byhour=[6, 12, 18, 24])
            minor_formatt = mtick.NullFormatter()

        elif 45 >= day_interval > 10 :
            major_locator = mdates.DayLocator(interval=7)
            major_formatt = mdates.DateFormatter('%d')
            minor_locator = mdates.DayLocator(interval=1)
            minor_formatt = mtick.NullFormatter()

        else:
            major_locator = mdates.WeekdayLocator(interval=2)
            major_formatt = mdates.DateFormatter('%d-%m')
            minor_locator = mdates.DayLocator(interval=7)
            minor_formatt = mtick.NullFormatter()
    
    else:
        major_locator = mtick.LinearLocator(10)
        major_formatt = mtick.FormatStrFormatter('%i')
        minor_locator = mtick.AutoMinorLocator(2)
        minor_formatt = None
    
    return (major_locator, major_formatt), (minor_locator, minor_formatt)
