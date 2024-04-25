

import math
import numpy as np
import os
import shutil
import time as pytime
import datetime as dt
import subprocess as sp
from obspy import UTCDateTime
from obspy.clients.fdsn.client import Client
from obspy.taup import TauPyModel
from geopy.distance import geodesic
import multiprocessing as mp
import simplekml

nCPU = mp.cpu_count()

COMMAND_LIST = ["LET", "H71", "LST", "KPR", "TOP", "REP", "ZTR", "POS", "MIN", "JUN", "DIS", "RMS", "WET", "ERF", "STA", "PHS"]

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


def in_interval(st, et, time_interval):
    cond1 = st >= time_interval[0] and et <= time_interval[1] # start and end in
    cond2 = time_interval[1] > st > time_interval[0] and et > time_interval[1] # start in, end out
    cond3 = time_interval[0] > st and time_interval[0] < et < time_interval[1] # end out, start in
    cond4 = st < time_interval[0] and time_interval[1] < et # start and end out
    return cond1 or cond2 or cond3 or cond4


def prt2dict(prtfile):
    fout  = open(prtfile, "r")
    lines = fout.readlines()
    fout.close()

    # search the start line
    def get_init_line(lnlist):
        for n, line in enumerate(lnlist):
            if line:
                x = list(filter(lambda x : x not in ("", "\n"), line.split(" ")))
                if x:
                    if x[0] == "STA" and x[1] == "NET":
                        return n

    ni = get_init_line(lines) + 1
    sta_code = None
    phsdict  = {}
    for ln in lines[ni:]:
        x = list(filter(lambda x : x not in ("", "\n"), ln.split(" ")))
        if x[0] not in ("S", "P"):
            if sta_code:
                if wave == "P":
                    resP = res
                    resS = res2
                else:
                    resP = res2
                    resS = res

                phsdict[sta_code] = {
                    "epi":epi_dist,
                    "resP":resP,
                    "resS":resS,
                }

            sta_code = x[0]
            epi_dist = float(x[1])
            res      = float(ln[61:66])
            res2     = None
            wave     = x[4]
            wave2    = None
            if "P" in wave:
                wave = "P"
            else:
                wave = "S"
        else:
            res2 = float(ln[61:66])
            if wave == "P":
                wave2 = "S"
            else:
                wave2 = "P"
    
    return phsdict


def sum2dict(sumfile):
    fout  = open(sumfile, "r")
    line = fout.readlines()
    fout.close()

    try:
        x = list(filter(lambda x : x not in ("", "\n"), line[0].split(" ")))
    except:
        return None
    
    # replace year 19 to 20
    xtime = list(x[0]+x[1]+x[2])
    xtime[0] = '2'
    xtime[1] = '0'
    xtime = ''.join(xtime)
    xtime = dt.datetime.strptime(xtime, "%Y%m%d%H%M%S.%f")

    # lat
    latdeg = float(line[0][19:22]) #deg
    sign   = line[0][22]    #sign
    if sign == "S":
        sign = -1
    else:
        sign = 1
    latmin = float(line[0][23:29])/60 #minute decimal
    lat = sign*(latdeg + latmin)

    #lon
    londeg = float(line[0][29:32]) #deg
    sign = line[0][32]    #sign
    if sign =="W":
        sign = -1
    else:
        sign = 1
    lonmin = float(line[0][33:39])/60 #minute decimal
    lon = sign*(londeg + lonmin)

    # depth
    depth = float(line[0][39:45])

    # azimuth gap
    az_gap = float(line[0][56:59])

    rms = float(line[0][65:69])
    horizontal_error = float(line[0][70:74])
    vertical_error   = float(line[0][75:79])
    quality_code     = line[0][80]

    hypdict = {
        "origin_time":xtime,
        "lat":lat,
        "lon":lon,
        "depth":depth,
        "az_gap":az_gap,
        "rms":rms,
        "h_error":horizontal_error,
        "v_error":vertical_error,
        "quality":quality_code
    }

    return hypdict


class HILoc(object):
    """
        Object class for locate earthquakes using Hypoinverse.
        vel_model :: a string or a list of string that indicates the velocity models input files
        model_ref :: int with the highest altitude of seismic stations
    """

    def __init__(self, vel_model, model_ref, **kwargs):

        self._set_attr(**kwargs)
        self.mod_ref = model_ref
        self.mod     = []
        self.hypfile = None
        self.sta     = None
        self.phs     = None
        self._set_model(vel_model)


    def _set_attr(self, **kwargs):
        ## FILE FORMATS AND RELATED CONTROLS
        # Use new, longer SEED station codes
        self.let = kwargs.get("let", [4, 0, 0, 0, 0])
        # Use hypo71 summary format
        self.h71 = kwargs.get("h71", [2, 1, 2])
        
        ## PRINTED OUTPUT FORMAT
        # Station list or models in printfile
        self.lst = kwargs.get("lst",[2, 0, 1])
        # Medium print output each event
        self.kpr = kwargs.get("kpr", 3)
        # No page ejects
        self.top = kwargs.get("top","T")
        # Log events to terminal; don't print unweighted stations
        self.rep = kwargs.get("rep", ["T", "F"])
        ##

        ## TRIAL DEPTH, VELOCITY RATIO & ERRORS
        # Trial depth for locate
        self.ztr = kwargs.get("ztr", [5.0, "F"])
        # Velocity ratio Vp/Vs
        self.pos = kwargs.get("pos", 1.78)
        # Min number of phases to locate
        self.min = kwargs.get("min", 4)
        ##

        ## WEIGHTING OF ARRIVAL TIMES
        # Force location of junk events
        self.jun = kwargs.get("jun", "F")
        # Main Distance weighting
        self.dis = kwargs.get("dis", [4, 50, 1, 3])
        # Residual weighting
        self.rms = kwargs.get("rms", [4, 0.16, 1.5, 3])
        # Weights for P&S weight codes 0-3
        self.wet = kwargs.get("wet", [1, 0.75, 0.2, 0.1])
        ##

        # Send error messages to terminal
        self.erf = kwargs.get("erf","T")

    
    def _check_file(self, file):
        ok = True

        if isinstance(file, str):
            if not os.path.isfile(file):
                ok = False
        
        elif isinstance(file, (list, tuple)):
            for f in file:
                if not os.path.isfile(f):
                    ok = False
        
        else:
            ok = False
        
        return ok


    def _set_model(self, model):
        if isinstance(model, (list, tuple)):
            for mod in model:
                if self._check_file(mod):
                    self.mod.append(mod)

        else:
            if self._check_file(model):
                self.mod = [model]


    def _check_input(self):
        if self.sta and self.phs and self.mod:
            return True
        else:
            return False
        

    def _check_commands(self, command_list):
        ok = True

        for com in command_list:
            if com not in COMMAND_LIST:
                ok = False

        return ok


    def _get_lines(self):
        lines = "* Hypoinverse startup file [SEISVO]\n* automatic input"

        for com in COMMAND_LIST:
            com_val = getattr(self, com.lower())
            line    = f"\n{com} "
            
            if isinstance(com_val, (list, tuple)):
                for cv in com_val:
                    if isinstance(cv, (str, int)):
                        line += f"{cv} "                        
                    else:
                        line += f"{cv:.2f} "

            if isinstance(com_val, (str, int)):

                if com in ("PHS", "STA"):
                    line += f"'{com_val}' "
                else:
                    line += f"{com_val} "
                
            if isinstance(com_val, float):
                line += f"{com_val:.2f} "
            
            lines += line
        
        lines += "\nFIL " # check phase file

        # add models
        for n, mod in enumerate(self.mod):
            lines += f"\nCRE {n+1} '{mod}' {self.mod_ref:.1f} T"
        
        return lines


    def load_stafile(self, sta):
        if self._check_file(sta):
            self.sta = sta
        
        else:
            print("error reading station file")


    def load_phsfile(self, phs):
        if self._check_file(phs):
            self.phs = phs
        
        else:
            print("error reading phase file")
            self.phs = None


    def _write(self, hypfile, outfile):
        if not self._check_input():
            print(" error with input files")
            return False

        # create directory for config file
        hypdir = os.path.dirname(hypfile)
        if not os.path.isdir(hypdir):
            os.makedirs(hypdir)

        # create directory for output files
        outdir = os.path.dirname(outfile)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)

        # init config lines
        try:
            lines = self._get_lines()
        except:
            return False

        # add output lines
        lines += "\n* output"
        lines += f"\nPRT '{outfile}.prt'"
        lines += f"\nSUM '{outfile}.sum'"
        lines += f"\nARC '{outfile}.arc'"
        lines += "\n* run " # locate the event
        lines += "\nLOC " # locate the event
        lines += "\nSTO " # stop the program

        fout = open(hypfile, "w")
        fout.write(lines)
        fout.close()
        
        return True
    

    def run(self, idname, hypfile, outfile):
        ok = self._write(hypfile, outfile)
        if ok:
            # execute hypoinverse
            p = sp.Popen(['hyp1.40'], stdin=sp.PIPE)
            s = "@{}".format(hypfile)
            p.communicate(s.encode())

            # wait for a second before reading output SUM file
            pytime.sleep(0.5)
            loc = sum2dict(outfile+".sum")

            # make kml file from sum file
            kml = get_kml(idname, loc)
            kml.save(outfile+".kml")

            return True
        
        return False
        
            
def dd2ddm(dd, which):
    """
    Create the string used fo Hyp71 file format
    
    Input degree decimal (DD)
    Output degree decimal minutes (DDM) 

    """

    assert which in ("lat", "lon")
    
    deg = int(dd)

    if deg < 0:
        if which == "lat":
            let = "S"
        else:
            let = "W"
    else:
        if which == "lon":
            let = "E"
        else:
            let = "N"

    dmin = int((dd - deg) * 60)
    dsec = round(((dd - deg) - (dmin / 60)) * 3600, 2)/60
    dmin += dsec
    deg  = abs(deg)
    dmin = abs(dmin)
    
    return f"{deg}{dmin:5.2f}{let}"


def get_kml(title, sumdict):

    def create_circle(lat, lon, erh):
        plts = []
        for th in np.pi*np.arange(0, 360, 1)/180:
            x = erh * np.cos(th) / 111.132
            y = erh * np.sin(th) / 111.132
            plts += [(lon + x, lat + y)]
        return(plts)

    kml = simplekml.Kml(name=title)
    loc = sumdict
                
    pnt = kml.newpoint(name="")
    pnt.coords = [(loc["lon"], loc["lat"])]
    pnt.style.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png'
    pnt.style.iconstyle.color = simplekml.Color.red
    pnt.description = f"{loc['origin_time'].strftime('%d %b %Y %H:%M:%S')}\nDepth {loc['depth']} km\nRMS {loc['rms']}\nH error {loc['h_error']} km\nV error {loc['v_error']} km\nQuality {loc['quality']}"

    pol = kml.newpolygon(name="", description='LOC Error', 
            outerboundaryis=create_circle(loc["lat"], loc["lon"], loc['h_error']))
    pol.style.linestyle.color = simplekml.Color.black
    pol.style.polystyle.outline = 1
    pol.style.polystyle.fill = 0
    pol.linestyle.width = 2

    return kml





    




    