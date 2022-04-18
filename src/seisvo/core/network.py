#!/usr/bin/python3
# coding=utf-8

from seisvo import __seisvo__, Station, Array
from .imports import get_network

class Network(object):
    def __init__(self, net_code):
        self.info = get_network(net_code)
        self.station = [Station(x) for x in self.info.stations]

    def __str__(self):
        return self.info.__str__()

    def __len__(self):
        return len(self.station)

    def __getitem__(self, item):
        return self.station[item]
    
    def get_array(self, list_sta_code, exclude_component=[]):
        return True

    def get_sta(self, sta_code, loc=''):
        """
        Get a station object
        """
        for sta in self.station:
            if sta_code == sta.info.code and loc ==sta.info.loc:
                return sta
        return None


    def plot_map(self, zoom_scale=1, arcgis_map='World_Shaded_Relief', epsg=4326, pixel=1500, dpi=100, save=False):
        """
        Plot a map for geograpich network
        :param zoom_scale: how much to zoom from coordinates (in degrees)
        :param map: http://server.arcgisonline.com/arcgis/rest/services
        """
        
        from seisvo.utils.maps import get_map

        # desired coordinates
        coord = [self.info.latOrig, self.info.lonOrig]
        title = "Network: %s" % self.info.code
        
        fig, ax = get_map(coord, arcgis_map, zoom_scale, epsg, pixel=pixel, dpi=dpi, title=title)

        # draw station
        for station in self.station:
            (lat, lon) = station.get_latlon(degree=True)
            ax.scatter(lon, lat, marker='^', label=station.info.id, transform=proj)
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        if save:
            fig.savefig('map.png', format='png', dpi=500)
        else:    
            fig.show()





