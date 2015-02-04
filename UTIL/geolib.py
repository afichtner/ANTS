import os
import INPUT.MODELS.wgs84 as wgs84
from math import exp, pi, cos, sin, sqrt
from geographiclib import geodesic,geodesicline

def len_deg_lon(lat):
    
    # This is the length of one degree of longitude 
    # approx. after WGS84, at latitude lat
    # in m
    lat = pi/180*lat
    dlon = (pi*wgs84.a*cos(lat))/180*sqrt((1-wgs84.e_2*sin(lat)**2))
    return round(dlon,2)

def len_deg_lat(lat):
    # This is the length of one degree of latitude 
    # approx. after WGS84, between lat-0.5deg and lat+0.5 deg
    # in m
    lat = pi/180*lat
    dlat = 111132.954 - 559.822 * cos(2*lat) + 1.175*cos(4*lat)
    return round(dlat,2)
    

def area_of_sqdeg(lat):
    # Give back the approx. area of a square degree at latitude lat
    # The sphericity of the Earth is not (yet) taken into account
    # This is a piecewise flat earth
    # in m^2
    l_lat = len_deg_lat(lat)
    l_lon = len_deg_lon(lat)
    
    return round(l_lat*l_lon,2)
        



def get_gcsegs(lat1,lon1,lat2,lon2,num,atten=False,\
                sta_dist=None,freq=None,Q=None,v=None):
    """
    Obtain coordinates of locations that separate the great circle b/w lat1lon1 
    and lat2lon2 into num equal segments
    
    """
    
    p=geodesic.Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2) 
    l=geodesic.Geodesic.WGS84.Line(p['lat1'],p['lon1'],p['azi1'])
    newcoords = list()

    if atten == False:
        for i in range(num+1):
            b=l.Position(i*p['s12']/(num))
            newcoords.append((b['lat2'],b['lon2']))
    
    else: 
       #determine omega
       w = 2 * pi * freq
       for i in range(num+1):
           #determine x
           x = sta_dist / 2. + i * p['s12'] / (num)
           # determine K on the ray
           k = exp(-w * x / (v * Q))

           b=l.Position(i*p['s12']/(num))
           newcoords.append((b['lat2'],b['lon2'],k)) 
        
    return newcoords
    
#def get_dec_segs(lat1,lon1,lat2,lon2,num,sta_dist,freq,v,Q):
    
    #p=geodesic.Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2); 
    #l=geodesic.Geodesic.WGS84.Line(p['lat1'],p['lon1'],p['azi1'])
    #coord_k = list()
    
    #return coord_k    
    

def get_midpoint(lat1,lon1,lat2,lon2):
    
    """
    Obtain the coordinates of the point which is halfway between 
    point 1 (lat1, lon1) and point 2 (lat2, lon2) on great circle
    """
    
    p=geodesic.Geodesic.WGS84.Inverse(lat1,lon1,lat2,lon2)
    l=geodesic.Geodesic.WGS84.Line(p['lat1'],p['lon1'],p['azi1'])
    b=l.Position(0.5*p['s12'])

    return (b['lat2'],b['lon2'])
    

def get_antipode(lat,lon):
    
    if lon <= 0.:
        lon_a = 180. + lon
    else:
        lon_a = lon - 180.
    
    lat_a = -1. * lat
    
    return(lat_a,lon_a)