import os
from geographiclib import geodesic,geodesicline

def get_gcsegs(lat1,lon1,lat2,lon2,num):
    """
    Obtain coordinates of locations that separate the great circle b/w lat1lon1 
    and lat2lon2 into num equal segments
    
    """
    
    p=geodesic.Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2); 
    l=geodesic.Geodesic.WGS84.Line(p['lat1'],p['lon1'],p['azi1'])
    newcoords = list()
    
    for i in range(num+1):
        b=l.Position(i*p['s12']/(num))
        newcoords.append((b['lat2'],b['lon2']))
    return newcoords
    

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