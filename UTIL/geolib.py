import os
import INPUT.MODELS.wgs84 as wgs84
from math import exp, pi, cos, sin, sqrt
from geographiclib import geodesic,geodesicline
# Can approximate pieces of Earth surface area by spherical earth surface element or by square lat-lon boxes on the ellipsoid. Quite similar results. 


def approx_surf_el(dlat,dlon,lat):
    # determine colatitude:
    colat = abs(lat-90.)
    # Radians
    colat = colat/180.*pi
    dlat = dlat/180.*pi
    dlon = dlon/180.*pi
    
    return(dlat*dlon*sin(colat))

def area_surfel(dlat,dlon,lat,r):
    surf_el = approx_surf_el(dlat,dlon,lat)
    
    return(r**2*surf_el)

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
    if l_lat*l_lon == 0:
        area = 0.000001
    else:
        area = round(l_lat*l_lon,2)
    return area
        



def get_gcsegs(lat1,lon1,lat2,lon2,num,num_max=None,atten=False,\
                sta_dist=None,freq=None,Q=None,v=None,line_kern=False):
    """
    Obtain coordinates of locations that separate the great circle b/w lat1lon1 
    and lat2lon2 into num equal segments
    lat1, lon1, lat2, lon2: Coordinates of the two points that define the great circle ends
    num: number of segments the gc is split into
    atten: Use attenuation or not
    sta_dist, freq, Q, v: interstation distance, central frequency, quality fac-
    tor and group velocity; all needed only if atten==True
    
    """
    # Distance between the station and the antipode of the station station midpoint
    p=geodesic.Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2) 
    l=geodesic.Geodesic.WGS84.Line(p['lat1'],p['lon1'],p['azi1'])
    newcoords = list()
    
    # With num_max, can limit the number of segments we want to get at all.
    if num_max==None:
        num_max=num+1
    
    
    if atten == False:
        for i in range(num_max):
            b=l.Position(i*p['s12']/(num))
            newcoords.append((b['lat2'],b['lon2']))
            
    else: 
       #determine omega
       w = 2 * pi * freq
       for i in range(num_max):
           #determine x; the path obtained with geodesic is split in equal length segments.
           x = sta_dist / 2. + i * p['s12']/num
           #x = i * p['s12'] / (num)

           # determine K on the station-station line:
           if line_kern == True:
               k = exp(-w * x / (v * Q))/sqrt(x**2+sta_dist**2/4)
           # The following kernel is integrated along the y-direction
           else:
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