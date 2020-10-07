#!/usr/bin/env python
import numpy as np
from mpas_tools.ocean import build_spherical_mesh
from mpas_tools.mesh.creation.util import lonlat2xyz
from numpy import radians, cos, sin, arcsin, sqrt
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')


def cellWidthVsLatLon():
    """
    Create cell width array for this mesh on a regular latitude-longitude grid.

    Returns
    -------
    cellWidth : ndarray
        m x n array of cell width in km

    lon : ndarray
        longitude in degrees (length n and between -180 and 180)

    lat : ndarray
        longitude in degrees (length m and between -90 and 90)
    """
    createPlots = True
    ddeg = 10.0 # computed grid resolution, in degrees
    latCenter = 45.0 # center point in degrees
    lonCenter = 0.0 # center point in degrees

# Giacomo, after we get this working with an earth-sized sphere, we need to
# change these flags to be a unit sphere, or reduce it to a unit sphere at the
# end of this routine
    courseResolution = 240 # km
    refinementFactor = 2
    fineRadius = 1000 # km
    transitionWidth = 100 # km
    earthRadius = 6.371e3 # radius in km
    
    lat = np.arange(-90, 90.01, ddeg)
    lon = np.arange(-180, 180.01, ddeg)
    latCenter = np.deg2rad(latCenter)
    lonCenter = np.deg2rad(lonCenter)

    # create meshed 2D grid in radians
    lonGrid, latGrid = np.meshgrid(np.deg2rad(lon), np.deg2rad(lat))

    # Halversine formula for distance
    # better to compute it in code directly, but example functions are below
# Giacomo, I just got this formula from this post:
# https://stackoverflow.com/a/51722117/12656503
# please double check it against an original formula.
    distance = earthRadius * np.sin((latGrid - latCenter)/2)**2 \
        + np.cos(latCenter)*np.cos(latGrid) * np.sin((lonGrid - lonCenter)/2)**2
     
    tanhDistance = np.tanh(distance)
    # Next line is just a dummy line, fill in with actual cell width.
    cellWidth = courseResolution*np.ones(np.shape(distance))

    if createPlots:
        varList = ['latGrid','lonGrid','distance','tanhDistance','cellWidth']
        fig = plt.gcf()
        plt.clf()
        fig.set_size_inches(20.0, 20.0)
        iPlt = 1
        for varName in varList:
            plt.subplot(3,2,iPlt)
            plt.imshow(vars()[varName])
            iPlt += 1
            plt.title(varName)
            plt.xlabel('lon index')
            plt.ylabel('lat index')
            plt.colorbar()
        plt.savefig('cellWidth.png')
    return cellWidth, lon, lat


# code examples from
# https://stackoverflow.com/questions/4913349/haversine-formula-in-python-bearing-and-distance-between-two-gps-points
# this fuction can be removed later if not used.
def distance(s_lat, s_lng, e_lat, e_lng):

   # approximate radius of earth in km
   R = 6373.0

   s_lat = s_lat*np.pi/180.0                      
   s_lng = np.deg2rad(s_lng)     
   e_lat = np.deg2rad(e_lat)                       
   e_lng = np.deg2rad(e_lng)  

   d = np.sin((e_lat - s_lat)/2)**2 + np.cos(s_lat)*np.cos(e_lat) * np.sin((e_lng - s_lng)/2)**2

   return 2 * R * np.arcsin(np.sqrt(d))

# code examples from
# https://stackoverflow.com/questions/4913349/haversine-formula-in-python-bearing-and-distance-between-two-gps-points
# this fuction can be removed later if not used.
def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """

    #Convert decimal degrees to Radians:
    lon1 = np.radians(lon1.values)
    lat1 = np.radians(lat1.values)
    lon2 = np.radians(lon2.values)
    lat2 = np.radians(lat2.values)

    #Implementing Haversine Formula: 
    dlon = np.subtract(lon2, lon1)
    dlat = np.subtract(lat2, lat1)

    a = np.add(np.power(np.sin(np.divide(dlat, 2)), 2),  
            np.multiply(np.cos(lat1), 
            np.multiply(np.cos(lat2), 
            np.power(np.sin(np.divide(dlon, 2)), 2))))
    c = np.multiply(2, np.arcsin(np.sqrt(a)))
    r = 6371

    return c*r

def main():
    cellWidth, lon, lat = cellWidthVsLatLon()
    build_spherical_mesh(cellWidth, lon, lat, out_filename='base_mesh.nc')


if __name__ == '__main__':
    main()
