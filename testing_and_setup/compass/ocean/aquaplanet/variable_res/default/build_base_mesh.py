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
    ddeg = 1.0 # computed grid resolution, in degrees
    latCenter = 90.0 # center point in degrees
    lonCenter = 135.0 # center point in degrees

    # Giacomo, after we get this working with an earth-sized sphere, we need to
    # change these flags to be a unit sphere, or reduce it to a unit sphere at the
    # end of this routine
    coarseResolution = 240 # km
    fineResolution = 140 # km
    refinementFactor = 2
    fineRadius = 5000 # km
    transitionWidth = 100 # km
    earthRadius = 6.371e3 # radius in km
    
    lat = np.arange(-90, 90.01, ddeg)
    lon = np.arange(-180, 180.01, ddeg)
    latCenter = np.deg2rad(latCenter)
    lonCenter = np.deg2rad(lonCenter)

    # create meshed 2D grid in radians
    lonGrid, latGrid = np.meshgrid(np.deg2rad(lon), np.deg2rad(lat))

    # Halversine formula for distance
    distance =  np.sin((latGrid - latCenter)*0.5)**2 + np.cos(latCenter)*np.cos(latGrid) * np.sin((lonGrid - lonCenter)*0.5)**2
    distance = 2.0 * earthRadius * np.arctan2(np.sqrt(distance),np.sqrt(1.0 - distance)) 

    tanhDistance = np.tanh((fineRadius-distance)/transitionWidth)
     
    cellWidth = coarseResolution + (fineResolution - coarseResolution) * tanhDistance

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

def main():
    cellWidth, lon, lat = cellWidthVsLatLon()
    build_spherical_mesh(cellWidth, lon, lat, out_filename='base_mesh.nc')


if __name__ == '__main__':
    main()
