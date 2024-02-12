#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 13:56:19 2020

@author: manuelsanchezpiedra
"""

import numpy as np

def lla2ecef(lat, lon, alt):

    lat = np.radians(lat)
    lon = np.radians(lon)
    alt = alt
 
    # WGS84 ellipsoid constants:
    a = 6378137
    e = 8.1819190842622*10**(-2)
 
    # intermediate calculation
    # (prime vertical radius of curvature)
    N = a / np.sqrt(1 - e**2 * (np.sin(lat))**2)
 
    # results:
    X = (N + alt) * np.cos(lat) * np.cos(lon);
    Y = (N + alt) * np.cos(lat) * np.sin(lon);
    Z = ((1 - e**2)  * N + alt) * np.sin(lat);

    
    return(X, Y, Z)
    

def ecef2lla(x,y,z):

    #WGS84 ellipsoid constants:
    a = 6378137;
    e = 8.1819190842622e-2;

    #calculations:
    b   = np.sqrt(a**2*(1-e**2))
    ep  = np.sqrt((a**2-b**2)/b**2)
    p   = np.sqrt(x**2 + y**2)
    th  = np.arctan2(a*z,b*p)
    lon = np.arctan2(y,x)
    lat = np.arctan2((z+(ep**2*b*np.sin(th)**3)),(p-(e**2*a*np.cos(th)**3)))
    N   = a/np.sqrt(1-(e**2*np.sin(lat)**2))
    alt = p/np.cos(lat)-N

    
    lon = np.mod(lon,2*np.pi)

    # correct for numerical instability in altitude near exact poles:
    # (after this correction, error is about 2 millimeters, which is about
    # the same as the numerical precision of the overall function)

    return(lat, lon, alt)
    
    
def geo2topo(X_interp, Y_interp, Z_interp, station, Range_interp):
    
    X_st = station[0]
    Y_st = station[1]
    Z_st = station[2]
    
    lat_st  = station[3]
    long_st = station[4]
    alt     = station[5]
    
    X_ecef = []
    Y_ecef = []
    Z_ecef = []
    
    X_ecef_uni = []
    Y_ecef_uni = []
    Z_ecef_uni = []
    
    north = []
    east = []
    azimuth = []
    
    elevation = []
    
    for i in range(len(X_interp)):
        X_ecef.append(X_interp[i] - X_st)
        Y_ecef.append(Y_interp[i] - Y_st)
        Z_ecef.append(Z_interp[i] - Z_st)
        
        X_ecef_uni.append(X_ecef[i] / Range_interp[i])
        Y_ecef_uni.append(Y_ecef[i] / Range_interp[i])
        Z_ecef_uni.append(Z_ecef[i] / Range_interp[i])

        north.append(-np.cos(long_st) * np.sin(lat_st) * X_ecef_uni[i] - np.sin(long_st) * np.sin(lat_st) * Y_ecef_uni[i] + np.cos(lat_st) * Z_ecef_uni[i])
        east.append(-np.sin(long_st) * X_ecef_uni[i] + np.cos(long_st) * Y_ecef_uni[i])
        azimuth.append(np.cos(long_st) * np.cos(lat_st) * X_ecef_uni[i] + np.sin(long_st) * np.cos(lat_st) * Y_ecef_uni[i] + np.sin(lat_st) * Z_ecef_uni[i])
    
        elevation.append(np.degrees((np.pi)/2 - np.arccos(azimuth[i])))
        
        azimuth[i] = np.arctan(east[i] / north[i])
                          
    for i in range(len(X_interp)):
        if azimuth[i] < 0:
            azimuth[i] = azimuth[i] + 2 * np.pi
        azimuth[i] = np.degrees(azimuth[i])
    
    
    return azimuth, elevation