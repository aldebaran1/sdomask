#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 12:54:42 2017

@author: Sebastijan Mrak <smrak@bu.edu>
"""

import numpy as np
import datetime
import glob
import scipy.ndimage
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as colors
import matplotlib.pyplot as plt

#EUVDIR = '/home/smrak/Documents/eclipse/MapsSDOdisk/'
EUVDIR = '/home/smrak/Documents/eclipse/MapsSDOdisk300/'


def getEUVMaskCoordinates(latlim=[-89.5,89.5],lonlim=[-180,180],nlat=180,nlon=360):
    xgrid, ygrid = np.mgrid[lonlim[0]:lonlim[1]:nlon*1j, latlim[0]:latlim[1]:nlat*1j]
    return xgrid,ygrid
def getEUVMaskFolder(f,nlat=180,nlon=360):
    xgrid, ygrid = getEUVMaskCoordinates(nlat=nlat, nlon=nlon)
    npts = nlat*nlon
    data = np.fromfile(f,count=npts, dtype=np.float32).reshape((nlat,nlon))
    return xgrid, ygrid, data
def getEUVMask(time,nlat=180,nlon=360):
    """
    I: time in posix
    """
    global EUVDIR
    xgrid, ygrid = getEUVMaskCoordinates(nlat=nlat, nlon=nlon)
    npts = nlat*nlon
    #Import EUV mask files
    flist = np.sort(glob.glob(EUVDIR+'*.bin'))
    Tframe_full = datetime.datetime.utcfromtimestamp(time)
    if int(Tframe_full.strftime('%H')) > 16 and int(Tframe_full.strftime('%H')) < 22:
        # find right filename extension
        TframeHM = Tframe_full.strftime('%H%M')
        flist = np.sort(glob.glob(EUVDIR+'*'+TframeHM+'.bin'))
        # Get Mask
        data = np.fromfile(flist[0],count=npts, dtype=np.float32).reshape((nlat,nlon))
        return xgrid, ygrid, data
    else:
        return 0, 0, 0
    
def plotMap(latlim=[20, 65], lonlim=[-160, -70], center=[39, -86],
            parallels=[20,30,40,50], 
            meridians = [-120,-110, -100, -90, -80,-70],
            epoto=False, totality=True, time=None):
    """
    Plot the map and return handlers of the figure
    """
    (fig,ax) = plt.subplots(1,1,facecolor='w', figsize=(12,8))
    m = Basemap(lat_0=40, lon_0=-95,llcrnrlat=latlim[0],urcrnrlat=latlim[1],
                llcrnrlon=lonlim[0],urcrnrlon=lonlim[1],
                projection='merc')#, resolution='i', ax=ax)
    m.drawmapboundary(fill_color='white')
    
    m.drawcoastlines(color='black')
    m.drawstates(color='black')
    m.drawcountries(color='black')
    m.drawparallels(parallels, labels=[1,0,0,0], linewidth=0.01)
    m.drawmeridians(meridians, labels=[0,0,0,1], linewidth=0.01)
    
    return fig, ax, m
    
def plot(xgrid,ygrid,data, name=''):
#    gradient = 1
#   if isinstance(data, np.ndarray):
    fig, ax, m = plotMap(latlim=[-30,70], lonlim=[-160,0.1], parallels=[-25, -5, 5,25,45,65], meridians=[0, -30,-60,-90,-120,-160])
#    if gradient:
    
    if int(name[2:]) == 60:
        tmp = int(name[:2])+1
        title = datetime.datetime.strptime('2017 08 21 {} {} 0'.format(str(tmp),'00'), '%Y %m %d %H %M %S')
    else:
        title = datetime.datetime.strptime('2017 08 21 {} {} 0'.format(name[:2],name[2:]), '%Y %m %d %H %M %S')
    ax.set_title(title)
    
    x,y = m(xgrid,ygrid)
    
    levels = np.linspace(0.005,0.035,50)
    im = abs(scipy.ndimage.filters.laplace(data))
    m.contourf(x,y,im.T, levels, cmap='jet', norm=colors.PowerNorm(gamma=0.8))
    cbar = m.colorbar(ticks=[0.005, 0.015, 0.025, 0.035])
    cbar.set_label('EUV flux gradient')
    
    levels = np.linspace(0.2,1,30)
    m.contour(x,y,data.T, levels, colors='r', alpha=0.7, linewidths=1)
#    plt.xlim([-150,0])
#    plt.ylim([0,70])
    plt.savefig('euv_mask_raw300/'+name+'.png', dpi=200)
    plt.close(fig)

def plotAll():
    flist = np.sort(glob.glob(EUVDIR+'*.bin'))
    for f in flist:
        xgrid, ygrid, data = getEUVMaskFolder(f)
        plot(xgrid, ygrid, data, name=f[-8:-4])
        
def getDMSPCoord(dmspfn = 'f16_rl172331909_eclipse_only.txt'):
    data = np.loadtxt(dmspfn, skiprows=3)# usecols=None, unpack=False, ndmin=0)[source]
    yydd = int(data[0][0])
    secinday = data[:,1]
    dmsp_dt = [datetime.datetime.strptime(str(yydd)+str(datetime.timedelta(seconds=t)), '%Y%j%H:%M:%S') for t in secinday]
    dmsp_lat = data[:,5]
    dmsp_lon = data[:,6]
    dmsp_Te = data[:,-2]
    
    return dmsp_dt, dmsp_lat, dmsp_lon, dmsp_Te

def plotSingle(td = datetime.datetime(2017,8,21,18,35,0), gradient=0, mask=0, dmsp=0):
    time = td.replace(tzinfo=datetime.timezone.utc).timestamp()
    xgrid, ygrid, data = getEUVMask(time)
    if isinstance(data, np.ndarray):
#        plt.figure()
        fig, ax, m = plotMap(latlim=[0,60], lonlim=[-140,-20])
        title = td.strftime('%Y-%m-%d %H:%M:%S UT')
        ax.set_title(title)
        if gradient:
            x,y = m(xgrid,ygrid)
            levels = np.linspace(-.015,0.015,50)
            levels = np.linspace(0.005,0.035,20)
            im = abs(scipy.ndimage.filters.laplace(data))
            m.contourf(x,y,im.T, levels, cmap='terrain', norm=colors.PowerNorm(gamma=0.8))
            cbar = m.colorbar(ticks=[0.005, 0.015, 0.025, 0.035])
            cbar.set_label('EUV flux gradient')
        if mask:
            x,y = m(xgrid,ygrid)
            levels = np.linspace(0,1,50)
            m.contour(x,y,data.T, levels, colors='r', linewidths=0.5)
        if dmsp:
            dmsp_dt, dmsp_lat, dmsp_lon, dmsp_Te = getDMSPCoord()
            x,y = m(180-dmsp_lon,dmsp_lat)
            m.plot(x,y, 'r', lw=3)
            
            plt.figure()
            plt.plot(dmsp_dt, dmsp_Te, 'b')

#plotAll()

td = datetime.datetime(2017,8,21,19,15,0)
plotSingle(td, gradient=1, mask=1, dmsp=1)

# DMSP


#
#time = td.replace(tzinfo=datetime.timezone.utc).timestamp()
#xgrid, ygrid, data = getEUVMask(time)
#flist = np.sort(glob.glob(EUVDIR+'*.bin'))
#gradient = 1
#for f in flist:
#    xgrid, ygrid, data = getEUVMaskFolder(f)
#    
#    plot(xgrid, ygrid, data, name=f[-8:-4])

#if isinstance(data, np.ndarray):
#    plt.figure()
#    if gradient:
#        levels = np.linspace(-.015,0.015,50)
##    levels = [-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5]
#        im = scipy.ndimage.filters.laplace(data)
#        plt.contour(xgrid,ygrid,im.T, levels, cmap='jet')
#        plt.colorbar()
#        levels = np.linspace(0.2,1,30)
#        plt.contour(xgrid,ygrid,data.T, levels, colors='k', alpha=0.5)
#        plt.xlim([-130,-50])
#        plt.ylim([15,55])
#        plt.savefig('euv_mask_raw/'+)
#    else:
#        levels = np.linspace(0,1,50)
#        plt.contour(xgrid,ygrid,data.T, levels, colors='k')
    
