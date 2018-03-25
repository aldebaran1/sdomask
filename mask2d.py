# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 13:21:32 2018

@author: smrak
"""
from os import path
from glob import glob
from h5py import File
from datetime import datetime
from numpy import fromfile, sort, mgrid, float32, where, nan, copy, ones
from scipy.ndimage import filters
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.colors as colors


EUVDIR = 'C:\\Users\\smrak\\Google Drive\\BU\\software\\sdomask\\HiResFull300\\'

def getEUVMaskCoordinates(latlim=[-89.5,89.5],lonlim=[-180,180],nlat=180,nlon=360):
    xgrid, ygrid = mgrid[lonlim[0]:lonlim[1]:nlon*1j, latlim[0]:latlim[1]:nlat*1j]
    return xgrid,ygrid
def getEUVMaskFolder(f,nlat=180,nlon=360):
    xgrid, ygrid = getEUVMaskCoordinates(nlat=nlat, nlon=nlon)
    npts = nlat*nlon
    data = fromfile(f,count=npts, dtype=float32).reshape((nlat,nlon))
    return xgrid, ygrid, data
def getEUVMask(time,nlat=180,nlon=360):
    """
    I: time in posix
    """
    global EUVDIR
    xgrid, ygrid = getEUVMaskCoordinates(nlat=nlat, nlon=nlon)
    npts = nlat*nlon
    #Import EUV mask files
    flist = sort(glob(EUVDIR+'*.bin'))
    if isinstance(time,float) or isinstance(time,int):
        Tframe_full = datetime.utcfromtimestamp(time)
    elif isinstance(time, datetime):
        Tframe_full = time
    if int(Tframe_full.strftime('%H')) > 15 and int(Tframe_full.strftime('%H')) < 23:
        # find right filename extension
        TframeHM = Tframe_full.strftime('%H%M')
        flist = sort(glob(EUVDIR+'*'+TframeHM+'.bin'))
        # Get Mask
        data = fromfile(flist[0],count=npts, dtype=float32).reshape((nlat,nlon))
        return xgrid, ygrid, data
    else:
        return 0, 0, 0
        
def getMaskKeogram(nlat=180, nlon=360,X=39,lim=[-130, -60],direction='lat'):
    global EUVDIR
    npts = nlat*nlon
    xgrid, ygrid = getEUVMaskCoordinates(nlat=nlat, nlon=nlon)
    if direction == 'lat':
        idx = where( (xgrid[:,0] >= lim[0])  &  (xgrid[:,0] <= lim[1]))[0]
        idy = abs(ygrid[0,:] - X).argmin()
        Y = xgrid[idx,0]
    elif direction == 'lon':
        idx = abs(xgrid[:,0] - X).argmin()
        idy = where( (ygrid[0,:] >= lim[0])  &  (ygrid[0,:] <= lim[1]))[0]
        Y = ygrid[0,idy]
    
    flist = glob(EUVDIR+'*.bin')
    
    keogram = nan*ones((len(flist), idx.shape[0]))
    keograml = copy(keogram)
    dt = []
    i=0
    for f in flist:
        head, tail = path.split(f)
        hh = int(tail[6:8])
        mm = int(tail[8:10])
        if mm == 60:
            hh = int(hh)+1
            mm = 00
        dt.append(datetime.strptime('20170821'+str(hh)+str(mm), '%Y%m%d%H%M'))
        data = fromfile(f,count=npts, dtype=float32).reshape((nlat,nlon))
        laplace = filters.laplace(data)
        imslice = data[idy,idx] 
        lapslice = laplace[idy,idx]
        keogram[i,:] = imslice
        keograml[i,:] = lapslice
        i+=1

    return dt, Y, keogram, keograml
        
def makeBasemapMap(latlim=[20, 65], lonlim=[-160, -70], center=[39, -86],
            parallels=[20,30,40,50], lw=0.1,
            meridians = [-120,-110, -100, -90, -80,-70],
            fill_color='white', line_color='black',title='',
            epoto=False, totality=True, time=None,states=False,
            countries=False):
    """
    Plot the map and return handlers of the figure
    """
    (fig,ax) = plt.subplots(1,1,facecolor='w', figsize=(12,8))
    ax.set_title(title)
    m = Basemap(lat_0=40, lon_0=-95,llcrnrlat=latlim[0],urcrnrlat=latlim[1],
                llcrnrlon=lonlim[0],urcrnrlon=lonlim[1],
                projection='merc')#, resolution='i', ax=ax)
    m.drawmapboundary(fill_color=fill_color)
    
    m.drawcoastlines(color=line_color)
    if states: m.drawstates(color=line_color)
    if countries: m.drawcountries(color=line_color)
    m.drawparallels(parallels, labels=[1,0,0,0], linewidth=lw)
    m.drawmeridians(meridians, labels=[0,0,0,1], linewidth=lw)
    
    return fig, ax, m
    
def plotSDOMask(xgrid,ygrid,data,m,plot_type='contour',levels=[],
                colors='b',lw=1, cmap='jet',alpha=1, clabel=False,
                cbar=False,cbar_ticks=None,cbar_label=''):

    x,y = m(xgrid,ygrid)
    if plot_type == 'contour':
        
        if cmap is not None:
            CS = m.contour(x,y,data.T, levels, cmap=cmap, linewidths=lw,alpha=alpha)
        else:
            CS = m.contour(x,y,data.T, levels, colors=colors, linewidths=lw,alpha=alpha)
        if clabel:
                plt.clabel(CS, inline=1, fontsize=10)
    elif plot_type == 'contourf':
        m.contourf(x,y,data.T, levels, cmap=cmap)
    if cbar:
        if cbar_ticks is not None:
            cbar = m.colorbar(ticks=cbar_ticks)
        else:
            cbar = m.colorbar()
        if cbar_label != '':
            cbar.set_label(cbar_label)
    
    return m
def plotSatelliteTrack(m='',lat=[],lon=[],c='r',lw=1):
    """
    """
    x,y = m(lon,lat)
    m.plot(x,y,c=c,lw=lw)
    return m
    
def plotTotality(m,c='w',lw=0.5,alpha=1,
                     filepath='C:\\Users\\smrak\\Google Drive\\BU\\Projects\\Eclipse2017\\data\\totality.h5'):
    """
    Get the totality coordinates. Remark: this is a totality on the ground!
    Reference: NASA web page
    https://eclipse.gsfc.nasa.gov/SEpath/SEpath2001/SE2017Aug21Tpath.html
    """
    totality_path = File(filepath, 'r')
    north_lat = totality_path['path/north_lat'].value
    north_lon = totality_path['path/north_lon'].value
    south_lat = totality_path['path/south_lat'].value
    south_lon = totality_path['path/south_lon'].value

    X1,Y1 = m(north_lon, north_lat)
    X2,Y2 = m(south_lon, south_lat)
    m.plot(X1,Y1,c=c,lw=lw,alpha=alpha)
    m.plot(X2,Y2,c=c,lw=lw,alpha=alpha)
    return m

def plotLineMap(lon=[],lat=[],z=[], m='', lw=1, c='b'):
    
    x,y = m(lon-360,lat)
    m.plot(x,y, c=c, lw=lw)
    return m
def plotDMSPTemp(fig,ax,m='',lon=[],lat=[],z=[], scale=1, gamma=2,alpha=1,
                 cmap='Reds',cbar=False,cbar_label='', cbar_ticks=None):
    x,y = m(lon,lat)
    cs = m.scatter(x,y, s=z*scale, marker='_', c=z,cmap=cmap, 
                   norm=colors.PowerNorm(gamma=gamma), alpha=alpha)
    cs.set_clim(cbar_ticks[0], cbar_ticks[-1])
    if cbar:
        if cbar_ticks is None:
            cbar = m.colorbar(cs,location='bottom',pad="5%")
        else:
            cbar = m.colorbar(cs,location='bottom',pad="5%", ticks=cbar_ticks)
            cbar.set_clim(cbar_ticks[0], cbar_ticks[-1])
        if cbar_label != '':
            cbar.set_label(cbar_label)
    return m
    
def plotDMSPpos(lon=[],lat=[],m=''):
    x1,y1 = m(lon, lat)
    m.plot(x1,y1, 'h', ms=10, color='magenta')
    return m
#    levels = np.linspace(0.005,0.035,50)
#    im = abs(scipy.ndimage.filters.laplace(data))
#    m.contourf(x,y,im.T, levels, cmap='jet', norm=colors.PowerNorm(gamma=0.8))
#    cbar = m.colorbar(ticks=[0.005, 0.015, 0.025, 0.035])
#    cbar.set_label('EUV flux gradient')
#    
#    levels = np.linspace(0.2,1,30)
#    m.contour(x,y,data.T, levels, colors='r', alpha=0.7, linewidths=1)