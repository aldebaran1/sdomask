# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 13:11:52 2018

@author: ~Sebastijan Mrak "smrak@bu.edu"
"""

from numpy import arange, array, linspace, loadtxt
from scipy.ndimage import filters
from datetime import datetime, timezone, timedelta
import mask2d
import matplotlib.pyplot as plt

def getDMSPCoord(sv='F16'):
    if sv == 'F16':
        dmspfn = 'f16_rl172331909_extended_eclipse.txt'
    elif sv == 'F17':
        dmspfn = 'f17_rl172331842_extended_eclipse.txt'
    else:
        return
    #f17_rl172331842_extended_eclipse.txt
    #f16_rl172331909_extended_eclipse.txt
    #f16_rl172331909_eclipse_only.txt
    data = loadtxt(dmspfn, skiprows=3)# usecols=None, unpack=False, ndmin=0)[source]
    yydd = int(data[0][0])
    secinday = data[:,1]
    dmsp_dt = [datetime.strptime(str(yydd)+str(timedelta(seconds=t)), '%Y%j%H:%M:%S') for t in secinday]
    dmsp_lat = data[:,5]
    dmsp_lon = data[:,6]
    dmsp_Te = data[:,-2]
    
    return dmsp_dt, dmsp_lat, dmsp_lon, dmsp_Te


SAVEDIR = 'C:\\Users\\smrak\\Google Drive\\BU\\software\\sdomask\\maskdmsp17\\'
SAVEDIR = 'C:\\Users\\smrak\\Google Drive\\BU\\Projects\\Eclipse2017\\GRL2\\rev2\\'
sv = 'F17'
mins = arange(0,30)
mins = [0]
for minute in mins:
    dt = datetime(2017,8,21,19,minute,0)
    title_dt = dt.strftime('%m/%d/%Y -- %H:%M')
    figname = dt.strftime('%H%M')
    time = dt.replace(tzinfo=timezone.utc).timestamp()
    xgrid, ygrid, im = mask2d.getEUVMask(time)
    
    
    dmsp_dt, dmsp_lat, dmsp_lon, dmsp_Te = getDMSPCoord(sv=sv)
    
    if sv == 'F17':
        latlim = [-35, 60]
        lonlim = [-140,10.01]
        parallels = [-30,-10,10,30,50]
        meridians=[-140,-110,-80,-50,-20,10]
    elif sv =='F16':
        latlim = [-15, 60]
        lonlim = [-140,-19.99]
        parallels = [-10,10,30,50]
        meridians=[-140,-110,-80,-50,-20]
    # Map
    fig, ax, m = mask2d.makeBasemapMap(title='Time: '+title_dt+' UT',#+' DMSP: '+sv,
                                       latlim=latlim, lonlim=lonlim, parallels=parallels,
                                       meridians=meridians, fill_color='grey')
    # Plot totality path
    m = mask2d.plotTotality(m,c='r',lw=2,alpha=0.5)
    # Penumbra contours
    levels = linspace(0,1.1,40)
    m = mask2d.plotSDOMask(xgrid,ygrid,im,m, levels=levels, plot_type='contour',
                           colors='w',lw=0.4)
    # Laplacian
    im =abs(filters.laplace(im))
    levels = linspace(0.0035,0.03,10)
    im[im>=levels[-1]] = levels[-1]
    cbar_ticks_Te = [500, 1000, 1500, 2000, 2500]
    cbar_ticks_euv = [0.0035, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03]
    scale = 0.3
    # Plot Laplacian
    m = mask2d.plotSDOMask(xgrid,ygrid,im,m, levels=levels, plot_type='contourf',cmap='viridis',
                           cbar=True, cbar_label='grad(grad(EUV mask))', cbar_ticks=cbar_ticks_euv)
    # Plot Temp
#    m = mask2d.plotDMSPTemp(fig,ax,m=m,lon=dmsp_lon,lat=dmsp_lat,z=dmsp_Te,scale=scale,
#                            cmap='Reds',cbar=1,gamma=2,cbar_ticks=cbar_ticks_Te,cbar_label='Te @ DMSP '+sv)
#    
#    tdd = array([(tt-dt).total_seconds() for tt in dmsp_dt])
#    idt = abs(tdd).argmin()
#    m = mask2d.plotDMSPpos(lon=dmsp_lon[idt],lat=dmsp_lat[idt],m=m)
    
    plt.savefig(SAVEDIR+'map.png', dpi=600)
#    plt.close(fig)

