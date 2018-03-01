# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 14:22:13 2018

@author: smrak
"""

import mask2d
from numpy import linspace
from datetime import datetime, timezone
from scipy.ndimage import filters

t = datetime(2017,8,21,18,0,0)
time = t.replace(tzinfo=timezone.utc).timestamp()
xgrid, ygrid, data = mask2d.getEUVMask(time)

im =abs(filters.laplace(data))
lap_levels = linspace(0.0035,0.03,10)

latlim = [-15, 60]
lonlim = [-140,-19.99]
parallels = [-10,10,30,50]
meridians=[-140,-110,-80,-50,-20]


fig, ax, m = mask2d.makeBasemapMap(title='Obscuration mask',
                                       latlim=latlim, lonlim=lonlim, parallels=parallels,
                                       meridians=meridians, fill_color='grey')
m = mask2d.plotSDOMask(xgrid,ygrid,im,m, levels=lap_levels, plot_type='contourf',cmap='viridis',
                       cbar_label='grad(grad(EUV mask))')