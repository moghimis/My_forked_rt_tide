#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Main waterlevel2nc script ####
"""
__author__ = "Saeed Moghimi"
__copyright__ = "Copyright 2013, Oregon State University"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "moghimis@gmail.com"

import string
import netCDF4
import netcdftime
import numpy as np
import pylab as pl
import datetime as datetime
import sys,os
from dateutil.parser import parse

from scipy.signal import firwin,lfilter

### EDT – Eastern Daylight Time = Time zone offset: UTC - 4 hours
### here dates are in GMT


inpdir= '../inp/wl_8658163_Wrightsville_Beach_NC/'
#inpdir= '../inp/wl_8656483_Beaufort_NC/'



date_first = datetime.datetime(2012,5,1)
date_part=date_first.isoformat()[:10]
hour_part=date_first.isoformat()[11:]
units_start="seconds since "+date_part+" "+hour_part
utime=netcdftime.utime(units_start)

####################################################
finfo  = open(inpdir+'info')
readme = finfo.readlines()
finfo.close()

finfo = open(inpdir+'info')
for line in finfo:
    words = string.split(line)
    if 'lat ' in line: lat = float(words[1])
    if 'lon ' in line: lon = float(words[1])
finfo.close()
    



#filename2 = 'inp/bouy-cdip190_noaa41109may2012.txt'
filename2 = inpdir + 'data.csv'
print 'Read observation at: ', filename2

fp = open(filename2, "r")
line=''
line = fp.readline()

obs = []
obs_dates = []
for line in fp:
    words = string.split(line,',')
    #Date Time, Water Level, Sigma, I, L 
    try:
       obs_dates.append(parse(words[0]))
    except:
       break
    wl =float(words[1])            # m
    obs.append(wl)

fp.close()

obs=pl.array(obs)
obs_dates=pl.array(obs_dates)
obs_time=utime.date2num(obs_dates)

outfile=filename2[:-3]+'nc'
outnc=netCDF4.Dataset(outfile,'w',format='NETCDF3_CLASSIC')
outnc.history  = 'SM measur2netCDF 2012'
outnc.latitude =  lat
outnc.longitude = lon
outnc.depth = obs.mean()
outnc.createDimension('time',None)

timea = outnc.createVariable('time','f8',('time',))
timea.units = units_start
timea[:]= obs_time
p0 = outnc.createVariable('elev','f8',('time',))
p0.units = 'm'
p0.missing_value = -99.0
p0[:]=obs

p2 = outnc.createVariable('lat','f8')
p2.units = 'deg'
p2.missing_value = -99.0
p2[:]=lat

p3 = outnc.createVariable('lon','f8')
p3.units = 'deg'
p3.missing_value = -99.0
p3[:]=lon


outnc.close()
print 'end > ' , filename2[:-3]+'nc', lat,lon

