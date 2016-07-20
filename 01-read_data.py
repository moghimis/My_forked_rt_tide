#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
"""
script to read water elevation data and make them ready for tidal analysis

"""

__author__ = "Saeed Moghimi"
__copyright__ = "Copyright 2016, Portland State University"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "moghimis@gmail.com"


import netCDF4
import netcdftime
from datetime import datetime,timedelta
from dateutil.parser import *
import scipy.io as sio
from   collections  import defaultdict
import numpy as np
import os,sys
import pandas as pd
import cPickle as pickle
import string
import matplotlib.pyplot as plt
import  pynmd.plotting.plot_settings as ps
from  pynmd.plotting.plot_settings import mat2py_datenum
import glob

data_dir = '/data01/01-projects/05-nasa-altimeter/01-progs/01-data/01-nri-pressure-gauges/02-bin/data/'
#
lon_orig = -77.3382
lat_orig =  34.5279
#
def projection_nri ():
    from mpl_toolkits.basemap import Basemap
    proj    = Basemap(resolution=None,projection='npstere',lon_0=-45,boundinglat=70)
    return proj
#
def lonlat2xy(proj,lon,lat):
    xu, yu = proj(lon, lat)
    xu_orig,yu_orig=proj(lon_orig, lat_orig)
    x = xu - xu_orig
    y = yu - yu_orig
    return x,y
#
def read_Elgar_pressure_data(flow_data):
    print '   > Elgar pressure data;'
    dir1 = data_dir + '/elgar/*_filter_abs.nc'
    flist    = glob.glob(dir1)
    flist.sort()

    for filename in flist[:]:
        ncf   = netCDF4.Dataset(filename,'r')
        ncvar = ncf.variables
        time       = ncvar['time']
        utime      = netcdftime.utime(time.units)
        dates      = utime.num2date(time[:])
        sta_name   = filename.split('/')[-1][:3]
        if 'p' in sta_name or 'q' in sta_name:
            print '    > read  > Station name: ', sta_name
    
            water_depth = ncvar['water_depth'][:]
    
            data  = pd.DataFrame(data = water_depth, columns = ['water_depth'], index = dates)    
            data  = data.dropna()
            data  = data.resample('H')  #hourly mean
            
            date_tmp =  ps.datetime64todatetime(data.index.values)
            datenum  = [ps.datetime2datenum(datei) for datei in  date_tmp]

            
            flow_data[sta_name]['name_long'] = 'Bottom pressure '+sta_name
            flow_data[sta_name]['lat']       = ncvar['lat'][:]
            flow_data[sta_name]['lon']       = ncvar['lon'][:]
            flow_data[sta_name]['elev']      = (data['water_depth'] - data['water_depth'].mean()).values
            flow_data[sta_name]['date']      = date_tmp
            flow_data[sta_name]['datenum_mat']   = datenum

#
def read_wl_noaa_station(flow_data,inp_dir,sta_name):
    print  inp_dir
    dir1  = data_dir + inp_dir + '/'
    finfo = open(dir1 + 'info')
    for line in finfo:
        words = string.split(line)
        if 'lat ' in line: lat = float(words[1])
        if 'lon ' in line: lon = float(words[1])
    finfo.close()
    
    filename2 = dir1 + 'data.csv'
    print 'Read observation at: ', filename2
    
    fp   = open(filename2, "r")
    line = ''
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
        wl = float(words[1])            # m
        obs.append(wl)
    
    fp.close()
    
    obs       = np.array(obs)
    obs_dates = np.array(obs_dates)
    
    print '    > read  > Station name: ', sta_name
    
    water_depth = obs
    data  = pd.DataFrame(data = water_depth, columns = ['water_depth'], index = obs_dates)    
    data  = data.dropna()
    data  = data.resample('H')  #hourly mean
    
    date_tmp =  ps.datetime64todatetime(data.index.values)
    datenum  = [ps.datetime2datenum(datei) for datei in  date_tmp]
    
    flow_data[sta_name]['name_long'] = inp_dir
    flow_data[sta_name]['lat']       = lat
    flow_data[sta_name]['lon']       = lon
    flow_data[sta_name]['elev']      = (data['water_depth'] - data['water_depth'].mean()).values
    flow_data[sta_name]['date']      = date_tmp
    flow_data[sta_name]['datenum_mat']   = datenum

#
def read_ncom(flow_data):
    print '  > NCOM;'
    filename  = data_dir + '/ncom/ncom_glb_regp01_2012-3-4-5.nc'

    ncf   = netCDF4.Dataset(filename,'r')
    ncvar = ncf.variables
    time       = ncvar['time']
    utime      = netcdftime.utime(time.units)
    dates      = utime.num2date(time[:])
    sta_name   = filename.split('/')[-1][:3]
    
    lona,lata = np.meshgrid(ncvar['lon'][:],ncvar['lat'][:])
    dist = ((lona - lon_orig)**2 + (lata - lat_orig)**2)
    [i,j] = np.where(dist==dist.min())
    
    elev  = ncvar['surf_el'][:,i,j].flatten()
    data  = pd.DataFrame(data = elev, columns = ['elev'], index = dates)    
    data  = data.dropna()
    data  = data.resample('H')  #hourly mean
    
    date_tmp =  ps.datetime64todatetime(data.index.values)
    datenum  = [ps.datetime2datenum(datei) for datei in  date_tmp]
    
    flow_data[sta_name]['name_long'] = 'Elevation from NCOM model '
    flow_data[sta_name]['lat']       = lata[i,j].item()
    flow_data[sta_name]['lon']       = lona[i,j].item()
    flow_data[sta_name]['elev']      = data['elev'].values
    flow_data[sta_name]['date']      = date_tmp
    flow_data[sta_name]['datenum_mat']   = datenum
#
#main code
#
flow_data = defaultdict(dict)
read_Elgar_pressure_data(flow_data = flow_data)
#read_wl_noaa_station( flow_data = flow_data, inp_dir = 'wl_8658163_Wrightsville_Beach_NC', sta_name = 'wri')
#read_wl_noaa_station( flow_data = flow_data, inp_dir = 'wl_8656483_Beaufort_NC'          , sta_name = 'bea')
#read_ncom(flow_data)

proj = projection_nri ()
#add xy
for sta_name in flow_data.keys():
    flow_data[sta_name]['x'],flow_data[sta_name]['y'] = lonlat2xy(proj,flow_data[sta_name]['lon'],flow_data[sta_name]['lat'])


pick_name = data_dir + '/tidal_data.pickle'
print ' > Write pickle > ', pick_name
pickle.dump( flow_data, open(pick_name , "wb" ) )




