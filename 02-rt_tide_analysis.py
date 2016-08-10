#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
"""
script to process water elevation data and do tidal analysis with rt_tide

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

import matlab.engine
eng = matlab.engine.start_matlab()
eng.addpath(r'/data01/01-projects/01-passaic/03-modeling/01-delft3d/02-pass-ideal/pycodes/post/06-r_t_tide/mcodes/r_t_tide',nargout=0);

data_dir = '/data01/01-projects/05-nasa-altimeter/01-progs/01-data/01-nri-pressure-gauges/02-bin/data/'
out_dir  = data_dir+'/out/'

pick_name = data_dir + '/tidal_data.pickle'
print 'Read pickle > ', pick_name
flow_data = pickle.load( open( pick_name , "r" ) )        


def do_r_t_tide_analysis(flow_data):
    """
    only works for hourly data
    use resample function in dataframe to resample it to 'H'
    
    TODO: remove hourly constraint
    """
       
    for sta_name in flow_data.keys():
        print sta_name
        r_t_tide_out_dir = out_dir + '/r_t_tide_out_txt_files/' 
        os.system('mkdir -p '+ r_t_tide_out_dir)
        
        #prepare pickle name !
        pick_dir = r_t_tide_out_dir + '/r_t_tide_pickle/' 
        os.system('mkdir -p '+ pick_dir)
        pick_name = pick_dir + 'rt_tide_data.pickle'
        
        t_sta  = np.array(flow_data[sta_name]['datenum_mat'])
        h_sta  = np.array(flow_data[sta_name]['elev'])
        [rk]   = np.where( ~np.isnan(h_sta) )
        t_tmp  = t_sta [rk].tolist()
        h_tmp  = h_sta [rk].tolist()
        out_file = r_t_tide_out_dir +flow_data[sta_name]['name_long'].replace(' ','_')+'.rtt'
        
        #run matlab r_t_tide
        
        try:
            lat =  flow_data[sta_name]['lat'].item()
        except:
            lat = flow_data[sta_name]['lat']
            
        
        s = eng.r_t_tide_saeed    (matlab.double(t_tmp),matlab.double(h_tmp),'latitude',
            lat ,'nodalcorrflag','true','greenwichcorrflag','true','method','cauchy',
            'output',out_file,'rayleigh',0.9 );
        
        
        #vectors to pick the requested constituents
        sta_cnam   = []
        sta_amp    = []
        sta_amp_ci = []
        sta_pha    = []
        sta_pha_ci = []
        
        
        #read t_tide_text file
        f_rtt = open(out_file)
        nam    = []
        amp    = []
        amp_ci = []
        pha    = []
        pha_ci = []
        start_read = False
        for line in f_rtt.readlines():
            if 'pha_err' in line: 
                start_read = True
                continue
            if start_read:
                tmp = line.split()
                nam.append(line[1:5].rstrip())
                amp.append(   float(tmp[2]))
                amp_ci.append(float(tmp[3]))
                pha.append(   float(tmp[4]))
                pha_ci.append(float(tmp[5]))
        f_rtt.close()
        
        nam    = np.array(nam)
        amp    = np.array(amp)
        amp_ci = np.array(amp_ci)
        pha    = np.array(pha)
        pha_ci = np.array(pha_ci)
        
        ##    
        for const in  constits:
            sta_cnam.append   (const)
            if const in nam:
                [ic] = np.where(const == nam)
                sta_amp.append   (amp   [ic])
                sta_amp_ci.append(amp_ci[ic])
                sta_pha.append   (pha   [ic])
                sta_pha_ci.append(pha_ci[ic])
            else:
                cmask = np.array([9.999e-12])
                sta_amp.append   (cmask)
                sta_amp_ci.append(cmask)
                sta_pha.append   (cmask)
                sta_pha_ci.append(cmask)                    #if len(const) > 2:
        
        sta_cnam    = np.array(sta_cnam).squeeze()  
        sta_amp     = np.array(sta_amp).squeeze() 
        sta_amp_ci  = np.array(sta_amp_ci).squeeze() 
        sta_pha     = np.array(sta_pha).squeeze() 
        sta_pha_ci  = np.array(sta_pha_ci).squeeze() 
        sta_all = np.array(zip(sta_cnam,sta_amp,sta_amp_ci,sta_pha,sta_pha_ci))                    
        sta_df  = pd.DataFrame(data=sta_all, columns=['Const','amp','amp_err','pha','pha_err'])
        
        sta_df['amp'        ] = sta_df['amp'    ].convert_objects(convert_numeric=True)
        sta_df['amp_err'    ] = sta_df['amp_err'].convert_objects(convert_numeric=True)
        sta_df['pha'        ] = sta_df['pha'    ].convert_objects(convert_numeric=True)
        sta_df['pha_err'    ] = sta_df['pha_err'].convert_objects(convert_numeric=True)
        
        flow_data[sta_name]['r_t_tide'] = sta_df
        
    
def do_tappy_tide_analysis(flow_data):
    from tappy_local import tappy
    ### Saeed tries to understand!  from here
    for sta_name in flow_data.keys():
        dates      = np.array(flow_data[sta_name]['date'])
        elev       = np.array(flow_data[sta_name]['elev'])
        lon        = np.array(flow_data[sta_name]['lon'])
        lat        = np.array(flow_data[sta_name]['lat'])
        
        tappy_tide_out_dir = out_dir + '/tappy_tide_out_txt_files/' 
        os.system('mkdir -p '+ tappy_tide_out_dir)
        out_file = tappy_tide_out_dir + flow_data[sta_name]['name_long'].replace(' ','_')+'.tappy'
        
        
        ### Saeed tries to understand!  from here
        data_filename    ='test'
        def_filename     = None
        config           = None
        quiet            = False
        debug            = False
        outputts         = False
        outputxml        = ''
        ephemeris        = False
        rayleigh         = 0.9
        print_vau_table  = False
        missing_data     = 'ignore'
        #missing_data    = 'fill'
        linear_trend     = False
        remove_extreme   = False
        zero_ts          = None
        filter           = None
        pad_filters      = None
        include_inferred = True
        xmlname          = flow_data[sta_name]['name_long']
        xmlcountry       = 'US'
        xmllatitude      = lat
        xmllongitude     = lon
        xmltimezone      = '0000'
        xmlcomments      = 'No comment'
        xmlunits         = 'm or ms-1'
        xmldecimalplaces = None
        
        ############## model
        x = tappy(
            outputts  = outputts,
            outputxml = 'model.xml',
            quiet     = quiet,
            debug     = debug,
            ephemeris = ephemeris,
            rayleigh  = rayleigh,
            print_vau_table = print_vau_table,
            missing_data = missing_data,
            linear_trend = linear_trend,
            remove_extreme = remove_extreme,
            zero_ts = zero_ts,
            filter  = filter,
            pad_filters = pad_filters,
            include_inferred = include_inferred,
            )
        
        x.dates      = dates
        x.elevation  = elev
        package      = x.astronomic(x.dates)
        (x.zeta, x.nu, x.nup, x.nupp, x.kap_p, x.ii, x.R, x.Q, x.T, x.jd, x.s, x.h, x.N, x.p, x.p1) = package
        ray = 1.0
        (x.speed_dict, x.key_list) = x.which_constituents(len(x.dates),package,rayleigh_comp = ray)
        
        x.constituents()
        x.print_con()
        x.print_con_file(filedat = out_file, lon = lon, lat = lat)
        only_const = False
        
        if not only_const:
            if x.missing_data == 'fill':
                    x.dates_filled, x.elevation_filled = x.missing(x.missing_data, x.dates, x.elevation)
                    x.write_file( x.dates_filled,
                                    x.elevation_filled,
                                    fname='outts_filled.dat')
            
            x.filter = 'usgs'
            #x.filter='doodson'
            
            if x.filter:
                    for item in x.filter.split(','):
                        if item in ['mstha', 'wavelet', 'cd', 'boxcar', 'usgs', 'doodson', 'lecolazet1', 'kalman', 'transform']:# 'lecolazet', 'godin', 'sfa']:
                            filtered_dates, result = x.filters(item, x.dates, x.elevation)
                            x.write_file(filtered_dates, result, fname='outts_filtered_%s.dat' % (item,))
                            x_dates_filter= filtered_dates
                            x_eleva_filter= result              
                    (x.speed_dict, x.key_list) = x.which_constituents(len(x.dates),package,rayleigh_comp = ray)
            
        #vectors to pick the requested constituents
        sta_cnam   = []
        sta_amp    = []
        sta_pha    = []
        
        #read t_tide_text file
        f_rtt = open(out_file)
        nam    = []
        amp    = []
        pha    = []
        start_read = False
        for line in f_rtt.readlines():
            if '=        =====' in line: 
                start_read = True
                continue
            
            if len(line) < 3:
                break
        
            if start_read:
                tmp = line.split()
                nam.append(line[9:14].strip())
                amp.append(   float(tmp[2]))
                pha.append(   float(tmp[3]))
        
        f_rtt.close()
        
        nam    = np.array(nam)
        amp    = np.array(amp)
        pha    = np.array(pha)
        
        ##    
        for const in  constits:
            sta_cnam.append   (const)
            if const in nam:
                [ic] = np.where(const == nam)
                sta_amp.append   (amp   [ic].item())
                sta_pha.append   (pha   [ic].item())
            else:
                cmask = np.array([9.999e-12])
                sta_amp.append   (cmask)
                sta_pha.append   (cmask)
        
        sta_cnam = np.array(sta_cnam).flatten()
        sta_amp  = np.array(sta_amp ).flatten()
        sta_pha  = np.array(sta_pha ).flatten()  
        
        sta_all  = np.array(zip(sta_cnam,sta_amp,sta_pha))                    
        sta_df   = pd.DataFrame(data = sta_all, columns=['Const','amp','pha'])
        
        sta_df['amp'] = sta_df['amp'].convert_objects(convert_numeric=True)
        sta_df['pha'] = sta_df['pha'].convert_objects(convert_numeric=True)
        
        flow_data[sta_name]['tappy_tide'] = sta_df
        print ' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> '


#####   ___ main ___
#NOW, Do piecwise HA of tide data from 2009 to 2010

#constits = ['M2','N2','S2','K2' ,'P1' ,'Q1' ,'K1' ,'O1',
#           'M4','MS4','MN4','MK4','MK3','SO3','MO3','S4',
#           'M8','J1' ,'T2' ,'L2' ,'2N2','MU2','NU2','M3',
#           'M6','SSA','SA']

constits = ['M2','N2','S2','Q1','K1','K2','P1','O1','M4']
do_r_t_tide_analysis   (flow_data = flow_data)
do_tappy_tide_analysis (flow_data = flow_data)
#
#

tide_file = data_dir +  '/tides_inner_v1_2012-05-01_v2.nc'
ncf       = netCDF4.Dataset(tide_file,'r')
ncvar     = ncf.variables

const1    = ncvar['tidal_constituents'][:]
lona      = ncvar['lon_rho'][:]
lata      = ncvar['lat_rho'][:]

nam     = [] 
for i in range(len(const1)):
    nam.append( string.capitalize (const1[i][0]+const1[i][1] ) )

nam = np.array(nam) 
   
for sta_name in flow_data.keys():
    lon   = flow_data[sta_name]['lon']
    lat   = flow_data[sta_name]['lat']
    dist  = ((lona - lon)**2 + (lata - lat)**2)
    [i,j] = np.where(dist==dist.min())
    
    amp   =  ncvar ['tide_Eamp'  ][:,i,j].squeeze()
    pha   =  ncvar ['tide_Ephase'][:,i,j].squeeze()
    
    #vectors to pick the requested constituents
    sta_cnam   = []
    sta_amp    = []
    sta_pha    = []
    
    ##    
    for const in  constits:
        sta_cnam.append(const)
        if const in nam:
            [ic] = np.where(const == nam)
            sta_amp.append (amp[ic].item())
            sta_pha.append (pha[ic].item())
        else:
            cmask = np.array([9.999e-12])
            sta_amp.append (cmask)
            sta_pha.append (cmask)
    
    sta_cnam = np.array(sta_cnam).flatten()
    sta_amp  = np.array(sta_amp ).flatten()
    sta_pha  = np.array(sta_pha ).flatten()  
    
    sta_all  = np.array(zip(sta_cnam,sta_amp,sta_pha))                    
    sta_df   = pd.DataFrame(data = sta_all, columns=['Const','amp','pha'])
    
    sta_df['amp'] = sta_df['amp'].convert_objects(convert_numeric=True)
    sta_df['pha'] = sta_df['pha'].convert_objects(convert_numeric=True)
    
    flow_data[sta_name]['otps_tide'] = sta_df

#sta      = flow_data['p88']
#sta_tap  = sta['tappy_tide']
#sta_rtd  = sta['r_t_tide']


pick_name = data_dir + '/tidal_data_const.pickle'
print ' > Write pickle > ', pick_name
pickle.dump( flow_data, open(pick_name , "wb" ) )




# os.system('mkdir -p '+outdirf+'/0scr')
# #Back_up scr
# args     = sys.argv
# scr_name = args[0]    
# scr_dir1 = os.getcwd()
# os.system('cp -fr  '+scr_name +'    '+outdirf+'/0scr')
# #os.system('cp -fr  '+inp_dir+'/*' +'    '+outdirf+'/0scr')


print ' End >>>'



