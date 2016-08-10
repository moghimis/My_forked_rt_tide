#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
"""
script to process water elevation data and do tidal analysis with rt_tide

"""
from compiler.ast import Const

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
from okean import pl_tools
from okean import calc
import matplotlib.tri as mtri
import matplotlib.mlab as ml


#cmap = pl_tools.cm_ncview.bright
#cmap = pl_tools.cm_ncview.3saw
#cmap = pl_tools.cm_ncview.3gauss

data_dir = '/data01/01-projects/05-nasa-altimeter/01-progs/01-data/01-nri-pressure-gauges/02-bin/data/'
out_dir  = data_dir+'/out/'

pick_name = data_dir + '/tidal_data_const.pickle'
print 'Read pickle > ', pick_name
flow_data = pickle.load( open( pick_name , "r" ) )        




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


def get_depth(x,y,all=False):
    bat_in   = 'data/nri_regional_grd_roms_v2.nc'
    tnc    = netCDF4.Dataset(bat_in)
    tvars  = tnc.variables
    xb     = tnc.variables['x_rho'][:]
    yb     = tnc.variables['y_rho'][:]
    bat    = tnc.variables['h'][:]
    bat    = np.ma.masked_where(bat==-10.0,bat)
    if False:
        dep    = calc.griddata(xb[~bat.mask], yb[~bat.mask],bat[~bat.mask] , x, y,extrap=True)
    else:
        import octant.csa as csa
        csa_interp = csa.CSA(xb[~bat.mask], yb[~bat.mask],bat[~bat.mask])
        dep = csa_interp(x,y)
    
    if all:
        return dep,xb,yb,bat
    else:
        return dep




for sta_name in flow_data.keys():
    sta      = flow_data[sta_name]
    sta_tap  = sta['tappy_tide']
    sta_rtt  = sta['r_t_tide']
    sta_otp  = sta['otps_tide']
    
    sta_tap.rename(columns={'amp': 'amp_tap', 'pha': 'pha_tap'}, inplace=True)
    
    sta_tap['amp_rtt'] = sta_rtt['amp']
    sta_tap['pha_rtt'] = sta_rtt.pha
    
    sta_tap['amp_otp'] = sta_otp.amp
    sta_tap['pha_otp'] = sta_otp.pha
    
    sta_tap = sta_tap.set_index('Const')
    
    flow_data[sta_name]['tides'] = sta_tap

lon = []
lat = []
m2  = []
k1  = []
stas = []

for sta_name in flow_data.keys():
    lon.append(flow_data[sta_name]['lon'])
    lat.append(flow_data[sta_name]['lat'])
    m2.append (flow_data[sta_name]['tides'].loc['M2'].values)
    k1.append (flow_data[sta_name]['tides'].loc['K1'].values)
    stas.append(sta_name)    
names =  flow_data[sta_name]['tides'].columns.values   
lon   = np.array(lon)
lat   = np.array(lat)
m2    = np.array(m2)
k1    = np.array(k1)
stas  = np.array(stas)

proj = projection_nri ()
x,y    = lonlat2xy(proj,lon,lat)

#otps tides        
tide_file = data_dir +  '/tides_inner_v1_2012-05-01_v2.nc'
ncf       = netCDF4.Dataset(tide_file,'r')
ncvar     = ncf.variables

const1    = ncvar['tidal_constituents'][:]
xa      = ncvar['x_rho'][:]
ya      = ncvar['y_rho'][:]

nam     = [] 
for i in range(len(const1)):
    nam.append( string.capitalize (const1[i][0]+const1[i][1] ) )

nam = np.array(nam) 


#read bathy local
bat_local = data_dir +  '/local_v1.nc'
ncb       = netCDF4.Dataset(bat_local,'r')
ncvb      = ncb.variables
xl        = ncvb['x_rho'][:]
yl        = ncvb['y_rho'][:]


jmin = 50  
jmax = 220

lim1 = 0.5
lim2 = 0.54
levels = np.linspace(lim1, lim2, 5)
fmt = '%.4g'

#tri = mtri.Triangulation(x, y)



for const in  ['M2']:
    [ic]    = np.where(const == nam)
    amp_all = ncvar['tide_Eamp']  [ic].squeeze() 
    pha_all = ncvar['tide_Ephase'][ic].squeeze() 
    
    fig,axgrid = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False,
            figsize=(7, 12),facecolor='w', edgecolor='k',dpi=450 )
    axgrid = np.array(axgrid)
    
    iax = 0
    
    titles = ['Tappy Amp','Tappy Pha']
    otps   = [amp_all,pha_all]
    lims   = [(0.5,0.5),(220,226)]
    
    
    for data1 in [m2[:,0],m2[:,1]]:
        print titles[iax]
        ax1 = axgrid[iax]
        ax1.set_title(titles[iax])
        #p1 = plt.contourf(xa[jmin:jmax,:],ya[jmin:jmax,:],amp_all[jmin:jmax,:],
        #             levels = levels,cmap = pl_tools.cm_ncview.banded,extend='both')
        
        p1 = ax1.pcolor(xa[jmin:jmax,:],ya[jmin:jmax,:],otps[iax][jmin:jmax,:],
                     cmap = pl_tools.cm_ncview.banded)    
        
        
        
        #plt.colorbar(p1)
        p1.set_clim(lims[iax])
        
        # set scale for vertical vector plots
        pos_ax   = np.array (ax1.get_position ())
        #
        xsub1=pos_ax[0][0]
        xsub2=pos_ax[1][0]
        ysub1=pos_ax[0][1]            
        ysub2=pos_ax[1][1]
        cbax    = fig.add_axes([xsub1+ 1.01 * (xsub2-xsub1), ysub1 + 0.14 * (ysub2-ysub1), 0.03, 0.2]) 
        cb      = plt.colorbar(p1,cax=cbax,format='%1.4g',orientation='vertical') #
        #cb.set_label(cbtext)
        

        
        
        p2 = ax1.scatter(x, y, s = 5, c = data1,edgecolors='None',cmap = pl_tools.cm_ncview.banded )
        #p2.set_clim(lims[iax])
        #ax1.colorbar(p2)

        for i in range(len(x)):
            ax1.text   (x[i], y[i], stas[i],fontsize = 6 )
        
        
        
        #plot contour line from obs
        #interp_lin = mtri.LinearTriInterpolator(tri,m2[:,0] )   #m2 amp tappy
        #m2_local   = interp_lin(xl, yl)
        #m2_local = ml.griddata(x, y,m2[:,0] , xl, yl)
        #m2_local = calc.griddata(x, y,data1 , xl, yl,extrap=True)
        #cl1 = ax1.contour(xl,yl,m2_local,colors='black')
        #ax1.clabel(cl1,fontsize=10,inline=True,fmt=fmt)
        iax += 1
        

plt.savefig('test.png',dpi=450)
# os.system('mkdir -p '+outdirf+'/0scr')
# #Back_up scr
# args     = sys.argv
# scr_name = args[0]    
# scr_dir1 = os.getcwd()
# os.system('cp -fr  '+scr_name +'    '+outdirf+'/0scr')
# #os.system('cp -fr  '+inp_dir+'/*' +'    '+outdirf+'/0scr')
plt.close('all')

dep,xb,yb,bat = get_depth(x,y,all=True)

xarg = np.argsort(x)

fig2 = plt.figure()
ax = fig2.add_subplot(111)
ax.set_title('M2 Amp. [m]')
#ax.plot(m2[:,0][xarg],'r',lw=1,label ='Tappy' )
ax.plot(m2[:,2][xarg],'b',lw=1,label ='RTT' )
ax.plot(m2[:,4][xarg],'g',lw=1,label ='OTPS' )
ax.xaxis.set_ticklabels(stas[xarg])
#ax.set_ylim(0.5,0.6)
ax.legend(frameon=False)


ax2 = ax.twinx()
ax2.plot(dep[xarg],'k')
ax2.set_ylabel('Depth [m]')
#ax2.set_ylim(5,8)

plt.savefig('m2_amp.png',dpi=450)


fig3 = plt.figure()
ax = fig3.add_subplot(111)
ax.set_title('M2 Pha. [deg]')
#ax.plot(m2[:,1][xarg],'r',lw=1,label ='Tappy' )
ax.plot(m2[:,3][xarg],'b',lw=1,label ='RTT' )
ax.plot(m2[:,5][xarg],'g',lw=1,label ='OTPS' )
ax.xaxis.set_ticklabels(stas[xarg])
#ax.set_ylim(220,240)

ax.legend()

ax2 = ax.twinx()
ax2.plot(dep[xarg],'k')
ax2.set_ylabel('Depth [m]')
#ax2.set_ylim(5,8)
plt.savefig('m2_pha.png',dpi=450)



fig2 = plt.figure()
ax = fig2.add_subplot(111)
ax.set_title('K1 Amp. [m]')
#ax.plot(k1[:,0][xarg],'r',lw=1,label ='Tappy' )
ax.plot(k1[:,2][xarg],'b',lw=1,label ='RTT' )
ax.plot(k1[:,4][xarg],'g',lw=1,label ='OTPS' )
ax.xaxis.set_ticklabels(stas[xarg])
#ax.set_ylim(0.5,0.6)

ax.legend()
ax2 = ax.twinx()
ax2.plot(dep[xarg],'k')
ax2.set_ylabel('Depth [m]')
#ax2.set_ylim(5,8)
plt.savefig('k1_amp.png',dpi=450)



fig3 = plt.figure()
ax = fig3.add_subplot(111)
ax.set_title('K1 Pha. [deg]')
#ax.plot(k1[:,1][xarg],'r',lw=1,label ='Tappy' )
ax.plot(k1[:,3][xarg],'b',lw=1,label ='RTT' )
ax.plot(k1[:,5][xarg],'g',lw=1,label ='OTPS' )
ax.xaxis.set_ticklabels(stas[xarg])
#ax.set_ylim(220,240)

ax.legend()
ax2 = ax.twinx()
ax2.plot(dep[xarg],'k')
ax2.set_ylabel('Depth [m]')
#ax2.set_ylim(5,8)

plt.savefig('k1_pha.png',dpi=450)




fig3 = plt.figure()
ax = fig3.add_subplot(111)
p1 = ax.contourf(xb,yb,bat,levels=np.linspace(0,8,11),extend='both')
plt.colorbar(p1)
p2 = ax.scatter(x, y, s = 5, c = 'k')
for i in range(len(x)):
    ax.text   (x[i], y[i], stas[i],fontsize = 12 )

ax.set_xlim(-2000,2000)
ax.set_ylim(-2000,2000)
ax.set_aspect(1)

plt.savefig('sta_map.png',dpi=450)
plt.show()
