
import sys, string
import os as os
import pylab as pl
import datetime
import netcdftime
import netCDF4

import struct
import glob


#read stations specifications
stations={}
filename1='inp/sensor_locations.inp'
fpt = open(filename1, "r")
fpt.readline()
fpt.readline()
fpt.readline()

for i in range(32):
    line=fpt.readline()
    words = string.split(line)
    #print  i,words
    id=words[0]
    lat   =float(words[1])   
    lon   =float(words[2])   
    qdist =float(words[3])  # in cm
    qdist=qdist/100.         # cm --> m
    
    data=pl.array([lat,lon,qdist])
    item={id:data}
    stations.update(item)

dir1=os.getcwd()
dir_cdf =dir1+'/../../../data/rev02_oct2012/RIVET_Britt_ADV_cdf/'

dircdf_cat_final=dir_cdf+'final_cat/'
comm0='mkdir -p '+dircdf_cat_final
os.system(comm0) 


#tmp=dir_cdf+'cat*'
#dirlist=glob.glob(tmp)
#dirlist.sort()

id_all=stations.keys()
id_all.sort()
for id1 in id_all:
    for param1 in ['u','v','w','p','q']:
        print 'cat  '+param1+id1
        comm1='ncrcat -O  '+dir_cdf+'cat*/*'+param1+id1+'.nc  '+dircdf_cat_final+param1+id1+'_0428-0531.nc'
        os.system(comm1)

