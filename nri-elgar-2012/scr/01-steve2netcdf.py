
import sys, string
import os as os
import pylab as pl
import datetime
import netcdftime
import netCDF4

import struct
import glob

dir1=os.getcwd()

##### reading bathymetry file
fngrd='inp/topo10m_nri.nc'  #f(time, yc, xc)
fgrd = netCDF4.Dataset(fngrd)    
nvargrd=fgrd.variables
imin=600
imax=3000
jmin=500
jmax=2000
latt=nvargrd['latc'][imin:imax]
lont=nvargrd['lonc'][jmin:jmax]
lat2,lon2=pl.meshgrid(latt, lont)
bat=nvargrd['bathymetry'][jmin:jmax,imin:imax]
fgrd.close()
##################################################3

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
    
    ### to find depth ####
    dist2=pl.sqrt ((lat2-lat)**2+(lon2-lon)**2)
    distm2=dist2.min()
    idx2=pl.where(dist2==distm2)    
    jnum,inum=idx2
    dep=bat[jnum,inum]
    print dep ,lat ,lat2[jnum,inum],  lon ,lon2[jnum,inum]
    ######################
    
    data=pl.array([lat,lon,qdist,dep])
    item={id:data}
    stations.update(item)

nt=6144
dir_data=dir1+'/../../../data/rev02_oct2012/RIVET_Britt_ADV/*'
dir_cdf =dir1+'/../../../data/rev02_oct2012/RIVET_Britt_ADV_cdf/'

dirlist=glob.glob(dir_data)
dirlist.sort()
for dir in dirlist[:]:
    dirin=dir+'/*'
    dircdf=dir_cdf+dirin[-6:-1]
    command='mkdir -p ' + dircdf
    os.system(command) 
    flist=glob.glob(dirin)
    flist.sort()
    for file in flist[:]:
        print '   ',file[-15:]
        
        infile=file
        if os.path.getsize(file)/1024 == 16: nt= 4096
        fid = open(infile, 'rb')
        data=struct.unpack('f'*nt, fid.read(4*nt))
        fid.close()
        data=pl.array(data)     # in cm/s or in cm
        data=data/100.          # cm --> m
        ### Constructing netCDF file stuff
        #start time:
        units ='seconds since 2012-04-01 00:00:00'
        base_date=datetime.datetime(2012,04,01) 
        utime=netcdftime.utime(units)
      
        year1=2012
        mont1=int(infile[-12:-10])
        day1 =int(infile[-10:-8]) 
        hour1=int(infile[-8:-6])
        min1 =int(infile[-6:-4])
        obs_start_time = datetime.datetime(year1,mont1,day1,hour1,min1) 
        time_dif=obs_start_time-base_date
        time_dif_sec=time_dif.total_seconds()
        #units ='seconds since '+year1+'-'+mont1+'-'+day1+' '+hour1+':'+min1+':'+'00'
        
        #file name
        cdfname=dircdf+'data_'+infile[-12:]+'.nc'
        outnc=netCDF4.Dataset(cdfname,'w',format='NETCDF3_CLASSIC')
        outnc.createDimension('time',None)
        timea = outnc.createVariable('time','f8',('time',))
        timea.units = units
        time1=pl.linspace(0,3072,nt)
        time1=time1+time_dif_sec
        timea[:]=time1
        
        param=infile[-3:-2]
        id_current=infile[-2:] 
        lat_current=stations[id_current][0]
        lon_current=stations[id_current][1]
        qdist_current=stations[id_current][2]
        depth_current=stations[id_current][3]

        pdist_current=44.5/100.
        veldist_current= 78./100.
        if   param=='q':
            param_name='water_depth'
            data=data + qdist_current
            dist_above_bottom=qdist_current
        elif param=='p':
            param_name='water_depth'
            data=data + pdist_current
            dist_above_bottom=pdist_current
        elif param=='u':
            param_name='u' 
            dist_above_bottom=veldist_current
        elif param=='v':
            param_name='v'
            dist_above_bottom=veldist_current
        elif param=='w':
            param_name='w'
            dist_above_bottom=veldist_current
        else:
            param_name=param
            dist_above_bottom=-9999.
        
        if   param in ['p','q']:
            param_unit='m'  
        else:
            param_unit='m/s'
        
        p0 = outnc.createVariable(param_name,'f8',('time',))
        p0.units = param_unit
        p0.missing_value = -9999.0
        p0[:]=data

        p2 = outnc.createVariable('lat','f8')
        p2.units = 'deg'
        p2.missing_value = -9999.0
        p2[:]=lat_current
        
        p3 = outnc.createVariable('lon','f8')
        p3.units = 'deg'
        p3.missing_value = -9999.0
        p3[:]=lon_current
    
        p4 = outnc.createVariable('dist_above_bottom','f8')
        p4.units = 'm'
        p4.missing_value = -9999.0
        p4[:]=dist_above_bottom

        p5 = outnc.createVariable('depth','f8')
        p5.units = 'm'
        p5.missing_value = -9999.0
        p5.att = 'From 10m resolution 2011 map'
        p5[:]=depth_current
        
        
        outnc.history  = 'moghimis@gmail.com  steve2netCDF   '+datetime.datetime.now().isoformat()
        outnc.att=  infile[-12:]
        outnc.close()
        
        
    id_all=stations.keys()
    id_all.sort()
    for id1 in id_all:
        for param1 in ['u','v','w','p','q']:
            dircdf_cat=dircdf+'../cat'+dir[-4:]+'/'
            comm0='mkdir -p '+dircdf_cat
            comm1='ncrcat -O  '+dircdf+ '/*.'+param1+id1+'.nc  '+dircdf_cat+dir[-4:]+param1+id1+'.nc'
            os.system(comm0) 
            os.system(comm1)
            
             
        
