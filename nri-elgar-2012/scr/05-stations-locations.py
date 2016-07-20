import pygmaps
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

for i in range(31):
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

fpt.close()

klm_pre='inp/pre_kml.inp'
f_klm_pre=open(klm_pre, "r")
content = f_klm_pre.readlines()
f_klm_pre.close()


klm_out='elgar_adv.kml'
fout=open(klm_out, "w")
fout.writelines(content)

id_all=stations.keys()
id_all.sort()
for id1 in id_all:
    lat1=str(stations[id1][0])
    lon1=str(stations[id1][1])
    title=''+id1
    
    
    
    #fout.write('   <name>Elgar data </name> \n')
    #fout.write('   <open>1</open>\n')
    fout.write('    <Placemark>\n')
    name=      '       <name>'+title+'</name>\n'
    fout.write(name)
    fout.write('        <LookAt>\n')
    lon_info=  '            <longitude>'+lon1+'</longitude>\n'
    lat_info=  '            <latitude>'+lat1+'</latitude>\n'
    fout.write(lon_info)
    fout.write(lat_info)
    fout.write('           <altitude>0</altitude>\n')
    fout.write('           <range>2866.08080091982</range>\n')
    fout.write('           <tilt>6.652947352409956</tilt>\n')
    fout.write('           <heading>-0.03141501526757225</heading>\n')
    fout.write('           <altitudeMode>relativeToGround</altitudeMode>\n')
    fout.write('            <gx:altitudeMode>relativeToSeaFloor</gx:altitudeMode>\n')
    fout.write('       </LookAt>\n')
    fout.write('      <styleUrl>#msn_shaded_dot</styleUrl>\n')
    #fout.write('      <styleUrl>#msh_placemark_circle_highlight</styleUrl>\n')
    fout.write('      <Point>\n')
    cord_info= '           <coordinates>'+lon1+','+lat1+','+'0</coordinates>\n'
    fout.write(cord_info)
    fout.write('       </Point>\n')
    fout.write('    </Placemark>\n')
fout.write('     </Folder>\n')
fout.write(' </Document>\n')
fout.write(' </kml>  \n ') 
    
fout.close()   



if False:
    #Somehow worksfor googlemaps
    
    ########## CONSTRUCTOR: pygmaps.maps(latitude, longitude, zoom) ##############################
    # DESC:         initialize a map  with latitude and longitude of center point  
    #               and map zoom level "15"
    # PARAMETER1:   latitude (float) latittude of map center point
    # PARAMETER2:   longitude (float) latittude of map center point
    # PARAMETER3:   zoom (int)  map zoom level 0~20
    # RETURN:       the instant of pygmaps
    #========================================================================================
    mymap = pygmaps.maps( 34.539667, -77.339933, 18)
    
    id_all=stations.keys()
    id_all.sort()
    for id1 in id_all:
        ########## FUNCTION:  addpoint(latitude, longitude, [color])#############################
        # DESC:         add a point into a map and dispaly it, color is optional default is red
        # PARAMETER1:   latitude (float) latitude of the point
        # PARAMETER2:   longitude (float) longitude of the point
        # PARAMETER3:   color (string) color of the point showed in map, using HTML color code
        #               HTML COLOR CODE:  http://www.computerhope.com/htmcolor.htm
        #               e.g. red "#FF0000", Blue "#0000FF", Green "#00FF00"
        # RETURN:       no return
        #========================================================================================
        lat1=stations[id1][0]
        lon1=stations[id1][0]
        title='Elg'+id1
        mymap.addpoint(lat1, lon1, "#0000FF") #, title)
    
    mymap.draw('./mymap.html')


 


