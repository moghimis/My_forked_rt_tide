from numpy import cos, sin, pi, absolute, arange
from scipy.signal import kaiserord, lfilter, firwin, freqz
from matplotlib.pylab import *

import datetime
import netcdftime
import netCDF4
import glob

#------------------------------------------------
# Create a signal for demonstration.
#------------------------------------------------

#sample_rate = 100.0  # Hz
#nsamples = 400       # total number
#t = arange(nsamples) / sample_rate
#x = cos(2*pi*0.5*t) + 0.2*sin(2*pi*2.5*t+0.1) + \
#        0.2*sin(2*pi*15.3*t) + 0.1*sin(2*pi*16.7*t + 0.1) + \
#            0.1*sin(2*pi*23.45*t+.8)


dr_plot=True
###############################################################
#If you like to cleaning which is arranged based on jun 
#rev of the data to be applied you need to change this to True
Sa_clean=True
###############################################################


n=1
every= n *(60 * 2)   # to skip every n min  for 2 Hz data sampling  


#tmp="/home/server/pi/homes/moghimi/Desktop/work/00-projs/01-muri/01-data/03-wave-data/00-NRI-specific/steve/data/RIVET_Britt_ADV_cdf/final_cat/*"
tmp="/home/server/pi/homes/moghimi/work/00-projs/01-muri/01-data/03-wave-data/00-NRI-specific/steve/data/rev02_oct2012/RIVET_Britt_ADV_cdf/final_cat/*"

#tmp="/home/server/pi/homes/moghimi/local/working/tmp/final_cat/*1.nc"
flist=glob.glob(tmp)
flist.sort()

for file in flist[:]:
    print file[-16:]
    print '>  reading signal ..'
    filtering='low'
    sample_rate = 2.0  # Hz
    obsncfile=file
    ncobs=netCDF4.Dataset(obsncfile)
    ncvarobs=ncobs.variables    
    timen='time'
    utime_obs=netcdftime.utime(ncvarobs[timen].units)
    #dates=utime_obs.num2date(ncvarobs[timen][:])
    param_name=ncvarobs.items()[1][0]
    param_namef=param_name+'_f'
    param_unit=ncvarobs[param_name].units
    
    elevobs=ncvarobs[param_name][:].squeeze()
    time=ncvarobs['time'][:].squeeze()
    
    x= 1. * elevobs
    t= 1. * time
    nsamples = len(x)
    
    print '>  doing filter ..'

    
    if filtering=='band':
       # total number
        
        
        #------------------------------------------------
        # Create a simple filter and apply it to x.
        #------------------------------------------------
        
        # The Nyquist rate of the signal.
        nyq_rate = sample_rate / 2.0
        
        #numtaps : int
        #    Length of the filter (number of coefficients, i.e. the filter
        #    order + 1).  `numtaps` must be even if a passband includes the
        #    Nyquist frequency.
        N=101
        
        # The cutoff frequency of the filter.
        cutoff_hz =   1/(3600.)
        low_limit_hz= 1/(48 * 3600.)
        
        taps = firwin(N, [low_limit_hz, cutoff_hz], pass_zero=False)
        # Use lfilter to filter x with the FIR filter.
        filtered_x = lfilter(taps, 1.0, x)
    elif filtering=='low':
        FS = sample_rate   #Hz                 # sampling rate
        FC = (1/60.)/(0.5*FS)                               # cutoff frequency at 1/60 Hz
        N = 101                                             # number of filter taps
        a = 1                                               # filter denominator
        b = firwin(N, cutoff=FC, window='hamming')          # filter numerator
        
        M = len(x)                                         # number of samples (60 seconds)
        n = arange(M)                                       # time index
        filtered_x = lfilter(b, a, x)

    

    ################################  Got Ride of Some part of the data  #####################
    if Sa_clean:
        if file[-16:-13]=='q00':
            print 'improve   q00'
    
            n1=976000
            n2=985000
            
            mean_part1=filtered_x[:n1].mean()
            mean_part2=filtered_x[n2:].mean()
            filtered_x[:n1]=filtered_x[:n1]+mean_part2-mean_part1
            
            filtered_x[n1:n2]=-999.99
            filtered_x[:1000]   =-999.99
            filtered_x[4750000:]=-999.99    
            filtered_x = ma.masked_where(filtered_x==-999.99,filtered_x)
            
        
        
        if file[-16:-13]=='q04':
            print 'improve   q04'
            n1=939800
            n2=944300
            filtered_x[n1:n2]   =-999.99
    
            ##################################################
            n3=4.94e6
            filtered_x[n3:]   =-999.99
            ##################################################
            n3=1000
            filtered_x[:n3]   =-999.99
            ##################################################        
            n1=976800
            n2=977600
            filtered_x[n1:n2]   =-999.99
            
            filtered_x = ma.masked_where(filtered_x==-999.99,filtered_x)
            
            mean_part1=filtered_x[:n1].mean()
            mean_part2=filtered_x[n2:].mean()
            filtered_x[:n1]=filtered_x[:n1]+4.6-3.9
    
        
        if file[-16:-13]=='q06':
            print 'improve   q06'
    
            ##################################################
            n3=4.64e6
            filtered_x[n3:]   =-999.99
            ##################################################
            n3=500
            filtered_x[:n3]   =-999.99
            ##################################################        
            n1=820000
            n2=940000
            filtered_x[n1:n2]   =-999.99
            
            filtered_x = ma.masked_where(filtered_x==-999.99,filtered_x)
            
            mean_part1=filtered_x[:n1].mean()
            mean_part2=filtered_x[n2:].mean()
            filtered_x[:n1]=filtered_x[:n1]+mean_part2-mean_part1
    
        
        if file[-16:-13]=='q52':
            print 'improve   q52'
    
            ##################################################
            n3=4972000
            filtered_x[n3:]   =-999.99
            ##################################################
            n3=500
            filtered_x[:n3]   =-999.99
            ##################################################        
            n1=1.2418e6+100
            n2=1.2418e6+450
            filtered_x[n1:n2]   =-999.99
            
            filtered_x = ma.masked_where(filtered_x==-999.99,filtered_x)
            
            mean_part1=filtered_x[:n1].mean()
            mean_part2=filtered_x[n2:].mean()
            filtered_x[:n1]=filtered_x[:n1]+(3.74-3.265)
        
        if file[-16:-13]=='q53':
            print 'improve   q53'
    
            ##################################################
            n3=5e6
            filtered_x[n3:]   =-999.99
            ##################################################
            n3=500
            filtered_x[:n3]   =-999.99
            ##################################################        
            n1=1.2404e6+200
            n2=1.2404e6+600
            filtered_x[n1:n2]   =-999.99
            
            filtered_x = ma.masked_where(filtered_x==-999.99,filtered_x)
            
            mean_part1=filtered_x[:n1].mean()
            mean_part2=filtered_x[n2:].mean()
            filtered_x[:n1]=filtered_x[:n1] +(3.75-3.15)    
        
        if file[-16:-13]=='q56':
            print 'improve   q56'
    
            ##################################################
            n3=5e6
            filtered_x[n3:]   =-999.99
            ##################################################
            n3=500
            filtered_x[:n3]   =-999.99
            ##################################################        
            n1=1.231e6+500
            n2=1.231e6+4000
            filtered_x[n1:n2]   =-999.99
            
            filtered_x = ma.masked_where(filtered_x==-999.99,filtered_x)
            
            mean_part1=filtered_x[:n1].mean()
            mean_part2=filtered_x[n2:].mean()
            filtered_x[:n1]=filtered_x[:n1] +(3.18-2.55)    
            x[:n1]=x[:n1] +(3.18-2.55)   
        
        
    
    #########################################################################################   
    
    
    print '>  create new Netcdf file ..'
    
    #file name
    
    cdfname=obsncfile[:-3]+'_filter_abs.nc'
    outnc=netCDF4.Dataset(cdfname,'w',format='NETCDF3_CLASSIC')
    outnc.createDimension('time',None)
    timea = outnc.createVariable('time','f8',('time',))
    timea.units = ncvarobs[timen].units
    timea[:]=t[::every]
    p0 = outnc.createVariable(param_name,'f8',('time',))
    p0.units = param_unit
    p0.missing_value = -999.99
    p0[:]=x[::every]
    
    p2 = outnc.createVariable('lat','f8')
    p2.units = 'deg'
    p2.missing_value = -9999.0
    p2[:]=ncvarobs['lat'][0]
    
    p3 = outnc.createVariable('lon','f8')
    p3.units = 'deg'
    p3.missing_value = -9999.0
    p3[:]=ncvarobs['lon'][0]
    
    p4 = outnc.createVariable('dist_above_bottom','f8')
    p4.units = 'm'
    p4.missing_value = -9999.0
    p4[:]=ncvarobs['dist_above_bottom'][0]
    
    p5 = outnc.createVariable('depth','f8')
    p5.units = 'm'
    p5.missing_value = -9999.0
    p5.att = 'From 10m resolution 2011 map'
    p5[:]=ncvarobs['depth'][0]

    p6 = outnc.createVariable(param_namef,'f8',('time',))
    p6.units = param_unit
    p6.missing_value = -9999.0
    p6[:]=filtered_x[::every]


    
    outnc.history  = 'moghimis@gmail.com  steve2netCDF_filtered   '+datetime.datetime.now().isoformat()
    outnc.att=  file
    outnc.close()
    ncobs.close()
    
   
    
    
    if dr_plot:
        #------------------------------------------------
        # Plot the FIR filter coefficients.
        #------------------------------------------------
        print '>  doing plots ..'
        if False:
            close('all')
            figure(1)
            plot(taps, 'bo', linewidth=2)
            title('Filter Coefficients (%d taps)' % N)
            grid(True)
            
            #------------------------------------------------
            # Plot the magnitude response of the filter.
            #------------------------------------------------
            
            figure(2)
            clf()
            w, h = freqz(taps, worN=8000)
            plot((w/pi)*nyq_rate, absolute(h), linewidth=2)
            xlabel('Frequency (Hz)')
            ylabel('Gain')
            title('Frequency Response')
            #ylim(-0.05, 1.05)
            grid(True)
            
            # Upper inset plot.
            ax1 = axes([0.42, 0.6, .45, .25])
            plot((w/pi)*nyq_rate, absolute(h), linewidth=2)
            #xlim(0,8.0)
            #ylim(0.9985, 1.001)
            grid(True)
            
            # Lower inset plot
            ax2 = axes([0.42, 0.25, .45, .25])
            plot((w/pi)*nyq_rate, absolute(h), linewidth=2)
            #xlim(12.0, 20.0)
            #ylim(0.0, 0.0025)
            grid(True)
            
            #------------------------------------------------
            # Plot the original and filtered signals.
            #------------------------------------------------
            
            # The phase delay of the filtered signal.
            delay = 0.5 * (N-1) / sample_rate
        
        fig3=figure(3,figsize=(30,8))
        # Plot the original signal.
        plot(x[::every],'b-', linewidth=0.25)
        # Plot thbe filtered signal, shifted to compensate for the phase delay.
        plot(filtered_x[::every], 'g-', linewidth=0.75)
        # Plot just the "good" part of the filtered signal.  The first N-1
        # samples are "corrupted" by the initial conditions.
        #plot(t[N-1:]-delay, filtered_x[N-1:], 'g', linewidth=4)
        #ylim(5, 7)
        
        xlabel('t')
        grid(True)
        pngfile='png/'+file[-16:-13]+'.png'
        savefig(pngfile,dpi=300)
        close('all')
        clf()
        #show()
