See attached files.  May not be all the ones needed.  

The file tides_v_river defines a for-loop for each of the different data sites, sequentially.  Would probably have been good to define a function, but for the quick answer a cut and paste technique seemed faster.  Hence the long file.







	
Attachments9:37 AM (13 minutes ago)
		
I've attached a copy of Stefan's r_t_tide fix to run a limited constituent set.  Place the the script in the Matlab directory with r_t_tide.  Use the argument 'con' and provided an array with the constituent set you would like to use (code and matlab file). Below is the command needed to run the script.  



[name,freq,tidecon,xout]=  r_t_tide_sat(t,h, 'latitude',46.2,  'nodalcorrflag','true','greenwichcorrflag''true','start time',tm(1), 'con', lmt_con);
    
















Thanks--good idea to send out.  Actually, think you sent out your own version of r_t_tide, though the command for calling the 'constituent' flag is correct.  The (slight) problem with the version you sent is that the 'constituent' flag is not in the header/help section.  The other issues is that the problem with the Rayleigh flag is still there.

So that there is no confusion, I've attached the full version of r_t_tide that I've been using.  The main difference with what you sent out (I think) is that I've included directions for both the 'constit' flag and the 'rayeigh' flag in the header.  And, I've fixed the error that happened when the Rayleigh flag was called.  There is also an example of using a kaiser filter in the file labeled 'r_t_tide_kaiser'.   It's definitely a hack, but works.

One of the troubles with r_t_tide is that there are many versions floating around, some of which don't work as intended.  This is hopefully a clean version--In the past I've checked the nodal correction functionality to make sure it worked, and compared against Ed's HA program. 

Here's an example of the format I use to apply a nodal correction.  It may not be necessary to include the start time, but (if memory serves) I included it to make sure that r_t_tide went into the correct part of the program (in t_tide, start time is necessary).   But it's possible that this is just a work-around left over from an early, un-fixed version of the program.  To make sure that nodal correction is working, or even just get a feel for what it is doing, it's a good idea to try running your program with and without the correction (or even programming your own simple HA program). 

[nameu,fu,tidecon,xout]=r_t_tide(t,WL,'starttime',stime,'latitude',latitude,'nodalcorrflag','true','greenwichcorrflag','true','method','cauchy');

One other note--I have an experimental, 'hacked' version of r_t_tide that alters the basis functions to include a derivative constraint (i.e., constrains the slope at high and low water to be zero).  If anyone is interested, send me an email and I'll pull it together. 


thanks
