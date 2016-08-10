Hi all,

I've fixed another bug in r_t_tide having to do with the nodal correction in the prediction file, i.e., r_t_predict.   Basically, nodal corrections were not being applied due to issues with passing the user preference or flag (e.g., 'nodalcorrflag'  'true') into the program.  It did work sometimes if one passed the flag in a totally undocumented way.  Why it was done this way, have no idea.

Anyhow, now it should be fixed.  See attached.  To be sure that a nodal correction is happening, the nodal factor applied to each constituent is now output to the screen.  This is interesting, among other things, because it makes the inner workings of the program less opaque.  You can now see, for example, that the M2 and K1 corrections are out of phase, and that the magnitude of K1 change is much more than M2.  Interestingly, the program applies corrections to the shallow water constituents.  Since these constituents are not astronomic, there are likely some major assumptions behind this correction. In other words, handle with care and don't believe the correction unless you validate it.

If you are doing a prediction for a long time period, the other important thing to remember is that r_t_tide only applies one nodal correction for each prediction.  This correction is applied using the correction from the middle of the data set.   In other words, use only 1 year chunks of data (or smaller) in your predictions, and increment forward year by year.   A prediction based off a 20 year data set of time will not have the correct nodal correction.

thanks,

Stefan
