%function [ output_args ] = writetestdata( input_args )
function writetestdata(fname,invals,outvals)
%WRITETESTDATA Summary of this function goes here
%   Detailed explanation goes here
    
    outfile = ['~/psu/code/tidal_ha/r_t_tide/testdata/' fname '.mat'];    
    save(outfile,'invals','outvals');
    
end

