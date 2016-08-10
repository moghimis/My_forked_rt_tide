%-------------------------------------------------
function [fx,Pxr,Pxi,Pxc]= psd_spec(x,varargin)
    % RESIDUAL_SPECTRUM: Computes statistics from an input spectrum over
    % a number of bands, returning the band limits and the estimates for
    % power spectra for real and imaginary parts and the cross-spectrum.          
    %
    % Mean values of the noise spectrum are computed for the following 
    % 8 frequency bands defined by their center frequency and band width:
    % M0 +.1 cpd; M1 +-.2 cpd; M2 +-.2 cpd; M3 +-.2 cpd; M4 +-.2 cpd; 
    % M5 +-.2 cpd; M6 +-.21 cpd; M7 (.26-.29 cpd); and M8 (.30-.50 cpd). 

    % S. Lentz  10/28/99
    % R. Pawlowicz 11/1/00
    % Version 1.0
    
    % ###############################################################
    % Setup code for saving test data
    global psd_spec_saved
    savetestdata = false;
    if nargin > 1
        varargs = parse_varargin(varargin);
        f = fieldnames(varargs);
        if ~isempty(intersect('savetestdata',f)), savetestdata = varargs.savetestdata;  end;
         invals.x = x;
         invals.varargs = varargs;    invals.varargs.savetestdata = false;
    end
    % ################################################################
    dt = 1;
    %n = length(x);
    
    % Spectral estimate (takes real time series only).
    idx = find(isfinite(x));
    n = length(idx);
    [Pxr,fx]=psd(real(x(idx)),n,1/dt); %#ok<NASGU> % Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme. If you have an error here you are probably missing this toolbox
    [Pxi,fx]=psd(imag(x(idx)),n,1/dt); %#ok<NASGU> % Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme.
    try
        [Pxc,fx]=cpsd(real(x(idx)),imag(x(idx)),n,1/dt); % Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme.
    catch
       %junk = 0; 
       rethrow(lasterror)
    end
    
    Pxr = (2/n)*Pxr;
    Pxi = (2/n)*Pxi;
    Pxc = (2/n)*Pxc;
   
    % ################################################################
    % Code for saving test data
    if savetestdata && ~psd_spec_saved
        psd_spec_saved = true;
        outvals.fx = fx;
        outvals.Pxr = Pxr;
        outvals.Pxi = Pxi;
        outvals.Pxc = Pxc;
        writetestdata('psd_spec',invals,outvals);
    end
    % ###############################################################
end

