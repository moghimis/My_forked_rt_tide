   function [ercx,eicx]=r_t_noise_stats(xres,fu,dt,varargin)
    % NOISE_STATS: Computes statistics of residual energy for all 
    % constituents (ignoring any cross-correlations between real and
    % imaginary parts).

    % S. Lentz  10/28/99
    % R. Pawlowicz 11/1/00
    % Version 1.0
    specmethod = 'psd';
    
    % ####################################################################
    % Set up code for saving test data
    savetestdata = false;
    global r_t_noise_stats_saved
    
    if nargin>3
        varargs = parse_varargin(varargin);
        f = fieldnames(varargs);
        if intersect('specmethod',f),   specmethod = varargs.specmethod;end;
        if intersect('savetestdata',f), savetestdata = varargs.savetestdata;end;
        invals.xres = xres;
        invals.fu = fu;
        invals.dt = dt;
        invals.varargs = varargs;   invals.varargs.savetestdata = false;
    end
    % #####################################################################
    if strcmpi(specmethod,'ls')
      [fx,Pxr,Pxi,Pxc] = ls_spec(t,xres,'savetestdata',savetestdata);
    else
      [fx,Pxr,Pxi,Pxc] = psd_spec(xres,'savetestdata',savetestdata);
    end
    [fband,Pxrave,Pxiave,Pxcave] = fband_avg(fx,fu,Pxr,Pxi,Pxc,'savetestdata',savetestdata);
    Pxcave=zeros(size(Pxcave));  %% For comparison with other technique!
    %fprintf('**** Assuming no covariance between u and v errors!*******\n');
    
    nfband=size(fband,1);
    mu=length(fu);

    % Get the statistics for each component.
    ercx=zeros(mu,1);
    eicx=zeros(mu,1);
    for k1=1:nfband;
       k=find(fu>=fband(k1,1) & fu<=fband(k1,2));
       ercx(k)=sqrt(Pxrave(k1));
       eicx(k)=sqrt(Pxiave(k1));
    end
    % #####################################################################
    % Save test data
    if savetestdata && ~r_t_noise_stats_saved
        r_t_noise_stats_saved = true;
        %ercx,eicx
        outvals.ercx = ercx;
        outvals.eicx = eicx;

        writetestdata('r_t_noise_stats',invals,outvals);

    end
    % ####################################################################
   end