function [ fband,Pxrave,Pxiave,Pxcave ] = fband_avg(fx,fu,Pxr,Pxi,Pxc,varargin )
    
    % ###############################################################
    % Setup code for saving test data.
    global fband_avg_saved
    savetestdata = false;
    if nargin > 5
        varargs = parse_varargin(varargin);
        
        f = fieldnames(varargs);
        if ~isempty(intersect('savetestdata',f)), savetestdata = varargs.savetestdata;  end;
        invals.fx = fx;
        invals.fu = fu;
        invals.Pxr = Pxr;
        invals.Pxi = Pxi;
        invals.Pxc = Pxc;
        invals.varargs = varargs; invals.varargs.savetestdata = false;
    end
    % ##############################################################
    
   % Mean values of the noise spectrum are computed for the following 
    % 8 frequency bands defined by their center frequency and band width:
    % M0 +.1 cpd; M1 +-.2 cpd; M2 +-.2 cpd; M3 +-.2 cpd; M4 +-.2 cpd; 
    % M5 +-.2 cpd; M6 +-.21 cpd; M7 (.26-.29 cpd); and M8 (.30-.50 cpd). 

    df=fx(3)-fx(2);
    if ~isempty(Pxr),   Pxr(round(fu./df)+1)=NaN; end; % Sets Px=NaN in bins close to analyzed frequencies
    if ~isempty(Pxi),   Pxi(round(fu./df)+1)=NaN; end; % (to prevent leakage problems?).
    if ~isempty(Pxc),   Pxc(round(fu./df)+1)=NaN; end;

    % Loop downwards in frequency through bands (cures short time series
    % problem with no data in lowest band).
    %
    % Divide by nx to get power per frequency bin, and multiply by 2
    % to account for positive and negative frequencies.
    % Define frequency bands for spectral averaging.
    fband =[.00010 .00417;
            .03192 .04859;
            .07218 .08884;
            .11243 .12910;
            .15269 .16936;
            .19295 .20961;
            .23320 .25100;
            .26000 .29000;
            .30000 .50000];
    
    % If we have a sampling interval> 1 hour, we might have to get
    % rid of some bins.
    %fband(fband(:,1)>1/(2*dt),:)=[];

    nfband=size(fband,1);
    
    Pxrave=zeros(nfband,1);
    Pxiave=zeros(nfband,1);
    Pxcave=zeros(nfband,1);

    for k=nfband:-1:1,
       b1 = find(fx>=fband(k,1));
       b2 = find(fx<=fband(k,2));
       b3 = find(isfinite(Pxr));

       jband=intersect(intersect(b1,b2),b3);

       if any(jband),
         if ~isempty(Pxr), Pxrave(k)=mean(Pxr(jband)); end;
         if ~isempty(Pxi),  Pxiave(k)=mean(Pxi(jband)); end;
         if ~isempty(Pxc),  Pxcave(k)=mean(Pxc(jband)); end;
       elseif k<nfband,
         Pxrave(k)=Pxrave(k+1);   % Low frequency bin might not have any points...
         Pxiave(k)=Pxiave(k+1);   
         Pxcave(k)=Pxcave(k+1);   
       end;
    end

    % #############################################################
    % Code to save test data
    if savetestdata && ~fband_avg_saved
        fband_avg_saved = true;
        outvals.fband = fband;
        outvals.Pxrave = Pxrave;
        outvals.Pxiave = Pxiave;
        outvals.Pxcave = Pxcave;
        writetestdata('fband_avg',invals,outvals);
    end
    % #############################################################
end       