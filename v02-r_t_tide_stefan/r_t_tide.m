function [nameu,fu,tidecon,xout,stats]=r_t_tide(t,xin,varargin)

% EXPERIMENTAL VERSION !!  K. Leffler 3 February, 2009.

% This version uses Matlab's implementation of iteratively reweighted 
% least squares to reduce the influence of high-residual data.
% 
% r_t_tide Harmonic analysis of a time series
% [NAME,FREQ,TIDECON,XOUT]=r_t_tide(t,XIN) computes the tidal analysis
% of the (possibly complex) time series XIN, recorded at time t.

% Further inputs are optional, and are specified as property/value pairs
% [...]=r_t_tide(XIN,property,value,property,value,...,etc.)
%
% These properties are:

%       'method'        Robust fit weighting method.  This can be any of
%                       methods supported by the robustfit method.  Default
%                       value is 'ols' (Ordinary Least Squares), the method
%                       implemented by t-tide.
%
%       'interval'       Sampling interval (hours), default = 1.
%
%   The next two are required if nodal corrections are to be computed,
%   otherwise not necessary. If they are not included then the reported
%   phases are raw constituent phases at the central time.
%       'start time'     [year,month,day,hour,min,sec]
%                        - min,sec are optional OR
%                        decimal day (matlab DATENUM scalar)
%
%                       KEL : start time is depracated in r_t_tide.

%       'latitude'       decimal degrees (+north) (default: none).
%
% IMPORTANT : K. Leffler finds the "behavior by inference" by passing in
% the start time and latitude, or not, to be really annoying.  Start time
% is depracated by the requirement of passing in the vector of observed
% times.  The following parameters have been added :
%      
%      'greenwichcorrflag' :  A boolean (true/false) flag indicating
%           whether to standardize the phase.  Default value is false
%
%      'nodalcorrflag'     : A boolean (true/false) flag indicating 
%           whether to compute nodal corrections.  Default value is false.

%   Where to send the output.
%       'output'         where to send printed output:
%                        'none'    (no printed output)
%                        'screen'  (to screen) - default
%                        FILENAME   (to a file)
%
%   Correction factor for prefiltering.
%       'prefilt'        FS,CORR
%                        If the time series has been passed through
%                        a pre-filter of some kind (say, to reduce the
%                        low-frequency variability), then the analyzed
%                        constituents will have to be corrected for
%                        this. The correction transfer function
%                        (1/filter transfer function) has (possibly
%                        complex) magnitude CORR at frequency FS (cph).
%                        Corrections of more than a factor of 100 are
%                        not applied; it is assumed these refer to tidal
%                        constituents that were intentionally filtered
%                        out, e.g., the fortnightly components.
%
%   Inference of constituents.
%       'inference'      NAME,REFERENCE,AMPRAT,PHASE_OFFSET
%                        where NAME is an array of the names of
%                        constituents to be inferred, REFERENCE is an
%                        array of the names of references, and AMPRAT
%                        and PHASE_OFFSET are the amplitude factor and
%                        phase offset (in degrees)from the references.
%                        NAME and REFERENCE are Nx4 (max 4 characters
%                        in name), and AMPRAT and PHASE_OFFSET are Nx1
%                        (for scalar time series) and Nx2 for vector
%                        time series (column 1 is for + frequencies and
%                        column 2 for - frequencies).
%
%   Shallow water constituents
%       'shallow'        NAME
%                        A matrix whose rows contain the names of
%                        shallow-water constituents to analyze.
%
%   Resolution criterions for least-squares fit.
%       'rayleigh'       scalar - Rayleigh criteria, default = 1.
%                        cell array of strings - names of constituents to
%                                   use, e.g.,  constit_name{1} = 'M2'; constit_name{2} = 'S2'; and so on.).
%                          
%
%   Apply subset of constituents
%       'constit'   %                cell array of strings - names of constituents to
%                                   use, e.g.,  constit_name{1} = 'M2'; constit_name{2} = 'S2'; and so on.).  
%
%   Calculation of confidence limits.
%       'error'          'wboot'  - Boostrapped confidence intervals
%                                   based on a correlated bivariate
%                                   white-noise model.
%                        'cboot'  - Boostrapped confidence intervals
%                                   based on an uncorrelated bivariate
%                                   coloured-noise model (default).
%                        'linear' - Linearized error analysis that
%                                   assumes an uncorrelated bivariate
%                                   coloured noise model.
%
%   Computation of "predicted" tide (passed to r_t_predic, but note that
%                                    the default value is different).
%       'synthesis'      0 - use all selected constituents
%                        scalar>0 - use only those constituents with a
%                                   SNR greater than that given (1 or 2
%                                   are good choices, 2 is the default).
%                              <0 - return result of least-squares fit
%                                   (should be the same as using '0',
%                                   except that NaN-holes in original
%                                   time series will remain).
%
%   Robust fitting - added for version 1.03.  K. Leffler, 20 Feb, 2007.
%   
%       'method'        One of the valid weighting methods for robustfit
%                       The default is 'ols', which is ordinary least
%                       squares.
%
%       'tconst'        The tuning constant of associated with method.
%                       The %Time only default is 1, which is really only valid
%                       for 'ols'.  It's highly recommended that you 
%                       pass the proper value, found in the robustfit 
%                       documentation.
%
%       It is possible to call r_t_tide without using property names,
%       in which case the assumed calling sequence is
%
%          r_t_tide(XIN,INTERVAL,START_TIME,LATITUDE,RAYLEIGH)
%
%
%  OUTPUT:
%
%    nameu=list of constituents used
%    fu=frequency of tidal constituents (cycles/hr)
%    tidecon=[fmaj,emaj,fmin,emin,finc,einc,pha,epha] for vector xin
%           =[fmaj,emaj,pha,epha] for scalar (real) xin
%       fmaj,fmin - constituent major and minor axes (same units as xin)
%       emaj,emin - 95% confidence intervals for fmaj,fmin
%       finc - ellipse orientations (degrees)
%       einc - 95% confidence intervals for finc
%       pha - constituent phases (degrees relative to Greenwich)
%       epha - 95% confidence intervals for pha
%    xout=tidal prediction
%    retval.stats = output structure from robustfit (Matlab Statistics toolbox).
%      
%    Note on 'robustfit' : The index entry from robustfit doesn'data.t work
%       in R2006a (Windows).  Search for robustfit to get the proper help page.

% Note: Although missing data can be handled with NaN, it is wise not
%       to have too many of them. If your time series has a lot of
%       missing data at the beginning and/or end, then truncate the
%       input time series.  The Rayleigh criterion is applied to
%       frequency intervals calculated as the inverse of the input
%       series length.
%
% A description of the theoretical basis of the analysis and some
% implementation details can be found in:
%
% Pawlowicz, R., B. Beardsley, and S. Lentz, "Classical Tidal
%   "Harmonic Analysis Including Error Estimates in MATLAB
%    using T_TIDE", Computers and Geosciences, 2002.
%
% (citation of this article would be appreciated if you find the
%  toolbox useful).



% R. Pawlowicz 11/8/99 - Completely rewritten from the transliterated-
%                        to-matlab IOS/Foreman fortran code by S. Lentz
%                        and B. Beardsley.
%              3/3/00  - Redid errors to take into account covariances
%                        between u and v errors.
%              7/21/00 - Found that annoying bug in error calc!
%              11/1/00 - Added linear error analysis.
%              8/29/01 - Made synth=1 default, also changed behavior
%                        when no lat/time given so that phases are raw
%                        at central time.
%              9/1/01  - Moved some SNR code to r_t_predic.
%              9/28/01 - made sure you can'data.t choose Z0 as constituent.
%              6/12/01 - better explanation for variance calcs, fixed
%                        bug in typed output (thanks Mike Cook).
%
% Version 1.03

% K Leffler 

% Robust fit additions are described in 
%
% Leffler, K, DA Jay, "Enhancing tidal harmonic analysis: Robust (hybrid
% L1/L2) solutions." Submitted to Continental Shelf Research, Special
% Volume of the Physics of Estuaries and Coastal Seas Conference, Astoria,
% OR, USA, September 2006.

% Questions / comments on the robust fit modifications can be sent to 
% leffler@cecs.pdx.edu

%  --------------------------------------------------------------------
% Initialization stuff

% t = get(self,'time');
% xin = get(self,'height');

% badidx = union(find(isnan(t)),find(isnan(xin)));
% t(badidx) = [];
% xin(badidx) = [];

% ######################################################################
% Setup code for saving test data
savetestdata = false;
global r_t_tide_saved 
if nargin > 2
    varargs = parse_varargin(varargin);
    f = fieldnames(varargs);
    if ~isempty(intersect('savetestdata',f)), savetestdata = varargs.savetestdata;  end;
    %t,xin,varargin
    invals.t = t;
    invals.xin = xin;
    invals.varargs = varargs;   invals.varargs.savetestdata = false;
end
% #####################################################################
nameu = [];
fu = [];
tidecon = [];
xout = [];
retval.stats = []; 

[inn,inm]=size(xin);
if ~(inn==1 || inm==1), error('Input time series is not a vector'); end;

data.xin=xin(:); % makes xin a column vector
%clear('xin');
data.nobs=length(data.xin);

% --------------------------------------------------------------------
% Main calculation 

opt = r_t_init(varargin);     % Set options from inputs or default.
[data.nobsu,data.t,retval.centraltime,opt.stime] = r_t_get_standardtime(t,'savetestdata',savetestdata);
%self = set(self,'centraltime',retval.centraltime);
data.xin = data.xin(1:data.nobsu);

% -------Get the frequencies to use in the harmonic analysis-----------
%minres = opt.ray/(data.nobsu);
minres = opt.ray/(max(data.t)-min(data.t));
[nameu,fu,ju,namei,fi,jinf,jref]=r_t_get_constituents(t,minres,opt.constitnames,...
    opt.shallownames,opt.inf.iname,opt.inf.irefname,retval.centraltime,'savetestdata',savetestdata);

if isempty(nameu),  return; end;

mu=length(fu); % # base frequencies
mi=length(fi); %#ok<NASGU> % # inferred

% Find the good data points (here I assume that in a complex time
% series, if u is bad, so is v).

gd=find(isfinite(data.xin(1:data.nobsu)));
data.ngood=length(gd);
fprintf('   Points used: %d of %d\n',data.ngood,data.nobs)


%----------------------------------------------------------------------
% Now solve the analysis. Instead of solving
% for + and - frequencies using exp(i*f*data.t), I use sines and cosines to
% keep tc real.  If the input series is real, than this will
% automatically use real-only computation. However, for the analysis,
% it's handy to get the + and - frequencies ('ap' and 'am'), and so
% that's what we do afterwards.

lastwarn('');

% t = get(self,'time');
% h = get(self,'height');
% h_aux = get(self,'height_extrema'); 

h = xin;
h_aux = [];
[retval.basis retval.basis_extrema] = r_t_get_basis(t,retval.centraltime,fu,'savetestdata',savetestdata);

retval.msl = mean(h(~isnan(h)));
retval.mslerr = 0; 

h = h-retval.msl;


% robustfit requires the Statistics toolbox.  An error here might
% mean you don'data.t have the toolbox installed.

%TC = [taper(self,retval.basis);taper(self,retval.basis_extrema)];
%H = [taper(self,h);taper(self,h_aux)];
TC = retval.basis;
H = h;


[coef stats]=robustfit(TC,H,opt.method,opt.tconst,'on');

% if ~isempty(lastwarn)
%     % A warning was generated.  Return with all values empty.
%     % If you're getting a maximum iterations reached error, it is possible
%     % to increase the maximum iterations by changing
%     % the robustfit m file.  Follow the stack trace to get to the
%     % right line.  Default is 50.  150 or 250 is a better choice.
% 
%     return;
% end

coef(1) = [];


ap=(coef(1:mu)-1i*coef(mu+1:end))/2;  % a+ amplitudes
am=(coef(1:mu)+1i*coef(mu+1:end))/2;  % a- amplitudes

%----------------------------------------------------------------------
% Check variance explained (but do this with the original fit).

xout = retval.basis(1:data.nobsu,:)*coef;
xres=data.xin-(xout+retval.msl); % and the residuals!

if isreal(data.xin),    % Real time series
    data.varx=cov(data.xin(gd));
    data.varxp=cov(xout(gd));
    data.varxr=cov(xres(gd));
    fprintf('   percent of var residual after lsqfit/var original: %5.2f %\n',100*(data.varxr/data.varx));
else               % Complex time series
    data.varx=cov(real(data.xin(gd)));
    data.varxp=cov(real(xout(gd)));
    data.varxr=cov(real(xres(gd)));
    fprintf('   percent of X var residual after lsqfit/var original: %5.2f %\n',100*(data.varxr/data.varx));

    vary=cov(imag(data.xin(gd)));
    varyp=cov(imag(xout(gd))); %#ok<NASGU>
    varyr=cov(imag(xres(gd)));
    fprintf('   percent of Y var residual after lsqfit/var original: %5.2f %\n',100*(varyr/vary));
end;


%---------- Correct for prefiltering-----------------------------------

corrfac=interp1(opt.corr_fs,opt.corr_fac,fu);
% To stop things blowing up!
corrfac(corrfac>100 | corrfac <.01 | isnan(corrfac))=1;

ap=ap.*corrfac;
am=am.*conj(corrfac);


%---------------Nodal Corrections--------------------------------------
% Generate nodal corrections and calculate phase relative to Greenwich.
% Note that this is a slightly weird way to do the nodal corrections,
% but is 'traditional'.  The "right" way would be to change the basis
% functions used in the least-squares fit above.

if opt.nodalcorrflag && opt.greenwichcorrflag    % Get nodal corrections at midpoint time.
    [v,u,f]=r_t_vuf(t,retval.centraltime,[ju;jinf],'lat',opt.lat,'savetestdata',savetestdata);

    vu=(v+u)*360; % total phase correction (degrees)
    data.nodcor='Greenwich phase computed with nodal corrections applied to amplitude \n and phase relative to center time';
elseif opt.greenwichcorrflag
    % Get nodal corrections at midpoint time
    [v,u,f]=r_t_vuf(t,retval.centraltime,[ju;jinf],'savetestdata',savetestdata);
    vu=(v+u)*360; % total phase correction (degrees)
    data.nodcor='Greenwich phase computed, no nodal corrections';
else   % No time, no latitude
    vu=zeros(length(ju)+length(jinf),1);
    f=ones(length(ju)+length(jinf),1);
    data.nodcor='Phases at central time';
end
fprintf(['   ',data.nodcor,'\n']);

%---------------Inference Corrections----------------------------------
% Once again, the "right" way to do this would be to change the basis
% functions.
ii=find(isfinite(jref));
if ii,
    fprintf('   Do inference corrections\n');
    snarg=data.nobsu*pi*(fi(ii)   -fu(jref(ii)) )*opt.dt;
    scarg=sin(snarg)./snarg;

    if size(opt.inf.amprat,2)==1,    % For real time series
        pearg=     2*pi*(vu(mu+ii)-vu(jref(ii))+opt.inf.ph(ii))/360;
        pcfac=opt.inf.amprat(ii).*f(mu+ii)./f(jref(ii)).*exp(1i*pearg);
        pcorr=1+pcfac.*scarg;
        mcfac=conj(pcfac);
        mcorr=conj(pcorr);
    else                          % For complex time series
        pearg=     2*pi*(vu(mu+ii)-vu(jref(ii))+opt.inf.ph(ii,1))/360;
        pcfac=opt.inf.amprat(ii,1).*f(mu+ii)./f(jref(ii)).*exp(1i*pearg);
        pcorr=1+pcfac.*scarg;
        mearg=    -2*pi*(vu(mu+ii)-vu(jref(ii))+opt.inf.ph(ii,2))/360;
        mcfac=opt.inf.amprat(ii,2).*f(mu+ii)./f(jref(ii)).*exp(1i*mearg);
        mcorr=1+mcfac.*scarg;
    end;

    ap(jref(ii))=ap(jref(ii))./pcorr;   % Changes to existing constituents
    ap=[ap;ap(jref(ii)).*pcfac];        % Inferred constituents

    am(jref(ii))=am(jref(ii))./mcorr;
    am=[am;am(jref(ii)).*mcfac];

    fu=[fu;fi(ii)];
    nameu=[nameu;namei(ii,:)];
end;

% --------------Error Bar Calculations---------------------------------
%
% Error bar calcs involve two steps:
%      1) Estimate the uncertainties in the analyzed amplitude
%         for both + and - frequencies (i.e., in 'ap' and 'am').
%         A simple way of doing this is to take the variance of the
%         original time series and divide it into the amount appearing
%         in the bandwidth of the analysis (approximately 1/length).
%         A more sophisticated way is to assume "locally white"
%         noise in the vicinity of, e.g., the diurnal consistuents.
%         This takes into account slopes in the continuum spectrum.
%
%      2) Transform those uncertainties into ones suitable for ellipse
%         parameters (axis lengths, angles). This can be done
%         analytically for large signal-to-noise ratios. However, the
%         transformation is non-linear at lows SNR, say, less than 10
%         or so.
%

weight = stats.w;
retval.ci_resid = weight(1:length(xres)).*xres;

if strmatch(opt.errcalc(2:end),'boot'),
    fprintf('   Using nonlinear bootstrapped error estimates\n');

    % "noise" matrices are created with the right covariance structure
    % to add to the analyzed components to create 'nreal' REPLICATES.
    %

    nreal=300;             % Create noise matrices
    [NP,NM,retval.Pxr]=r_t_noise_realizations(data.t,retval.ci_resid,fu,nreal,opt.errcalc,opt.specmethod,'savetestdata',savetestdata);

    % All replicates are then transformed (nonlinearly) into ellipse
    % parameters.  The computed error bars are then based on the std
    % dev of the replicates.

    AP=ap(:,ones(1,nreal))+NP;        % Add to analysis (first column
    AM=am(:,ones(1,nreal))+NM;        % of NM,NP=0 so first column of
    % AP/M holds ap/m).

    ap=abs(AP);                       % Angle/magnitude form:
    am=abs(AM);
    epsp=angle(AP)*180/pi;
    epsm=angle(AM)*180/pi;
elseif strmatch(opt.errcalc,'linear'),
    fprintf('   Using linearized error estimates\n');
    %
    % Uncertainties in analyzed amplitudes are computed in different
    % spectral bands. Real and imaginary parts of the residual time series
    % are treated separately (no cross-covariance is assumed).
    %
    % Noise estimates are then determined from a linear analysis of errors,
    % assuming that everything is uncorrelated. This is OK for scalar time
    % series but can fail for vector time series if the noise is not
    % isotropic.

    [ercx,eicx]=r_t_noise_stats(retval.ci_resid,fu,opt.dt,'specmethod',opt.specmethod,'savetestdata',savetestdata);
    % Note - here we assume that the error in the cos and sin terms is
    % equal, and equal to total power in the encompassing frequency bin.
    % It seems like there should be a factor of 2 here somewhere but it
    % only works this way! <shrug>
    [emaj,emin,einc,epha]=r_t_errell(ap+am,1i*(ap-am),ercx,ercx,eicx,eicx,'savetestdata',savetestdata);

    epsp=angle(ap)*180/pi;
    epsm=angle(am)*180/pi;
    ap=abs(ap);
    am=abs(am);
else
    error(['Unrecognized type of error analysis: ''' opt.errcalc ''' specified!']);
end;

%-----Convert complex amplitudes to standard ellipse parameters--------

aap=ap./f(:,ones(1,nreal));	% Apply nodal correretval.xoutctions and
aam=am./f(:,ones(1,nreal));	% compute ellipse parameters.

fmaj=aap+aam;                   % major axis
fmin=aap-aam;                   % minor axis

gp=mod( vu(:,ones(1,nreal))-epsp ,360); % pos. Greenwich phase in deg.
gm=mod( vu(:,ones(1,nreal))+epsm ,360); %#ok<NASGU> % neg. Greenwich phase in deg.

finc= (epsp+epsm)/2;
finc(:,1)=mod( finc(:,1),180 ); % Ellipse inclination in degrees
% (mod 180 to prevent ambiguity, i.e.,
% we always ref. against northern
% semi-major axis.

finc=r_t_cluster(finc,180,'savetestdata',savetestdata); 	% Cluster angles around the 'true'
% angle to avoid 360 degree wraps.

%pha = gp+finc;
pha=mod( gp+finc ,360); 	% Greenwich phase in degrees.

pha=r_t_cluster(pha,360,'savetestdata',savetestdata);		% Cluster angles around the 'true' angle
% to avoid 360 degree wraps.

%----------------Generate 95% CI---------------------------------------
% For bootstrapped errors, we now compute limits of the distribution.
if strmatch(opt.errcalc(2:end),'boot'),
    % std dev-based estimates.
    % The 95% CI are computed from the sigmas
    % by a 1.96 fudge factor (infinite degrees of freedom).
    % emaj=1.96*std(fmaj,0,2);
    % emin=1.96*std(fmin,0,2);
    % einc=1.96*std(finc,0,2);
    % epha=1.96*std(pha ,0,2);
    % Median-absolute-deviation (MAD) based estimates.
    % (possibly more stable?)
    emaj=median(abs(fmaj-median(fmaj,2)*ones(1,nreal)),2)/.6375*1.96;
    emin=median(abs(fmin-median(fmin,2)*ones(1,nreal)),2)/.6375*1.96;
    einc=median(abs(finc-median(finc,2)*ones(1,nreal)),2)/.6375*1.96;
    epha=median(abs( pha-median( pha,2)*ones(1,nreal)),2)/.6375*1.96;
else
    % In the linear analysis, the 95% CI are computed from the sigmas
    % by this fudge factor (infinite degrees of freedom).
    emaj=1.96*emaj;
    emin=1.96*emin;
    einc=1.96*einc;
    epha=1.96*epha;
end;

if isreal(data.xin),
    tidecon=[fmaj(:,1),emaj,pha(:,1),epha];
else
    tidecon=[fmaj(:,1),emaj,fmin(:,1),emin, finc(:,1),einc,pha(:,1),epha];
end;

% Sort results by frequency (needed if anything has been inferred since
% these are stuck at the end of the list by code above).
if any(isfinite(jref)),
    [fu,I]=sort(fu);
    nameu=nameu(I,:);
    tidecon=tidecon(I,:);
end;

data.snr=(tidecon(:,1)./tidecon(:,2)).^2;  % signal to noise ratio

%--------Generate a 'prediction' using significant constituents----------
if opt.synth>=0,
    if (opt.nodalcorrflag) && (opt.greenwichcorrflag) %~isempty(opt.lat) && ~isempty(opt.stime),
        fprintf('   Generating prediction with nodal corrections, SNR is %f\n',opt.synth);
        va{1} = 'synth'; va{2} = opt.synth;
        va{3} = 'lat';   va{4} = opt.lat;
        xout=r_t_predict(t(1:data.nobsu),nameu,fu,tidecon,retval.msl,opt.nodalcorrflag,opt.greenwichcorrflag,'lat',opt.lat,'savetestdata',savetestdata);
    elseif opt.greenwichcorrflag %~isempty(opt.stime),
        fprintf('   Generating prediction without nodal corrections, SNR is %f\n',opt.synth);
        xout=r_t_predict(t(1:data.nobsu),nameu,fu,tidecon,retval.msl,opt.nodalcorrflag,opt.greenwichcorrflag,'savetestdata',savetestdata);
        %xout=r_t_predic(opt.stime+(0:data.nobsu-1)*opt.dt/24.0,nameu,fu,tidecon,'synth',opt.synth);
    else
        fprintf('   Generating prediction without nodal corrections, SNR is %f\n',opt.synth);
        xout=r_t_predict(t(1:data.nobsu),nameu,fu,tidecon,retval.msl,opt.nodalcorrflag,opt.greenwichcorrflag,'savetestdata',savetestdata);
        %xout=r_t_predic(data.t/24.0,nameu,fu,tidecon,'synth',opt.synth);
    end;
else
    fprintf('   Returning fitted prediction\n');
end;

%----------------------------------------------------------------------
% Check variance explained (but now do this with the synthesized fit).
xres=data.xin(:)-xout(:); % and the residuals!

%error;

if isreal(data.xin),    % Real time series
    data.varx=cov(data.xin(gd));data.varxp=cov(xout(gd));data.varxr=cov(xres(gd));
    fprintf('   percent of var residual after synthesis/var original: %5.2f %\n',100*(data.varxr/data.varx));
else               % Complex time series
    data.varx=cov(real(data.xin(gd)));data.varxp=cov(real(xout(gd)));data.varxr=cov(real(xres(gd)));
    fprintf('   percent of X var residual after synthesis/var original: %5.2f %\n',100*(data.varxr/data.varx));

    vary=cov(imag(data.xin(gd)));
    varyp=cov(imag(xout(gd))); %#ok<NASGU>
    varyr=cov(imag(xres(gd)));
    fprintf('   percent of Y var residual after synthesis/var original: %5.2f %\n',100*(varyr/vary));
end;


%-----------------Output results---------------------------------------

if (opt.fid > 0)
    retval.fu = fu;
    retval.nameu = nameu;
    retval.tidecon = tidecon;
    r_t_writeoutfile(opt,data,retval);
end;
% #########################################################################
% Save test data
if savetestdata && ~r_t_tide_saved
    r_t_tide_saved = true;
    %nameu,fu,tidecon,xout
    outvals.nameu = nameu;
    outvals.fu = fu;
    outvals.tidecon = tidecon;
    outvals.xout = xout;
    writetestdata('r_t_tide',invals,outvals);    
end
% #########################################################################


end
%xout=reshape(xout,inn,inm);

% switch nargout,
%     case {0,3,4}
%     case {1}
%         nameu = struct('name',nameu,'freq',fu,'tidecon',tidecon);
%     case {2}
%         xout=reshape(xout,inn,inm);
%         nameu = struct('name',nameu,'freq',fu,'tidecon',tidecon);
%         fu=xout;
% end;
% -------------------------------------------------
function opt =  r_t_init(args)

%args = args{1};
opt.constitnames=[];
opt.corr_fac=[1  1];
opt.corr_fs=[0 1e6];
opt.dt = 1;
opt.errcalc='cboot';
opt.fid=1;
opt.greenwichcorrflag = false;
opt.inf.iname=[];
opt.inf.irefname=[];
opt.lat=[];
opt.method = 'ols';
opt.nodalcorrflag = false;
opt.ray=1;
opt.run = [];
opt.savetestdata = false;
opt.shallownames=[];
opt.specmethod = 'psd';
opt.synth=2;
opt.tconst = 1;
opt.tseries = 'hourly';

k=1;
while k<=length(args)
    if ischar(args{k}),
        assert(ischar(args{k}),'Assertion error: r_t_tide:r_t_init, non-character argument');
        dk = 2;
        switch lower(args{k}(1:3)),           
            case 'con',
                assert(iscell(args{k+1}),'Next argument after "constit" needs to be a cell array');
                opt.constitnames = args{k+1};
            case 'err',
                assert(ischar(args{k+1}),'Next argument after "error" needs to be a string');
                opt.errcalc=args{k+1};  
            case 'gre'
                assert(ischar(args{k+1}), 'Next argument after "greenwichcorrflag" needs to be a string');
                if strcmpi(args{k+1},'true'),   opt.greenwichcorrflag = true;  end;
            case 'inf',
                assert(ischar(args{k+1}) && ischar(args{k+2}) && ischar(args{k+3}) && ischar(args{k+4}),'Next four arguments after "inference" need to be strings');
                opt.inf.iname=args{k+1};
                opt.inf.irefname=args{k+2};
                opt.inf.amprat=args{k+3};
                opt.inf.ph=args{k+4};
                dk = 5;
            case 'lat',
                assert(isnumeric(args{k+1}),'argument after "latitude" needs to be a number');
                opt.lat=args{k+1};
            case 'met',
                assert(ischar(args{k+1}),'argument after "method" needs to be a string');
                opt.method = args{k+1};      
                
                switch lower(opt.method)
                    case 'andrews'
                        opt.tconst = 1.339;
                    case 'bisquare'
                        opt.tconst = 4.685;
                    case 'cauchy'
                        opt.tconst  = 2.385;
                    case 'fair'
                        opt.tconst = 1.400;
                    case 'huber'
                        opt.tconst = 1.345;
                    case 'logistic'
                        opt.tconst = 1.205;
                    case 'ols'
                        opt.tconst = 1;
                    case 'talwar'
                        opt.tconst = 2.795;
                    case 'welsch'
                        opt.tconst = 2.985;
                end
                
            case 'nod'
                assert(ischar(args{k+1}), 'argument after "nodalcorrflag" needs to be a string');
                if strcmpi(args{k+1},'true'),   opt.nodalcorrflag = true;   end;
            case 'out',                
                assert(ischar(args{k+1}),'argument after "output" needs to be a string');
                switch args{k+1}
                    case 'none',
                        opt.fid=-1;
                    case 'screen',
                        opt.fid=1;
                    otherwise
                        opt.filen=args{k+1};
                        [opt.fid,msg]=fopen(opt.filen,'w');
                        if opt.fid==-1, error(msg); end;
                end;
            case 'pre',
                assert(isnumeric(args{k+1}) && isnumeric(args{k+2}),'Next two arguments after "pre" need to be numeric');
                dk = 3;
                opt.corr_fs=args{k+1};
                opt.corr_fac=args{k+2};
            case 'ray',
                assert(isnumeric(args{k+1}) || iscell(args{k+1}),'argument after "rayleight" needs to be numeric or a cell array');
                if isnumeric(args{k+1}),
                    opt.ray=args{k+1};
                    %default_constits = opt.constitnames,
                    %disp('pause'); pause;
                else
                    opt.constitnames=args{k+1};
                    if iscellstr(opt.constitnames),% This line obsolete.  SAT 3/2016; %opt.constitnames=char(opt.constitnames); ss =opt.constitnames, 
                    end;
                end;
            case 'sav',
                assert(isa(args{k+1},'logical'),'argument after "savetestdata" needs to be true/false');
                opt.savetestdata = args{k+1};
            case 'sha',
                assert(iscell(args{k+1}),'argument after "shallow" needs to be a cell array');
                opt.shallownames=args{k+1};
            case 'spe',
                assert(ischar(args{k+1}),'argument after "specmethod" needs to be string');
                opt.specmethod = args{k+1};
            case 'sta',
                assert(isnumeric(args{k+1}),'argument after "starttime" needs to be numeric');
                opt.stime=args{k+1};
                if length(opt.stime)>1,
                    opt.stime=[opt.stime(:)' zeros(1,6-length(opt.stime))];
                    opt.stime=datenum(opt.stime(1),opt.stime(2),opt.stime(3),opt.stime(4),opt.stime(5),opt.stime(6));
                end;
            case 'syn',
                assert(isnumeric(args{k+1}),'argument after "synth" needs to be numeric');
                opt.synth=args{k+1};
            case 'tco',
                assert(isnumeric(args{k+1}),'argument after "tconst" needs to be numeric');
                opt.tconst=args{k+1};
            case 'tse'
                assert(iscell(args{k+1}),'argument after "tseries" needs to be a cell array');
                opt.tseries = args{k+1};
        end;
    end;
    k=k+dk;
end;
    global fband_avg_saved load_constits_saved ls_spec_saved psd_spec_saved r_t_astron_saved r_t_cluster_saved r_t_equilib_saved
    global r_t_errel_saved r_t_noise_realizations_saved r_t_noise_stats_saved r_t_predict_saved r_t_synth_saved r_t_tide_saved
    global r_t_vuf_saved r_t_xstat_saved r_t_xtide_saved r_t_get_basis_saved r_t_get_constituents_saved r_t_get_standardtime
 
    fband_avg_saved = false;
    load_constits_saved = false;
    ls_spec_saved = false;
    psd_spec_saved = false;
    r_t_astron_saved = false;
    r_t_cluster_saved = false;
    r_t_equilib_saved = false;
    
    r_t_errel_saved  = false;
    r_t_get_basis_saved = false;
    r_t_get_constituents_saved = false;
    r_t_get_standardtime = false;
    
    r_t_noise_realizations_saved = false;
    r_t_noise_stats_saved =false;
    r_t_predict_saved =false;
    r_t_synth_saved =false;
    r_t_tide_saved = false;
    
    r_t_vuf_saved =false;
    r_t_xstat_saved =false;
    r_t_xtide_saved = false;

    
end
% -------------------------------------------------
function   r_t_writeoutfile(opts,data,retval)

if opts.fid>1,
    fprintf(opts.fid,'\n%s\n',['file name: ',opts.filen]);
elseif opts.fid==1,
    fprintf(opts.fid,'-----------------------------------\n');
end

if opts.fid>0,
    fprintf(opts.fid,'date: %s\n',date);
    fprintf(opts.fid,'nobs = %d,  ngood = %d,  record length (days) = %.2f\n',data.nobs,data.ngood,length(data.xin)/24);
    if ~isempty(opts.stime); fprintf(opts.fid,'%s\n',['start time: ',datestr(opts.stime)]); end
    fprintf(opts.fid,'rayleigh criterion = %.1f\n',opts.ray);
    fprintf(opts.fid,'%s\n',data.nodcor);
    %  fprintf(opts.fid,'\n     coefficients from least squares fit of x\n');
    %  fprintf(opts.fid,'\n tide    freq        |a+|       err_a+      |a-|       err_a-\n');
    %  for k=1:length(retval.fu);
    %    if ap(k)>eap(k) | am(k)>eam(k), fprintf('*'); else fprintf(' '); end;
    %    fprintf(opts.fid,'%s  %8.5f  %9.4f  %9.4f  %9.4f  %9.4f\n',retval.namu(k,:),retval.fu(k),ap(k),eap(k),am(k),eam(k));
    %  end
    fprintf(opts.fid,'\nx0= %.3g, x trend= %.3g\n',real(retval.msl),real(0.00));
    fprintf(opts.fid,['\nvar(x)= ',num2str(data.varx),'   var(xp)= ',num2str(data.varxp),'   var(xres)= ',num2str(data.varxr) '\n']);
    fprintf(opts.fid,'percent var predicted/var original= %.1f %\n',100*data.varxp/data.varx);

    if isreal(data.xin)
        fprintf(opts.fid,'\n     tidal amplitude and phase with 95% CI estimates\n');
        fprintf(opts.fid,'\ntide   freq       amp     amp_err    pha    pha_err     data.snr\n');
        for k=1:length(retval.fu);
            if data.snr(k)>opts.synth, fprintf(opts.fid,'*'); else fprintf(opts.fid,' '); end;
            fprintf(opts.fid,'%s %9.7f %11.7f %11.7f %8.2f %8.2f %8.2g\n',retval.nameu(k,:),retval.fu(k),retval.tidecon(k,:),data.snr(k));
        end
    else
        fprintf(opts.fid,'\ny0= %.3g, x trend= %.3g\n',imag(z0),imag(dz0));
        fprintf(opts.fid,['\nvar(y)= ',num2str(vary),'    var(yp)= ',num2str(varyp),'  var(yres)= ',num2str(varyr) '\n']);
        fprintf(opts.fid,'percent var predicted/var original= %.1f %\n',100*varyp/vary);
        fprintf(opts.fid,'\n%s\n','ellipse parameters with 95% CI estimates');
        fprintf(opts.fid,'\n%s\n','tide   freq      major  emaj    minor   emin     inc    einc     pha    epha      data.snr');
        for k=1:length(retval.fu);
            if data.snr(k)>opts.synth, fprintf(opts.fid,'*'); else fprintf(opts.fid,' '); end;
            fprintf(opts.fid,'%s %9.7f %6.3f %7.3f %7.3f %6.2f %8.2f %6.2f %8.2f %6.2f %6.2g\n',...
                retval.namu(k,:),retval.fu(k),retval.tidecon(k,:),data.snr(k));
        end
        fprintf(opts.fid,['\ntotal var= ',num2str(data.varx+vary),'   pred var= ',num2str(data.varxp+varyp) '\n']);
        fprintf(opts.fid,'percent total var predicted/var original= %.1f %\n\n',100*(data.varxp+varyp)/(data.varx+vary));
    end

    if opts.fid~=1, fclose(opts.fid); end
end
end