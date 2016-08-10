function y=r_t_predict(tim,names,freq,tidecon,msl,nodalcorrflag,greenwichcorrflag,varargin)
% r_t_predic Tidal prediction
% YOUT=r_t_predic(TIM,NAMES,FREQ,TIDECON,MSL) makes a tidal prediction
% using the output of r_t_tide at the specified times TIM in decimal 
% days (from DATENUM). Optional arguments can be specified using
% property/value pairs: 
%
%       YOUT=r_t_predic(...,TIDECON,property,value,...)
%
% Available properties are:
%
%    In the simplest case, the tidal analysis was done without nodal
%    corrections, and thus neither will the prediction. If nodal 
%    corrections were used in the analysis, then it is likely we will
%    want to use them in the prediction too and these are computed 
%    using the latitude.
%
%     'latitude'        decimal degrees (+north) (default: none)
%
%    The tidal prediction may be restricted to only some of the 
%    available constituents:
%
%     'synthesis'    0 - Use all selected constituents.  (default)
%                    scalar>0 - Use only those constituents with a SNR
%                               greater than that given (1 or 2 are
%                               good choices).
%
%
%  It is possible to call r_t_predic without using property names, in
%  which case the assumed calling sequence is
%
%    YOUT=r_t_predic(TIM,NAMES,FREQ,TIDECON,LATITUDE,SYNTHESIS);
%
%  r_t_predic can be called using the tidal structure available as an 
%  optional output from r_t_tide
%
%    YOUT=r_t_predic(TIM,TIDESTRUC,...)
%
% R. Pawlowicz 11/8/99
% Version 1.0

% Do the synthesis.        

synth = 3;
lat = [];
% #########################################################################
% Setup code for saving test data
savetestdata = false;
global r_t_predict_saved
if nargin > 7
    varargs = parse_varargin(varargin);
    f = fieldnames(varargs);
    if ~isempty(intersect('synth',f)),  synth = varargs.synth;  end;
    if ~isempty(intersect('lat',f)),    lat = varargs.lat;  end;
    if ~isempty(intersect('savetestdata',f)), savetestdata = varargs.savetestdata;  end;
    %tim,names,freq,tidecon,msl,nodalcorrflag,greenwichcorrflag
    invals.tim = tim;
    invals.names = names;
    invals.freq = freq;
    invals.tidecon = tidecon;
    invals.msl = msl;
    invals.nodalcorrflag = nodalcorrflag;
    invals.greenwichcorrflag = greenwichcorrflag;
    invals.varargs = varargs;   invals.varargs.savetestdata = false;
end
% #########################################################################


snr=(tidecon(:,1)./tidecon(:,2)).^2;  % signal to noise ratio
if synth>0,
   I=snr>synth;
   if ~any(I),
     warning('No predictions with this SNR');
     y=NaN+zeros(size(tim));
     return;
   end;  
   tidecon=tidecon(I,:);
   names=names(I,:);
   freq=freq(I);  
end;    

    
if size(tidecon,2)==4,  % Real time series
  ap=tidecon(:,1)/2.*exp(-1i*tidecon(:,3)*pi/180);
  am=conj(ap);
else
  ap=(tidecon(:,1)+tidecon(:,3))/2.*exp( 1i*pi/180*(tidecon(:,5)-tidecon(:,7)));
  am=(tidecon(:,1)-tidecon(:,3))/2.*exp( 1i*pi/180*(tidecon(:,5)+tidecon(:,7)));
end;

% Mean at central point (get rid of one point at end to take mean of
% odd number of points if necessary).
jdmid=mean(tim(1:2*fix((length(tim)-1)/2)+1));

minres = 1/(tim(2)-tim(1));
const=load_constits(tim,0,minres,'savetestdata',savetestdata);
ju=zeros(size(freq));

% Check to make sure names and frequencies match expected values.

for k=1:size(names,1),
  ju(k)=strmatch(names(k,:),const.name);
end;
%if any(freq~=const.freq(ju)),
%  error('Frequencies do not match names in input');
%end;

% Get the astronical argument with or without nodal corrections.
if  ~isempty(lat) && abs(jdmid)>1,%nodalcorrflag && greenwichcorrflag %
  [v,u,f]=r_t_vuf(tim,jdmid,ju,'lat',lat,'savetestdata',savetestdata);
  disp('Stefan, a nodal correction has occurred');
elseif abs(jdmid)>1, % a real date greenwichcorrflag %
  [v,u,f]=r_t_vuf(tim,jdmid,ju,'lat',lat,'savetestdata',savetestdata);
else
   v=zeros(length(ju),1);
   u=v;
   f=ones(length(ju),1);  
end;

ap=ap.*f.*exp(+1i*2*pi*(u+v));
am=am.*f.*exp(-1i*2*pi*(u+v));


tim=tim-jdmid;

[n,m]=size(tim);
tim=tim(:)';

y=sum(exp( 1i*2*pi*freq*tim*24).*ap(:,ones(size(tim))),1)+ ...
     sum(exp(-1i*2*pi*freq*tim*24).*am(:,ones(size(tim))),1);

y = y + msl;     
y=reshape(y,n,m);

    % #######################################################################
    % Save test data
    if savetestdata && ~r_t_predict_saved
        r_t_predict_saved = true;
        outvals.y = y;   
        writetestdata('r_t_predict',invals,outvals);    
    end
    % #######################################################################
end

  
