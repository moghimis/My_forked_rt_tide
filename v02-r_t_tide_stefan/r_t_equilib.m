function [name,freq,amp]=r_t_equilib(lat,varargin)
% r_t_equilib Equilibrium amplitude of the tidal potential
% [NAME,FREQ,AMPLITUDE]=r_t_equilib(LAT) returns vectors with the
% NAME of tidal constituents, their FREQ (in cph), and their
% equilibrium AMPLITUDE in the tidal potential as a function of
% LATitude (degrees). If LAT is a vector, then AMPLITUDE is a
% matrix in which each column corresponds to a specific latitude.
%
% If no output arguments are specified, the equilibrium spectrum
% is plotted.

% R. Pawlowicz 9/11/99
% Version 1.0
    
    % ###############################################################
    % Setup code for saving test data
    savetestdata = false;
    global r_t_equilib_saved
    if nargin > 1
        varargs = parse_varargin(varargin);
        f = fieldnames(varargs);
        if ~isempty(intersect('savetestdata',f)), savetestdata = varargs.savetestdata;  end;
        invals.lat = lat;
        invals.varargs = varargin;  invals.varargs.savetestdata = false;
    end
    % ###############################################################


const=r_t_getconsts;

g=9.81;            % m/s^2;
erad=6365;         % km
earthmoond=3.84e5; % km
Mmoon=7.38e22;     % kg
Mearth=5.977e24;   % kg

G=3/4*Mmoon^2*(erad/earthmoond)^3/Mearth;

jk=isfinite(const.doodsonamp);


freq=const.freq(jk);
name=const.name(jk,:);

slat=sin(lat(:)'*pi/180);
clat=cos(lat(:)'*pi/180);

G1=zeros(6,length(clat));

% Latitude dependence of amplitude for various species -
% + for A, -for B (from Godin, 1972).

G1(3+0,:)=    0.5*G*(1-3*slat.^2);
G1(3-1,:)=      2*G*slat.*clat;
G1(3+1,:)= .72618*G*clat.*(1-5*slat.^2);
G1(3-2,:)=2.59808*G*slat.*clat.^2;
G1(3+2,:)=        G*clat.^2;
G1(3+3,:)=        G*clat.^3;


amp=abs(const.doodsonamp(jk,ones(1,length(clat)))/g.*G1(const.doodsonspecies(jk)+3,:));


if nargout==0,
    
    plot(24*[freq,freq]',[0;1]*amp');
    
    cnam=cellstr(name);
    for k=1:length(cnam),
        cnam{k}=deblank(cnam{k});
        ff=min([find(abs(cnam{k}(2:end))>=abs('0') & abs(cnam{k}(2:end))<=abs('9'))+1,length(cnam{k})+1]);
        cnam{k}=[ cnam{k}(1:ff-1) '_{' cnam{k}(ff:end) '}'];
    end;
    text(freq*24,amp,cnam,'vertical','bottom','horiz','center','fontangle','italic','fontweight','bold',...
        'clip','on');
    xlabel('Frequency (cpd)');
    ylabel('Potential');
    set(gca,'tickdir','out');
    
end;
    % ##################################################################
    % Save test data
    if savetestdata && ~r_t_equilib_saved

        r_t_equilib_saved = true;
        %[name,freq,amp]
        outvals.name = name;
        outvals.freq = freq;
        outvals.amp = amp;
        writetestdata('r_t_equilib',invals,outvals);

    end
    % ##################################################################
end



