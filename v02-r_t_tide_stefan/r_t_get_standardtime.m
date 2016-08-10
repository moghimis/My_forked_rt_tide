
function [ nobsu,t,centraltime,stime  ] = r_t_get_standardtime(tin,varargin )
    % ###################################################################
    % Setup code for saving test data
    savetestdata = false;
    global r_t_get_standardtime
    if nargin > 1
        varargs = parse_varargin(varargin);
        f = fieldnames(varargs);
        if ~isempty(intersect('savetestdata',f)), savetestdata = varargs.savetestdata;  end;
        invals.tin = tin;
        invals.varargs = varargs;   invals.varargs.savetestdata = false;
    end
    % #####################################################################
    n = length(tin);
    nobsu=n-rem(n-1,2);% makes series odd to give a center point
    tin = tin(1:nobsu);
    stime = datevec(min(tin));
    centraltime = mean([min(tin) max(tin)]);
    t = 24*(tin-centraltime);
    
    % ################################################################
    % Save test data
    if savetestdata && ~r_t_get_standardtime

        r_t_get_standardtime = true;
         % [nameu,fu,ju,namei,fi,jinf,jref]
        outvals.nobsu = nobsu;
        outvals.t = t;
        outvals.centraltime = centraltime; 
        outvals.stime = stime;       
        writetestdata('r_t_get_standardtime',invals,outvals);

    end
    % #################################################################
end


