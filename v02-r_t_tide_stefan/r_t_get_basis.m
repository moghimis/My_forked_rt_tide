%=================================================
function [tc tc_aux] = r_t_get_basis(t,centraltime,fu,varargin )
    % ###################################################################
    % Setup code for saving test data
    savetestdata = false;
    global r_t_get_basis_saved
    if nargin > 3
        varargs = parse_varargin(varargin);
        f = fieldnames(varargs);
        if ~isempty(intersect('savetestdata',f)), savetestdata = varargs.savetestdata;  end;
        invals.t = t;
        invals.centraltime = centraltime;
        invals.fu = fu;
        invals.varargs = varargs;   invals.varargs.savetestdata = false;
    end
    % #####################################################################
    % 333
    tc = [];     %#ok<NASGU>
    tc_aux = [];
    %t = get(self,'time');
    t_aux = [];%get(self,'time_extrema');

    if (size(t,2) > size(t,1)), t = t'; end;
    if (size(t_aux,2) > size(t_aux,1)), t_aux = t_aux'; end;

    t = 24*(t-centraltime);

    tc = [cos((2*pi*t)*fu') sin((2*pi*t)*fu') ];

    if ~isempty(t_aux)
        c = ones(length(t_aux),1)*(2*pi*fu');
        t_aux = 24*(t_aux-centraltime);
        tc_aux = [ -c.*(sin((2*pi*t_aux)*fu'))  c.*(cos((2*pi*t_aux)*fu'))];
    end
    
    % ################################################################
% Save test data
if savetestdata && ~r_t_get_basis_saved
    
    r_t_get_basis_saved = true;
    % [emaj,emin,einc,epha]
    outvals.tc = tc;
    outvals.tc_aux = tc_aux;    
    writetestdata('r_t_get_basis',invals,outvals);
    
end
% #################################################################

end

