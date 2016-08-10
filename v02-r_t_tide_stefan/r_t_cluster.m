    function ain=r_t_cluster(ain,clusang,varargin)
    % CLUSTER: Clusters angles in rows around the angles in the first 
    % column. CLUSANG is the allowable ambiguity (usually 360 degrees but
    % sometimes 180).
    
    % ##########################################################
    % Setup code for saving test data
    global r_t_cluster_saved
    savetestdata = false;
    if nargin > 2
        varargs = parse_varargin(varargin);
        f = fieldnames(varargs);
        if ~isempty(intersect('savetestdata',f)), savetestdata = varargs.savetestdata;  end;
        invals.ain = ain;
        invals.clusang = clusang;
        invals.varargs = varargs;   invals.varargs.savetestdata = false;
    end
    % ##########################################################
    
    ii=(ain-ain(:,ones(1,size(ain,2))))>clusang/2;
    ain(ii)=ain(ii)-clusang;
    ii=(ain-ain(:,ones(1,size(ain,2))))<-clusang/2;
    ain(ii)=ain(ii)+clusang;
    
    % #########################################################
    % Save test data
    if savetestdata && ~r_t_cluster_saved
       outvals.ain = ain;
       writetestdata('r_t_cluster',invals,outvals);
    end
    % #########################################################
end