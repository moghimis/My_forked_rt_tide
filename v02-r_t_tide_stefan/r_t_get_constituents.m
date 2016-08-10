%--------------------------------------------------
function [nameu,fu,ju,namei,fi,jinf,jref]=r_t_get_constituents(t,minres,constit,...
                                         shallow,infname,infref,centraltime,varargin)
    % [name,freq,kmpr]=constituents(minres,infname) loads tidal constituent
    % table (containing 146 constituents), then picks out only the '
    % resolvable' frequencies (i.e. those that are MINRES apart), base on 
    % the comparisons in the third column of constituents.dat. Only 
    % frequencies in the 'standard' set of 69 frequencies are actually used.
    % Also return the indices of constituents to be inferred.

    % If we have the mat-file, read it in, otherwise create it and read
    % it in!

    % R Pawlowicz 9/1/01 
    % Version 1.0
    %
    %    19/1/02 - typo fixed (thanks to  Zhigang Xu)

    % Compute frequencies from astronomical considerations.
        % ###################################################################
    % Setup code for saving test data
    savetestdata = false;
    global r_t_get_constituents_saved
    if nargin > 7
        varargs = parse_varargin(varargin);
        f = fieldnames(varargs);
        if ~isempty(intersect('savetestdata',f)), savetestdata = varargs.savetestdata;  end;
        invals.t = t;
        invals.minres = minres;
        invals.constit = constit;
        invals. shallow = shallow;
        invals.infname = infname;
        invals.infref = infref;
        invals.centraltime = centraltime;
        invals.varargs = varargs;   invals.varargs.savetestdata = false;
    end
    % #####################################################################

  
    [const,~,cshallow]=load_constits(t,centraltime,minres,'savetestdata',savetestdata);


    if isempty(constit),
        ju=find(const.df>=minres);
    else
        ju=[];
        for k=1:length(constit)
            j1=strmatch(constit{k},const.name);
            if isempty(j1),
                disp(['Can''t recognize name ' constit(k,:) ' for forced search']);
            elseif j1==1,
                disp('*************************************************************************');
                disp('Z0 specification ignored - for non-tidal offsets see ''secular'' property');
                disp('*************************************************************************');
            else
                ju=[ju;j1]; %#ok<AGROW>
            end;
        end;
        [~,II]=sort(const.freq(ju)); % sort in ascending order of frequency.
        ju=ju(II);
    end;


    disp(['   number of standard constituents used: ',int2str(length(ju))])

    if ~isempty(shallow),
     for k=1:size(shallow,1),
       j1=strmatch(shallow(k,:),const.name);
       if isempty(j1),
         disp(['Can''t recognize name ' shallow(k,:) ' for forced search']);
       else
         if isnan(const.ishallow(j1)),
           disp([shallow(k,:) ' Not a shallow-water constituent']);
         end;
         disp(['   Forced fit to ' shallow(k,:)]);
         ju=[ju;j1];
       end;
     end;

    end;

    nameu=const.name(ju,:);
    fu=const.freq(ju);


    % Check if neighboring chosen constituents violate Rayleigh criteria.
    jck=find(diff(fu)<minres);
    if ~isempty(jck)
       disp('  Warning! Following constituent pairs violate Rayleigh criterion');
       for ick=1:length(jck);
       disp(['     ',nameu(jck(ick),:),'  ',nameu(jck(ick)+1,:)]);
       end;
    end

    % For inference, add in list of components to be inferred.

    fi=[];namei=[];jinf=[];jref=[];
    if ~isempty(infname),
      fi=zeros(size(infname,1),1);
      namei=zeros(size(infname,1),4);
      jinf=zeros(size(infname,1),1)+NaN;
      jref=zeros(size(infname,1),1)+NaN;

      for k=1:size(infname,1),
       j1=strmatch(infname(k,:),const.name);
       if isempty(j1),
         disp(['Can''t recognize name' infname(k,:) ' for inference']);
       else
        jinf(k)=j1;
        fi(k)=const.freq(j1);
        namei(k,:)=const.name(j1,:);
        j1=strmatch(infref(k,:),nameu);
        if isempty(j1),
          disp(['Can''t recognize name ' infref(k,:) ' for as a reference for inference']);
        else
          jref(k)=j1;
          fprintf(['   Inference of ' namei(k,:) ' using ' nameu(j1,:) '\n']);
        end;
       end;
      end;    
      jinf(isnan(jref))=NaN;
    end;
    
    % ################################################################
% Save test data
if savetestdata && ~r_t_get_constituents_saved
   
    r_t_get_constituents_saved = true;
     % [nameu,fu,ju,namei,fi,jinf,jref]
    outvals.nameu = nameu;
    outvals.fu = fu;    
    outvals.ju = ju;
    outvals.namei = namei;
    outvals.fi = fi;
    outvals.jinf = jinf;
    outvals.jref = jref;
    writetestdata('r_t_get_constituents',invals,outvals);
    
end
% #################################################################
    end
