function writePreCompLayup(precomp,output_file)
%WRITEPRECOMPLAYUP  Write a PreComp layup input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   writePreCompLayup(precomp,'output_file')
%      Created for PreComp v1.00.03
%      Designed for use with NuMAD2PreComp()
%
%      precomp = PreComp data structure; view default structure 
%          in NuMAD2PreComp()
%      output_file = filename of file to be created.
%

% loop on all blade stations
for i=1:length(precomp.station)
    dot=max(strfind(output_file,'.'));
    ext=output_file(dot+1:end);
    fn=output_file(1:dot-1);
    num_output_file=[fn '_' num2str(i) '.' ext];
    
    if ~exist('num_output_file','var')  %print to command window if no output_file given
        fid=1;
    else
        fid=fopen(num_output_file,'wt');   %try to open output_file for Writing in Text mode
        if (fid == -1)
            error('Could not open file "%s"\n',num_output_file);
            return
        end
    end
    
    
    wip(fid,[],'Composite laminae lay-up inside the blade section')
    wip(fid,[],'')
    N_scts=precomp.station(i).lp.n_segments;
    wip(fid,[],'*************************** TOP SURFACE ****************************')
    wip(fid,N_scts,'N_scts(1):  no of sectors on top surface')
    wip(fid,[],'')
    wip(fid,[],'normalized chord location of  nodes defining airfoil sectors boundaries (xsec_node)')
    
    PropTable = precomp.station(i).lp.segmentloc;
    Frmt = '%5.4f   ';
    wiphlst(fid,Frmt,PropTable);

    for j=1:N_scts
        try%ble
        stackid=precomp.station(i).lp.stackid(j);
        wip(fid,[],['...........................................................................'])
        
        wip(fid,[],['Sect_num    no of laminae (N_laminas)          STACK:' precomp.stack( stackid ).name{1} ])
        
        N_layers=precomp.stack(stackid).n_layer;
        PropTable = [j N_layers];
        Frmt = '%i             ';
        wiphlst(fid,Frmt,PropTable);
        
        wip(fid,[],'')
        
        wip(fid,[],'lamina    num of  thickness   fibers_direction  composite_material ID')
        wip(fid,[],'number    plies   of ply (m)       (deg)               (-)')
        wip(fid,[],'lam_num  N_plies    Tply         Tht_lam            Mat_id')
        
        for k=1:N_layers
            PropTable = {k,...
                precomp.stack(stackid).layer(k).n, ...
                precomp.stack(stackid).layer(k).thickness, ...
                precomp.stack(stackid).layer(k).theta, ...
                precomp.stack(stackid).layer(k).matid, ...
                precomp.material.name( precomp.stack(stackid).layer(k).matid )};
            Frmt = {'  %i        ', '%i      ', '%5.3e       ', ...
                '%3.1f           ', '%i ', '(%s)'};
            
            wiptbl(fid,Frmt,PropTable,1);
        end
        catch%ble
            keyboard
        end
    end
    
    wip(fid,[],'')
    wip(fid,[],'')
    N_scts=precomp.station(i).hp.n_segments;
    wip(fid,[],'*************************** BOTTOM SURFACE ****************************')
    wip(fid,N_scts,'N_scts(1):  no of sectors on top surface')
    wip(fid,[],'')
    wip(fid,[],'normalized chord location of  nodes defining airfoil sectors boundaries (xsec_node)')
    
    PropTable = precomp.station(i).hp.segmentloc;
    Frmt = '%5.4f   ';
    wiphlst(fid,Frmt,PropTable);
    
    for j=1:N_scts
        stackid=precomp.station(i).hp.stackid(j);
        wip(fid,[],['...........................................................................'])
        
        wip(fid,[],['Sect_num    no of laminae (N_laminas)          STACK:' precomp.stack( stackid ).name{1} ])
        
        N_layers=precomp.stack(stackid).n_layer;
        PropTable = [j N_layers];
        Frmt = '%i             ';
        wiphlst(fid,Frmt,PropTable);
        
        wip(fid,[],'')
        
        wip(fid,[],'lamina    num of  thickness   fibers_direction  composite_material ID')
        wip(fid,[],'number    plies   of ply (m)       (deg)               (-)')
        wip(fid,[],'lam_num  N_plies    Tply         Tht_lam            Mat_id')
        
        for k=1:N_layers
            PropTable = {k,...
                precomp.stack(stackid).layer(k).n, ...
                precomp.stack(stackid).layer(k).thickness, ...
                precomp.stack(stackid).layer(k).theta, ...
                precomp.stack(stackid).layer(k).matid, ...
                precomp.material.name( precomp.stack(stackid).layer(k).matid )};
            Frmt = {'  %i        ', '%i      ', '%5.3e       ', ...
                '%3.1f           ', '%i ', '(%s)'};
            
            wiptbl(fid,Frmt,PropTable,1);
        end
    end
    
    wip(fid,[],'')
    wip(fid,[],'')
    try
        N_webs=length(precomp.station(i).sw.stackid);
    catch
        N_webs=0;
    end
    wip(fid,[],'**********************************************************************')
    wip(fid,[],'Laminae schedule for webs (input required only if webs exist at this section):')
    
    % write shear web layup information
    for j=1:N_webs
        wip(fid,[],'')
        stackid=precomp.station(i).sw.stackid(j);
        wip(fid,[],['web_num    no of laminae (N_weblams)    Name of stack: ' precomp.stack( stackid ).name{1} ])
        N_layers=precomp.stack(stackid).n_layer;
        PropTable = [j N_layers];
        Frmt = '%i            ';
        wiphlst(fid,Frmt,PropTable);
        wip(fid,[],'')
        wip(fid,[],'lamina    num of  thickness   fibers_direction  composite_material ID')
        wip(fid,[],'number    plies   of ply (m)       (deg)               (-)')
        wip(fid,[],'wlam_num N_Plies   w_tply       Tht_Wlam            Wmat_Id')
        
        for k=1:N_layers
            PropTable = {k,...
                precomp.stack(stackid).layer(k).n, ...
                precomp.stack(stackid).layer(k).thickness, ...
                precomp.stack(stackid).layer(k).theta, ...
                precomp.stack(stackid).layer(k).matid, ...
                precomp.material.name( precomp.stack(stackid).layer(k).matid )};
            Frmt = {'  %i        ', '%i      ', '%5.3e       ', ...
                '%3.1f           ', '%i ', '(%s)'};
            
            wiptbl(fid,Frmt,PropTable,1);
        end
    end
    
    % if there is no shear web, set j=0 in the case that there is a flatback AF 
    if isempty(j), j = 0; end

    % write flat back as shear web layup information
    if precomp.station(i).fb.flag
        wip(fid,[],'')
        stackid=precomp.station(i).fb.stackid;
        wip(fid,[],['web_num    no of laminae (N_weblams)    Name of flat back stack: ' precomp.stack( stackid ).name{1} ])
        N_layers=precomp.stack(stackid).n_layer;
        PropTable = [j+1 N_layers];
        Frmt = '%i            ';
        wiphlst(fid,Frmt,PropTable);
        wip(fid,[],'')
        wip(fid,[],'lamina    num of  thickness   fibers_direction  composite_material ID')
        wip(fid,[],'number    plies   of ply (m)       (deg)               (-)')
        wip(fid,[],'wlam_num N_Plies   w_tply       Tht_Wlam            Wmat_Id')
        
        for k=1:N_layers
            PropTable = {k,...
                precomp.stack(stackid).layer(k).n, ...
                precomp.stack(stackid).layer(k).thickness, ...
                precomp.stack(stackid).layer(k).theta, ...
                precomp.stack(stackid).layer(k).matid, ...
                precomp.material.name( precomp.stack(stackid).layer(k).matid )};
            Frmt = {'  %i        ', '%i      ', '%5.3e       ', ...
                '%3.1f           ', '%i ', '(%s)'};
            
            wiptbl(fid,Frmt,PropTable,1);
        end
    end
    
    if fid~=1, fclose(fid); end
    
end

end


%==========================================================================
%===== FUNCTION DEFINITIONS ===============================================
%==========================================================================
function wip(fid,param,descrip)
% write input file parameter
if ~any(size(param)) && ~ischar(param)
    % do nothing if param = []
    % note: used ~any(size(param)) rather than isempty(param)
    %       so that unset parameters (size = [1 0]) will still
    %       get through to the following elseif statements
elseif ischar(param)
    fprintf(fid,'%-16s ',param);  %output string
elseif isfloat(param)
    if numel(param)==1
        fprintf(fid,'%-16g ',param);  %output single number
    else
        str = sprintf(' %g',param);        %create list of numbers
        str = sprintf('"%s"',str(2:end));  %quote the list of numbers
        fprintf(fid,'%-16s ',str);          %output the quoted list
    end
end
fprintf(fid,'%s\n',descrip);
end

function wiphlst(fid,Frmt,param)
% write input file horizontal list of numbers (all on one line, separated by spaces)
for k = 1:length(param)
    fprintf(fid,Frmt,param(k));
end
fprintf(fid,'\n');
end

function wiptbl(fid,frmt,table,nrows)
for r=1:nrows
    for c=1:length(table)
        col = table{c};
        if iscell(col)
            param = col{r};
        else
            param = col(r);
        end
        fprintf(fid,frmt{c},param);
    end
    fprintf(fid,'\n');
end
end