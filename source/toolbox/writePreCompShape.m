function writePreCompShape(precomp,output_file)
%WRITEPRECOMPSHAPE  Write a PreComp shape input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   writePreCompShape(precomp,'output_file')
%      Created for PreComp v1.00.03
%      Designed for use with NuMAD2PreComp()
%
%      precomp = PreComp data structure; view default structure 
%          in NuMAD2PreComp()
%      output_file = filename of file to be created.
%

for i=1:length(precomp.station)
    numNodes=length(precomp.station(i).x);
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
    
    wip(fid,numNodes,'       N_af_nodes :no of airfoil nodes, counted clockwise starting')
    wip(fid,[],'                      with leading edge (see users'' manual, fig xx)')
    wip(fid,[],'')
    wip(fid,[],' Xnode      Ynode   !! chord-normalized coordinated of the airfoil nodes')
    
    PropTable = {precomp.station(i).x, precomp.station(i).y};
    
    Frmt = {'%5.4f   ', '%5.4f'};
    
    wiptbl(fid,Frmt,PropTable,numNodes);
    
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