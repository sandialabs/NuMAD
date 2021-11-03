function writePreCompMaterial(precomp,output_file)
%WRITEPRECOMPMATERIAL  Write a PreComp material input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   writePreCompMaterial(precomp,'output_file') 
%      Created for PreComp v1.00.03
%      Designed for use with NuMAD2PreComp()
%
%      precomp = PreComp data structure; view default structure 
%          in NuMAD2PreComp()
%      output_file = filename of file to be created.
%

if ~exist('output_file','var')  %print to command window if no output_file given
    fid=1;
else
    fid=fopen(output_file,'wt');   %try to open output_file for Writing in Text mode
    if (fid == -1)
        error('Could not open file "%s"\n',output_file);
        return
    end
end

wip(fid,[],'')
wip(fid,[],'Mat_Id     E1           E2          G12       Nu12     Density      Mat_Name')
wip(fid,[],' (-)      (Pa)         (Pa)        (Pa)       (-)      (Kg/m^3)       (-)')

PropTable = {precomp.material.id, precomp.material.e1, precomp.material.e2, ...
            precomp.material.g12, precomp.material.nu12, precomp.material.dens,...
            precomp.material.name};
Frmt = {'  %i   ', '%5.3e   ', '%5.3e   ', ...
        '%5.3e   ', '%5.3f   ', '%5.3e    ', ...
        '(%s)'};

wiptbl(fid,Frmt,PropTable,precomp.material.id(end));

if fid~=1, fclose(fid); end
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