function writeAirfoilData(polar,output_file)
%writeAirfoilData  Write an AeroDyn airfoil performance input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   writeAirfoilData(polar,'file_name') 
%      Created for AeroDyn V??
%
% Note! This does not work well with multiple airfoil performance tables, yet
%
%      polar = Airfoil performance data structure; view default structure 
%          by utilizing readAirfoilData()
%      output_file = name of file to be created.

if ~exist('output_file','var')  %print to command window if no output_file given
    fid=1;
else
    fid=fopen(output_file,'wt');   %try to open output_file for Writing in Text mode
    if (fid == -1)
        error('Could not open file "%s"\n',output_file);
        return
    end
end

wip(fid,[],polar.title{1});
wip(fid,[],polar.title{2});

wip(fid,polar.NTables,'Number of airfoil tables in this file');
wip(fid,polar.TableIDs,'Table ID parameter');
wip(fid,0,'Stall angle (deg)');
wip(fid,0,'No longer used, enter zero');
wip(fid,0,'No longer used, enter zero');
wip(fid,0,'No longer used, enter zero');
wip(fid,polar.ALPHAL,'Zero Cn angle of attack (deg)');
wip(fid,polar.CNA,'Cn slope for zero lift (dimensionless)');
wip(fid,polar.CNS,'Cn extrapolated to value at positive stall angle of attack');
wip(fid,polar.CNSL,'Cn at stall value for negative angle of attack');
wip(fid,polar.AOD,'Angle of attack for minimum CD (deg)');
wip(fid,polar.CDO,'Minimum CD value');
   
PropTable = {polar.AoA,polar.CL,polar.CD,polar.CM};
Frmt = {'%7.5f  ','%7.5f  ','%7.5f  ','%7.5f  '};
wiptbl(fid,Frmt,PropTable,length(polar.AoA));

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

function wiplst(fid,param,descrip)
% write input file parameter list (one per line)
for k = 1:length(param)
    fprintf(fid,'%-16s ',param{k});  %output string
    if k==1
        fprintf(fid,'%s',descrip);
    end
    fprintf(fid,'\n');
end
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