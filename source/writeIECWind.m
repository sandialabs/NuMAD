function writeIECWind(iecwind,output_file)
%WRITIECWIND  Write an IECWind input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   writeIECWind(iecwind,'output_file') 
%      Created for IECWind v5.01.01
%
%      iecwind = IECWind data structure; view default structure 
%          in this source code
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

wip(fid,[],'!HEADER: Sample input file for IECWind version 5.01.01');
wip(fid,[],'!Output file parameters');
wip(fid,iecwind.SI,'SI UNITS (True=SI or False=ENGLISH)');
wip(fid,iecwind.TStart,'Time for start of IEC transient condition, sec');
wip(fid,[],'!Wind Site parameters');
wip(fid,iecwind.Class,'IEC WIND TURBINE CLASS (1, 2 or 3)');
wip(fid,iecwind.Turb,'WIND TURBULENCE CATEGORY (A, B or C)');
wip(fid,iecwind.Slope,'Slope of the wind inflow (IEC specifies between -8 and +8), deg');
wip(fid,iecwind.Shear,'IEC standard used for wind shear exponent');
wip(fid,[],'!Turbine parameters');
wip(fid,iecwind.HubHt,'Wind turbine hub-height, m or ft');
wip(fid,iecwind.RotDia,'Wind turbine rotor diameter, m or ft');
wip(fid,iecwind.CutIn,'Cut-in wind speed, m/s or ft/s');
wip(fid,iecwind.Rated,'Rated wind speed, m/s or ft/s');
wip(fid,iecwind.CutOut,'Cut-out wind speed, m/s or ft/s');
wip(fid,[],'!List of Conditions to generate (one per line)');
wiplst(fid,iecwind.Conditions,'');

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
