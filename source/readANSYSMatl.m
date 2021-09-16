function mpdata = readANSYSMatl(filename)
%readANSYSMatl  Read ANSYS list of material properties
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   output = readANSYSMatl(filename)
%
%   Default input file name is Materials.txt
%
%   Output is a structure with material ID's and corresponding properties.
%

if ~exist('filename','var')
    filename = 'Materials.txt';
end

% Open the file and read the entire contents
fid = fopen(filename);
if (fid == -1)
    error('Could not open file "%s"',filename);
end
filecontents = fread(fid,inf,'uint8=>char')';
fclose(fid);
%assignin('base','filecontents',filecontents);

pat = 'MPDATA,(?<prop>\w*)\s*,\s*(?<matnum>\d+),\s*(?<stloc>\d+),\s*(?<value>[^,]+),';
data = regexp(filecontents,pat,'names');
%assignin('base','data',data);

for k = 1:numel(data)
    matnum = str2double(data(k).matnum);
    propname = data(k).prop;
    value = str2double(data(k).value);
    mpdata(matnum).(propname) = value;
end
%assignin('base','mpdata',mpdata);
end