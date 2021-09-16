function af = readNuMADAirfoil(filename)
%readNuMADAirfoil  Read a NuMAD airfoil file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   af = readNuMADAirfoil(filename) 
%
%   af =  data structure containing the contents of the NuMAD airfiol file.
%   filename = NuMAD airfoil file name
%

% Open the file and read the entire contents
fid = fopen(filename);
if (fid == -1)
    error('Could not open file "%s"',filename);
end
filecontents = fread(fid,inf,'uint8=>char')';
fclose(fid);

[pn,fn,ext] = fileparts(filename);
af.name = fn;

% The following regular expression pattern matches any number of characters
% found between the opening and closing "reference" tags
pattern = '<reference>(.*)</reference>';
t = regexp(filecontents, pattern, 'tokens');
% t is a cell containing a cell array
try
    % jcb: is there a better way to extract the contents of t?
    af.reference = cell2mat(t{1});  
catch me
    af.reference = '';
end

% The following regular expression pattern matches any number of characters
% found between the opening and closing "coords" tags
pattern = '<coords>(.*)</coords>';
t = regexp(filecontents, pattern, 'tokens');
% t is a cell containing a cell array
try
    % jcb: is there a better way to extract the contents of t?
    coord_text = cell2mat(t{1});  
catch me
    error('Airfoil <coords>..</coords> not found in file "%s".  Please check the file format.',filename);
end
% Convert coord_text to floating point; coords is a cell array  
% coords{1} are the x coordinates and coords{2} are the y coordinates
af.coords = textscan(coord_text,'%f %f');

