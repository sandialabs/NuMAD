function af = readAirfoil(filename)
%READAIRFOIL  Read a NuMAD airfoil file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   AF = readAirfoil(FILENAME)
%   Reads a NuMAD-format airfoil file and returns a structure containing
%   the airfoil name, reference, and coordinates
%
%   FILENAME is a string containing the name of the file to be opened.
%   The pathname (absolute or relative) can be included with the filename.
%
%   AF is a structure with fields: 
%     .name = filename, but without the extension
%     .reference = reference block of the airfoil file
%     .coords = Nx2 matrix of x- and y-coordinates
%                            (1st and 2nd columns, respectively)

% Open the file and read the entire contents
fid = fopen(filename);
if (fid == -1)
    error('Could not open file "%s"',filename);
end
filecontents = fread(fid,inf,'uint8=>char')';
fclose(fid);

[~,fn,~] = fileparts(filename);
af.name = fn;

% The following regular expression pattern matches any number of characters
% found between the opening and closing "reference" tags
pattern = '<reference>(.*)</reference>';
t = regexp(filecontents, pattern, 'tokens');
% t is a cell containing a cell array
try
    % jcb: is there a better way to extract the contents of t?
    af.reference = cell2mat(t{1});  
catch ME  %#ok
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
catch ME  %#ok
    error('Airfoil <coords>..</coords> not found in file "%s".  Please check the file format.',filename);
end
% Convert coord_text to floating point; coords is an Nx2 matrix  
% coords(:,1) are the x coordinates and coords(:,2) are the y coordinates
af.coords = cell2mat(textscan(coord_text,'%f %f'));

