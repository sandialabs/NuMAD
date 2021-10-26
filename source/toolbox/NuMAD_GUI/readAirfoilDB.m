function [afdb aflist] = readAirfoilDB(folder)
%READAIRFOILDB  Read a database (folder) of NuMAD airfoil files.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [AFDB AFLIST] = readAirfoilDB(FOLDER)
%
%   FOLDER is a string containing the name of the database folder to be 
%   read.  The pathname may be absolute or relative.
%
%   AFDB is an array of structures with fields:
%       .name = filename, but without the extension
%       .reference = reference block of the airfoil file
%       .xy = (x,y) coords of airfoil
%  Database entry i is accessed as AFDB(i)
%  The k-th x coord of entry i is AFDB(i).xy(k,1)
%
%  AFLIST is a cell array of airfoil names in the database.

%       .TEtype = trailing edge type as determined by simple rules
%       .LE = index of the leading edge coordinate
%       .xn = chord parameter used for interpolation
%       .sn = arc length parameter used for interpolation


% ASSUMPTIONS / REQUIREMENTS
%    It is assumed that an airfoil's leading edge is located
% at (0,0) and it's trailing edge is located at (1,0), 
% but note there might not necessary be coordinate pairs at 
% these locations due to point spacing or finite TE thickness.
%    This script will test some assumptions, but will not
% re-orient or re-normalize the airfoil -- the user must
% provide airfoil data in the appropriate format.

%folder = 'airfoils';
files = dir(folder);
kdb = 0;
for kf=1:numel(files)
    if files(kf).isdir == 0   % the file is not a directory
        filename = fullfile(folder,files(kf).name);
        try
            af = readAirfoil(filename);
            kdb = kdb + 1;  % this line is not executed if there is an
                            % error reading the airfoil file
        catch ME
            % notify the user of the error, then skip to next file 
            fprintf('%s: %s\n',mfilename,ME.message);
            continue
        end
    else
        continue  % the file is a directory, so skip it
    end
    aflist{kdb,1} = af.name;
    afdb(kdb) = af;
    
%     xy = af.coords;
%     afdb(kdb).name = af.name;
%     afdb(kdb).reference = af.reference;
%     afdb(kdb).xy = xy;
    

%     
%     
%     % find index of leading edge and make HP x-coord negative
%     xn = xy(:,1);
%     [~,LE] = min(xn);
%     if xy(LE,2) >= 0
%         xn(1:LE-1) = -1*xn(1:LE-1);
%     else
%         xn(1:LE) = -1*xn(1:LE);
%     end
%     afdb(kdb).LE = LE;
%     afdb(kdb).xn = xn;
%     
%     % calculate arc length of xy points clockwise from trailing edge
%     n_points=size(xy,1);
%     sn=zeros(n_points,1);
%     for i=2:n_points
%         sn(i) = hypot(xy(i,1)-xy(i-1,1),xy(i,2)-xy(i-1,2)) + sn(i-1);
%     end
%     afdb(kdb).sn = sn - sn(LE);  % arc length now measured from LE coord pair
    
end