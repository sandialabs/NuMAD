function [nargout] = writeNuMADAirfoil(coords,reftext,fname)
%WriteNuMADAirfoil  Write NuMAD airfoil files
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   WriteNuMADAirfoil(coords,reftext,fname)
% 
%     fname - full filename, incl extension, of NuMAD airfoil file to write
%     coords - Nx2 array of airfoil coordinate data.  First column contains
%     x-values, second column contains y-values.  Airfoil coordinates are in
%     order as specified by NuMAD (i.e. trailing edge = (1,0) and leading
%     edge = (0,0)
%     reftext = string representing reference text

fid=fopen(fname,'wt');

fprintf(fid,'<reference>\n%s</reference>\n',reftext);
fprintf(fid,'<coords>\n');
for i=1:length(coords)
    fprintf(fid,'%8.12f\t%8.12f\n',coords(i,:));
end
fprintf(fid,'</coords>');

fclose(fid);
end

