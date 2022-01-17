 function BladeDef_to_NuMADfile(blade,numad_file,matdb_file,airfoils_path)
%BladeDef_to_NuMADfile   Construct NuMAD files from BladeDef.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%
%   Usage:
%     BladeDef_to_NuMADfile(blade,'BladeDef_test.nmd','MatDBsi.txt');
%
%     BladeDef_to_NuMADfile(blade,'folder/BladeDef_test.nmd',...
%                                 'folder/MatDBsi.txt','folder/airfoils/');
%
%   See also xlsBlade, BladeDef

[nmd_pathstr, nmd_filename] = fileparts(numad_file);
matdb_pathstr = fileparts(matdb_file);

if ~isempty(nmd_pathstr) && exist(nmd_pathstr,'dir')~=7
    error('Path does not exist: %s\n',nmd_pathstr);
end

if ~isempty(matdb_pathstr) && exist(matdb_pathstr,'dir')~=7
    error('Path does not exist: %s\n',nmd_pathstr);
end


[data,matdb]=dataForNuMADobjects(blade,airfoils_path);

%Add NuMAD GUI specific data here:
     data.active.list = {'**UNSPECIFIED**'};
%     % ansys ==============================================================
    data.ansys.BoundaryCondition     = 'cantilevered';
    data.ansys.ElementSystem         = '181';
    data.ansys.MultipleLayerBehavior = 'distinct';
    data.ansys.meshing               = 'elementsize';
    data.ansys.smartsize             = 5;
    % make the element size 5 percent of max chord
    data.ansys.elementsize           = 0.05*max(blade.chord);
    data.ansys.shell7gen             = 1;
    data.ansys.dbgen                 = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Write the NuMAD files          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     assignin('base','data',data);
%     assignin('base','matdb',matdb);
writeNuMADinput(data,numad_file);
writeMatDB(matdb,matdb_file);
end