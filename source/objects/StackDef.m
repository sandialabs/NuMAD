classdef StackDef
%StackDef  A class definition for a stack of composite layers.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
% StackDef_obj
%   .name        Name of the stack / composite material used by NuMAD
%   .plygroups   Array of ply structures (see below), one for each ply
%   .indices     Indices of stack [ibSta,obSta,keypt1,keypt2]
%
%  ply = struct('component','',...  % parent component
%               'materialid',[],... % materialid of ply
%               'thickness',[],...  % thickness [mm] of single ply
%               'angle',[],...      % ply angle
%               'nPlies',[]);       % number of plies
%
%   Usage examples: 
%     stack = StackDef();
%     stack.name = '000000_HP_LE_PANEL';
%     stack.addPly(ply);
%
%   See also xlsBlade, BladeDef, BladeDef.updateBOM
    properties
        name
        plygroups   % plygroups = struct('component','','materialid',[],'thickness',[],'angle',[],'nPlies',[])
        indices
    end
    
    methods
        function obj = StackDef()
            if nargin > 0
                
            end
        end
        
        function obj = addply(obj,ply)
            %  StackDef.addply(ply)
            %  ply = struct('component','',...  % parent component
            %               'materialid',[],... % materialid of ply
            %               'thickness',[],...  % thickness [mm] of single ply
            %               'angle',[],...      % ply angle
            %               'nPlies',[]);       % number of plies
            if ( (numel(obj.plygroups) > 0) ...
                    && isequal(ply.component,obj.plygroups(end).component) ...
                    && isequal(ply.angle ,obj.plygroups(end).angle) )
                obj.plygroups(end).nPlies = obj.plygroups(end).nPlies + 1;
            else
                if numel(obj.plygroups) > 0
                    obj.plygroups(end+1) = ply;
                else
                    obj.plygroups = ply;
                end
            end
        end
    end
end