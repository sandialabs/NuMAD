%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Part of the SNL NuMAD Toolbox                    
%  Developed by Sandia National Laboratories Wind Energy Technologies 
%              See license.txt for disclaimer information             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef StackDef
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ``StackDef``  A class definition for a stack of composite layers.
%   
%
% Examples: 
% 
%     ``stack = StackDef();``
% 
% See also ``xlsBlade``, ``BladeDef``, ``BladeDef.updateBOM``
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        name = ''   % String: Name of the stack or composite material used by NuMAD, e.g. ``'000000_HP_LE_PANEL'``
        plygroups   % Array: Array of ``ply`` structures, one for each ``ply``, ``plygroups = struct('component','','materialid',[],'thickness',[],'angle',[],'nPlies',[])``
        indices     % Indices of stack, ``[in board station, out board station, 1st kepoint, 2nd keypoint]``, e.g. ``[ibSta,obSta,keypt1,keypt2]``
    end
    
    methods
        function obj = StackDef()
            if nargin > 0
                
            end
        end
        
        function obj = addply(obj,ply)
            %  This method adds ``ply`` to Stack
            % 
            % Example:
            % 
            %     ``StackDef.addply(ply)``
            %             
            % where 
            %             
            % ``ply = struct('component','',...  % parent component``
            %                         
            %               ``'materialid',[],... % materialid of ply``
            %             
            %               ``'thickness',[],...  % thickness [mm] of single ply``
            %             
            %               ``'angle',[],...      % ply angle``
            %             
            %               ``'nPlies',[]);       % number of plies`` 
            %             
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