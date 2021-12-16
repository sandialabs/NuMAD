%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Part of the SNL NuMAD Toolbox                    
%  Developed by Sandia National Laboratories Wind Energy Technologies 
%              See license.txt for disclaimer information             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef StationDef < handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ``StationDef``  A class definition for blade stations.
%
% Examples: 
% 
%	``sta = StationDef();``
% 
%	``sta = StationDef(af);`` where ``af`` = airfoil filename or ``AirfoilDef`` object
%
% See also ``xlsBlade``, ``BladeDef``, ``BladeDef.addStation``
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        airfoil                 % ``AirfoilDef`` object
        spanlocation            % Spanwise location where station is defined [m]
    end
    properties (Dependent = true, SetAccess = private)
        degreestwist
        chord
        percentthick
        coffset
        xyz
    end
    properties (Hidden)
        parent
        hgProfile
    end
    
    methods
        function obj = StationDef(af)
            if ~exist('af','var') || isempty(af)
                obj.airfoil = AirfoilDef();
            else
                if isa(af,'AirfoilDef')
                    obj.airfoil = af;
                elseif ischar(af)
                    obj.airfoil = AirfoilDef(af);
                end
            end
        end
        
        function degreestwist = get.degreestwist(obj)
            degreestwist = interp1(obj.parent.span,obj.parent.degreestwist,obj.spanlocation);
        end
        function chord = get.chord(obj)
            chord = interp1(obj.parent.span,obj.parent.chord,obj.spanlocation);
        end
        function percentthick = get.percentthick(obj)
            percentthick = interp1(obj.parent.span,obj.parent.percentthick,obj.spanlocation);
        end
        function coffset = get.coffset(obj)
%             coffset = interp1(obj.parent.span,obj.parent.coffset,obj.spanlocation);
            coffset = obj.airfoil.maxthick;
        end
        
        function xyz = get.xyz(obj)
            twistFlag = -1;
            tratio = obj.percentthick / obj.airfoil.percentthick;
            thick = obj.airfoil.thickness * tratio;
            hp = obj.airfoil.camber - 0.5*thick;
            lp = obj.airfoil.camber + 0.5*thick;
            c = obj.airfoil.c;
            x = [ c(end); flipud(c);   c(2:end);  c(end)];
            y = [hp(end); flipud(hp); lp(2:end); lp(end)];
            x = (x - obj.coffset) * obj.chord * twistFlag;
            y = (y              ) * obj.chord;
            twist = twistFlag * obj.degreestwist * pi/180;
            xyz = zeros(length(x),3);
            xyz(:,1) = cos(twist) * x - sin(twist) * y;
            xyz(:,2) = sin(twist) * x + cos(twist) * y;
            xyz(:,3) = obj.spanlocation;
        end
        
        function updateProfiles(obj)
            N = numel(obj);
            for k=1:N
                xyz = obj(k).xyz;
                if isempty(obj(k).hgProfile) || ~ishandle(obj(k).hgProfile)
                    obj(k).hgProfile = line(0,0,0);
                end
                set(obj(k).hgProfile,'XData',xyz(:,3),...
                                     'YData',xyz(:,1),...
                                     'ZData',xyz(:,2));
            end
        end
        
        function delete(obj)
            if ishandle(obj.hgProfile)
                delete(obj.hgProfile);
            end
        end
    end
end