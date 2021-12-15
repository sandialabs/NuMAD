%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Part of the SNL NuMAD Toolbox                    
%  Developed by Sandia National Laboratories Wind Energy Technologies 
%              See license.txt for disclaimer information             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MaterialDef < handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ``MaterialDef``  A class definition for blade materials.
%
% Examples: 
% 
%     ``mat_obj = ComponentDef();``
%
% See also ``xlsBlade``, ``BladeDef``, ``BladeDef.addMaterial``
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        name            % User selected name of the material
        type            % Two options: ‘isotropic’ or ‘orthotropic’
        layerthickness  % Layer thickness [mm]
        ex              % Longitudinal elastic modulus [Pa]
        ey              % Transverse elastic modulus [Pa]
        ez              % Through-the-thickness elastic modulus in the principal material coordinates [Pa]
        gxy             % In-plane shear modulus [Pa]
        gyz             % Transverse shear modulus [Pa]
        gxz             % Transverse shear modulus [Pa]
        prxy            % In-plane Poisson ratio [ ]
        pryz            % Transverse Poisson ratio [ ]
        prxz            % Transverse Poisson ratio [ ]
        density         % Cured mass density [kg/m2]
        drydensity      % Density of fabric
        uts             % 1 × 3 array of ultimate tensile strength design values. Sequence: SL , ST, Sz, 1 × 1 for isotropic.
        ucs             % 1 × 3 array of ultimate compressive strength design values. Sequence: SL , ST, Sz, 1 × 1 for isotropic.
        uss             % 1 × 3 array of ultimate shear strength design values. Sequence: SLT , STz, SLz, 1 × 1 for isotropic.
        xzit            % Lz tensile inclination parameter for Puck failure index
        xzic            % Lz compressive inclination parameter for Puck failure index
        yzit            % Tz tensile inclination parameter for Puck failure index
        yzic            % Tz compressive inclination parameter for Puck failure index
        g1g2            % Fracture toughness ratio between GI (mode I) and GII (mode II) [ ]
        alp0            % Fracture angle under pure transverse compression [degrees]
        etat            % Transverse friction coefficient for Larc [ ]
        etal            % Longitudinal friction coefficient for Larc [ ]
        m               % Fatigue slope exponent [ ]
        gamma_mf        % from DNL-GL standard, fatigue strength reduction factor
        gamma_ms        % from DNV-GL standard, short term strength reduction factor
        reference       % 
    end
    
    methods
        function obj = MaterialDef(mat)
            if nargin > 0
                obj.name           = mat.name;
                obj.type           = mat.type;
                obj.layerthickness = mat.layerthickness;
                obj.ex             = mat.ex;
                obj.ey             = mat.ey;
                obj.ez             = mat.ez;
                obj.gxy            = mat.gxy;
                obj.gyz            = mat.gyz;
                obj.gxz            = mat.gxz;
                obj.prxy           = mat.prxy;
                obj.pryz           = mat.pryz;
                obj.prxz           = mat.prxz;
                obj.density        = mat.density;
                obj.drydensity     = mat.drydensity;
                obj.uts            = mat.uts;
                obj.ucs            = mat.ucs;
                try
                    obj.uss            = mat.uss;
                catch
                    %Do nothing
                end
                try
                    obj.xzit           = mat.xzit;
                    obj.xzic           = mat.xzic;
                    obj.yzit           = mat.yzit;
                    obj.yzic           = mat.yzic;
                catch
                    %Do nothing
                end
                try
                    obj.g1g2           = mat.g1g2;                    
                    obj.alp0           = mat.alp0;
                    obj.etat           = mat.etat;
                    obj.etal           = mat.etal;
                catch
                    %Do nothing
                end
                try
                        obj.m          = mat.m;
                catch
                        %Do nothing
                end
                obj.reference      = mat.reference;
            end
        end
        
        function set.density(obj,density)
            if density < 0
                throw(MException('MaterialDef:DensityRangeError', ...
                                 'Density must be a positive number.'));
            else
                obj.density = density;
            end
        end
        
        function set.prxy(obj,prxy)
            if prxy < 0 || prxy > 1
                throw(MException('MaterialDef:PoissonRangeError', ...
                                 'Poisson''s Ratio must be in range (0,1).'));
            else
                obj.prxy = prxy; 
            end
        end          
        
    end
    
    
end % classdef