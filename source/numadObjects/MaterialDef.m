classdef MaterialDef < handle
%MaterialDef  A class definition for blade materials.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%
%   Usage examples: 
%     mat_obj = ComponentDef();
%     mat_obj = ComponentDef(mat_struct);
%
%   mat_struct:
%     .name
%     .type
%     .layerthickness
%     .ex
%     .ey
%     .ez
%     .gxy
%     .gyz
%     .gxz
%     .prxy
%     .pryz
%     .prxz
%     .density
%     .drydensity
%     .uts
%     .ucs
%     .uss
%     .g1g2
%     .alp0
%     .etat
%     .etal
%     .reference
%
%   See also xlsBlade, BladeDef, BladeDef.addMaterial
    properties
        name
        type
        layerthickness
        ex
        ey
        ez
        gxy
        gyz
        gxz
        prxy
        pryz
        prxz
        density
        drydensity
        uts
        ucs
        uss
        xzit
        xzic
        yzit
        yzic
        g1g2
        alp0
        etat
        etal
        m
        reference
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