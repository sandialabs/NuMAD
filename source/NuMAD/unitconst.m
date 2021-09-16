function units = unitconst()
%UNITCONST  Defines conversion constants between unit systems  
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%

% lbf2N = 4.4482216152605;
% in2m = 2.54e-2;
% ft2in = 12;
% slug2kg = 14.593902937206;
% lbfs2in_to_slug = 12;    % lbf-s^2/in -> slug
% 
% ft2m = ft2in*in2m;
% 
% Pa2psi = in2m^2 / lbf2N
% Pa2psf = ft2m^2 / lbf2N
% 
% %kg/m^3 -> slug/ft^3
% kgm3_to_slugft3 = ft2m^3 / slug2kg   
% 
% %kg/m^3 -> lbf-s^2/in^4
% kgm3_to_lbfs2in4 = in2m^3 / (lbfs2in_to_slug * slug2kg)

units.stress = {{'Pa','psf','psi'},...
                [1.0, 2.088543423315013e-002, 1.450377377302092e-004]};  % storage units: Pa
units.density = {{'kg/m^3','slug/ft^3','lbf-s^2/in^4'},...
                 [1.0, 1.940320331979764e-003, 9.357254687402411e-008]};  % storage units: kg/m^3       
units.length = {{'m','ft','in'},...
                [1.0, 1/0.3048, 1/0.0254]};  % storage units: m
            
units.angle = {{'deg','rad'},...
               [1.0, pi/180]};   % storage units: degrees
end