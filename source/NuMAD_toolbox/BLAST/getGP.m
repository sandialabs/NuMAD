function [xi,weight] = getGP(numGP)
%getGP Defines Gauss point information for numerical integration
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [xi,weight] = getGP(numGP)
%                    
%   This function defines gauss point coordinates in a local element frame
%   and the associated weights for Gaussian quadrature numerical
%   integration.
%
%      input:
%      numGP        = number of quad points used for integration
%
%      output:
%      xi           = list of quad point coordinates in local element frame
%      weight       = associated weights for quad point coordinate

    %define Gauss integration points
    if(numGP == 1)
        xi(1) = 0;
        weight(1) = 2.0;
    elseif(numGP == 2)
        xi(1) = -sqrt(1/3);
        xi(2) = sqrt(1/3);
        weight(1) = 1.0;
        weight(2) = 1.0;
    elseif(numGP == 3)
        xi(1) = -sqrt(3/5);
        xi(2) = 0.0;
        xi(3) = sqrt(3/5);
        weight(1) = 5/9;
        weight(2) = 8/9;
        weight(3) = 5/9;
    elseif(numGP == 4)
        xi(1) = sqrt((3.0-2*sqrt(6.0/5.0))/7.0);
        xi(2) = -sqrt((3.0-2*sqrt(6.0/5.0))/7.0);
        xi(3) = sqrt((3.0+2*sqrt(6.0/5.0))/7.0);
        xi(4) = -sqrt((3.0+2*sqrt(6.0/5.0))/7.0);
        
        weight(1) = (18+sqrt(30))/36;
        weight(2) = (18+sqrt(30))/36;
        weight(3) = (18-sqrt(30))/36;
        weight(4) = (18-sqrt(30))/36;
    end

end