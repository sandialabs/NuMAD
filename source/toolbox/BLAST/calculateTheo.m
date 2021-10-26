function [Theo] = calculateTheo(k)
%calculateTheo Calculates Theodorsen function value for unsteady aero
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [Theo] = calculateTheo(k)
%
%   This function accepts a reduced frequency value and calculates the
%   complex Theodorsen function value.
%
%   input:
%   k     = reduced frequency
%
%   output:
%   Theo  = Theodorsen  constant value

% Approx R.T. Jones 1938
% Theo = 1.0 - 0.165/(1-0.0455/k*1i) - 0.335/(1-0.3/k*1i);
% Exact
Theo = conj(1i.*besselh(1,k)./(besselh(0,k)+1i.*besselh(1,k)));


if(isinf(k))
    Theo = 0.5; % Was incorrectly 1 before. solution approaches 0.5 at k = infinity
end

if(isnan(k))
    Theo = 1; % nan returned for k = 0.
end

end

