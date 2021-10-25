function [N,p_N_x,Jac] = calculateShapeFunctions(elementOrder,xi,x)
%calculateShapeFunctions Calculates Lagrange shape functions
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [N,p_N_x,Jac] = calculateShapeFunctions(elementOrder,xi,x)
%                    
%   This function calculates the Lagrange shape function, shape 
%   function derivative, and Jacobian to map between the local element 
%   domain and physical length of the element. The shape function 
%   derivative is defined with respect to the physical length domain. The
%   shape functions may be linear or quadratic in order.
%
%      input:
%      elementOrder = order of element: 1 linear, 2 quadratic
%      xi           = guass point values to evaluate shape functions at
%      x            = nodal coordinates in physical length domain
%
%      output:
%      N            = shape function value at specified gauss points
%      p_N_x        = shape function derivative w.r.t physical length
%                     domain at specified gauss points
%      Jac          = Jacobian for mat between local element domain and
%                     physical length domain.

% N shape function
% p_N_xi partial derivative of shape function w.r.t. xi

%Linear interpolation functions
if(elementOrder == 1)
   
    N(1) = 0.5*(1.0 - xi);
    N(2) = 0.5*(1.0 + xi);
    
    p_N_xi(1) = -0.5;
    p_N_xi(2) = 0.5;
    
end

%Quadratic interpolation functions
if(elementOrder == 2)
   N(1) = 0.5*(xi-1.0)*xi;
   N(2) = 1.0-xi^2;
   N(3) = 0.5*(xi+1.0)*xi;
   
   p_N_xi(1) = xi - 0.5;
   p_N_xi(2) = -2.0*xi;
   p_N_xi(3) = xi + 0.5;
end

numNodesPerEl = length(N);
Jac=0.0;
for i=1:numNodesPerEl
   Jac = Jac + p_N_xi(i)*x(i);
end

for i=1:numNodesPerEl
   p_N_x(i) = p_N_xi(i)/Jac; 
end
end