function [H,p_H_x,p_H_xx] = calculateHermiteShapeFunctions(xi,xNode,factor)
%calculatHermiteShapeFunctions Calculates Hermite shape functions
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [H,p_H_x,p_H_xx] = calculateHermiteShapeFunctions(xi,xNode,factor)
%                    
%   This function calculates the Hermite shape function and shape 
%   function derivatives. The shape function derivatives are defined with 
%   respect to the physical length domain. 
%
%      input:
%      xi           = guass point values to evaluate shape functions at
%      xNode        = nodal coordinates in physical length domain
%      factor       = factor to multiply slope values of Hermite shape
%                     function and shape function derivatives by.  This is
%                     needed for differences in relations between 
%                     rotations and deflection slopes for Euler-Bernoulli
%                     beam theory.
%
%      output:
%      H            = shape function value at specified gauss points
%      p_H_x        = shape function derivative w.r.t physical length
%                     domain at specified gauss points
%      p_H_xx       = shape function second derivative w.r.t physical
%                     length domain at specified gauss points



h = xNode(2) - xNode(1);
x = xNode(1) + (xi + 1)*h/2;
xbar = x-xNode(1);


H(1) = 1-3*(xbar/h)^2 + 2*(xbar/h)^3;
H(2) = -xbar*(1-xbar/h)^2*factor;
H(3) = 3*(xbar/h)^2 - 2*(xbar/h)^3;
H(4) = -(xbar)*((xbar/h)^2 - xbar/h)*factor;

p_H_x(1) = -6*xbar/(h^2)*(1-xbar/h);
p_H_x(2) = -(1 + 3*(xbar/h)^2 -4*(xbar/h))*factor;
p_H_x(3) = -p_H_x(1);
p_H_x(4) = -xbar/h*(3*xbar/h - 2)*factor;

p_H_xx(1) = -6/h^2*(1-2*xbar/h);
p_H_xx(2) = -2/h*(3*xbar/h-2)*factor;
p_H_xx(3) = -p_H_xx(1);
p_H_xx(4) = -2/h*(3*xbar/h - 1)*factor;


end