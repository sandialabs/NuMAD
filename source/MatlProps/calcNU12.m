function y = calcNU12(Vf,NU12f,NUm)
%NU12 This function returns Poisson’s ratio NU12
% Its input are three values:
% Vf - fiber volume fraction
% NU12f - Poisson’s ratio NU12 of the fiber
% NUm - Poisson’s ratio of the matrix
% This function uses the simple rule-of-mixtures
% formula of equation (3.3)
%
% See Chapter 3 of 
%     George Z. Voyiadjis and Peter I. Kattan. Mechanics of Composite 
%       Materials with MATLAB. Springer-Verlag Berlin Heidelberg, 2005.
% for more information

Vm = 1 - Vf;
y = Vf*NU12f + Vm*NUm;
end