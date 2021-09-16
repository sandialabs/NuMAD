function y = calcE1(Vf,E1f,Em)
%E1 This function returns Young’s modulus in the
% longitudinal direction. Its input are three values:
% Vf - fiber volume fraction
% E1f - longitudinal Young’s modulus of the fiber
% Em - Young’s modulus of the matrix
% This function uses the simple rule-of-mixtures formula
% of equation (3.2)
%
% See Chapter 3 of 
%     George Z. Voyiadjis and Peter I. Kattan. Mechanics of Composite 
%       Materials with MATLAB. Springer-Verlag Berlin Heidelberg, 2005.
% for more information

Vm = 1 - Vf;
y = Vf*E1f + Vm*Em;
end