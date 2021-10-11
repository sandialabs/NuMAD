function y = calcGxy(E1,E2,NU12,G12,theta)
%Gxy This function returns the shear modulus
% Gxy in the global
% coordinate system. It has five arguments:
% E1 - longitudinal elastic modulus
% E2 - transverse elastic modulus
% NU12 - Poisson’s ratio
% G12 - shear modulus
% theta - fiber orientation angle
% The angle "theta" must be given in degrees.
% Gxy is returned as a scalar
%
% See Chapter 6 of 
%     George Z. Voyiadjis and Peter I. Kattan. Mechanics of Composite 
%       Materials with MATLAB. Springer-Verlag Berlin Heidelberg, 2005.
% for more information

m = cosd(theta);
n = sind(theta);
denom = n^4 + m^4 + 2*(2*G12*(1 + 2*NU12)/E1 + 2*G12/E2 - 1)
*n*n*m*m;
y = G12/denom;
end