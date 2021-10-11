function y = calcEx(E1,E2,NU12,G12,theta)
%Ex This function returns the elastic modulus
% along the x-direction in the global
% coordinate system. It has five arguments:
% E1 - longitudinal elastic modulus
% E2 - transverse elastic modulus
% NU12 - Poisson’s ratio
% G12 - shear modulus
% theta - fiber orientation angle
% The angle "theta" must be given in degrees.
% Ex is returned as a scalar
%
% See Chapter 6 of 
%     George Z. Voyiadjis and Peter I. Kattan. Mechanics of Composite 
%       Materials with MATLAB. Springer-Verlag Berlin Heidelberg, 2005.
% for more information

m = cosd(theta);
n = sind(theta);
denom = m^4 + (E1/G12 - 2*NU12)*n*n*m*m + (E1/E2)*n^4;
y = E1/denom;
end