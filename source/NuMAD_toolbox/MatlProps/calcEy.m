function y = calcEy(E1,E2,NU21,G12,theta)
%Ey This function returns the elastic modulus
% along the y-direction in the global
% coordinate system. It has five arguments:
% E1 - longitudinal elastic modulus
% E2 - transverse elastic modulus
% NU21 - Poisson’s ratio
% G12 - shear modulus
% theta - fiber orientation angle
% The angle "theta" must be given in degrees.
% Ey is returned as a scalar
%
% See Chapter 6 of 
%     George Z. Voyiadjis and Peter I. Kattan. Mechanics of Composite 
%       Materials with MATLAB. Springer-Verlag Berlin Heidelberg, 2005.
% for more information

m = cosd(theta);
n = sind(theta);
denom = m^4 + (E2/G12 - 2*NU21)*n*n*m*m + (E2/E1)*n^4;
y = E2/denom;
end