function y = calcG23(Vf,Gf,Gm,NUm)
%G23 This function returns the shear modulus G23
% Its input are five values:
% Vf - fiber volume fraction
% Gf - shear modulus G of the fiber
% Gm - shear modulus of the matrix
% NUm - Poisson's ratio of the matrix
%
% See Eqn.4.31 of
%     Ever J. Barbero, Introduction to Composite Material Design. Taylor &
%           Francis, 1999.
% for more information
eta23=(3-4*NUm+Gm/Gf) / (4*(1-NUm));
y = Gm * (Vf+eta23*(1-Vf)) / (eta23*(1-Vf)+Vf*Gm/Gf);
end