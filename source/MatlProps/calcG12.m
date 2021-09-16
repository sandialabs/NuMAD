function y = calcG12(Vf,G12f,Gm,EtaPrime,p)
%G12 This function returns the shear modulus G12
% Its input are five values:
% Vf - fiber volume fraction
% G12f - shear modulus G12 of the fiber
% Gm - shear modulus of the matrix
% EtaPrime - shear stress-partitioning factor
% p - parameter used to determine which equation to use:
% p = 1 - use equation (3.5)  (not typically accurate)
% p = 2 - use equation (3.13) (use if Etaprime is known)
% p = 3 - use equation (3.14) (use if etaprime=0.6; agrees with elasticity)
% Use the value zero for any argument not needed
% in the calculations.
%
% See Chapter 3 of 
%     George Z. Voyiadjis and Peter I. Kattan. Mechanics of Composite 
%       Materials with MATLAB. Springer-Verlag Berlin Heidelberg, 2005.
% for more information

Vm = 1 - Vf;
if p == 1
    y = 1/(Vf/G12f + Vm/Gm);
elseif p == 2
    if EtaPrime==0,error('EtaPrime must be nonzero'),end
    y = 1/((Vf/G12f + EtaPrime*Vm/Gm)/(Vf + EtaPrime*Vm));
elseif p == 3
    y = Gm*((Gm + G12f) - Vf*(Gm - G12f))/((Gm + G12f) + Vf*(Gm - G12f));
end
end