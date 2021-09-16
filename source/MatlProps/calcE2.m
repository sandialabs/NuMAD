function y = calcE2(Vf,E2f,Em,Eta,NU12f,NU21f,NUm,E1f,xi1,p)
%E2 This function returns Young’s modulus in the
% transverse direction. Its input are nine values:
% Vf - fiber volume fraction
% E2f - transverse Young’s modulus of the fiber
% Em - Young’s modulus of the matrix
% Eta - stress-partitioning factor
% NU12f - Poisson’s ratio NU12 of the fiber
% NU21f - Poisson’s ratio NU21 of the fiber
% NUm - Poisson’s ratio of the matrix
% E1f - longitudinal Young’s modulus of the fiber
% xi1 - Halpin-Tsai parameter; reinforcing efficiency factor for transverse
%       loading; varies from 1 to 2
% p - parameter used to determine which equation to use:
% p = 1 - use equation (3.4) (not typically accurate)
% p = 2 - use equation (3.9) (use if Eta is known; typically 0.4-0.6)
% p = 3 - use equation (3.10)(use if Eta is not known but Nu12,Nu21,Num,E1f are known)
% p = 4 - use eq(3.97)[2]    (use if Eta is not known but Num is known)
% p = 5 - use eq(3.98)[2]    (use if Halpin-Tsai parameter, xi1, is known)
% Use the value zero for any argument not needed
% in the calculations.
%
% See Chapter 3 of
%     George Z. Voyiadjis and Peter I. Kattan. Mechanics of Composite
%       Materials with MATLAB. Springer-Verlag Berlin Heidelberg, 2005.
% for more information
%
% Ref[2]: Isaac M. Daniel and Ori Ishai.  Engineering Mechanics of Composite
%           Materials. Oxford University Press, 1994.

if (Vf==0 || E2f==0 || Em==0),error('Vf,E2f,Em must be nonzero'),end
Vm = 1 - Vf;
if p == 1
    y = 1/(Vf/E2f + Vm/Em);
elseif p == 2
    if Eta==0,error('Eta must be nonzero'),end
    y = 1/((Vf/E2f + Eta*Vm/Em)/(Vf + Eta*Vm));
elseif p == 3
    if (NU12f==0 || NU21f==0 || NUm==0),error('NU''s must be nonzero'),end
    if E1f==0,error('E1f must be nonzero'),end
    deno = E1f*Vf + Em*Vm;
    etaf = (E1f*Vf + ((1-NU12f*NU21f)*Em + NUm*NU21f*E1f)*Vm) /deno;
    etam = (((1-NUm*NUm)*E1f - (1-NUm*NU12f)*Em)*Vf + Em*Vm) /deno;
    y = 1/(etaf*Vf/E2f + etam*Vm/Em);
elseif p == 4
    if NUm==0,error('NUm must be nonzero'),end
    Em_prime=Em/(1-NUm^2);
    y = 1/(Vf/E2f + Vm/Em_prime);
elseif p == 5
    if xi1==0,error('xi1 must be known'),end
    eta1=(E2f-Em)/(E2f+xi1*Em);
    y = Em*(1+xi1*eta1*Vf)/(1-eta1*Vf);
end
end