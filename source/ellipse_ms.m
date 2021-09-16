function MS = ellipse_ms(c,t,wt,rho,E,G)
% ELLIPSE_MS  compute mass and stiffness matrices for an elliptical, thin-
%       walled section
%                           Under construction
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
% MS = ellipse_ms(c,t,wt,rho,E,G)
% 
%   c = full length of major axis
%   t = 
%   wt = 
%   rho = density
%   E = young's modulus
%   G = shear modulus
% 
%   MS(1,1,1)=rho*A;
%   MS(2,2,1)=rho*A;
%   MS(3,3,1)=rho*A;
%   MS(4,4,1)=rho*Ix;
%   MS(5,5,1)=rho*Iy;
%   MS(6,6,1)=rho*J;
%   MS(1,1,2)=E*A;
%   MS(2,2,2)=E*Ix;
%   MS(3,3,2)=E*Iy;
%   MS(4,4,2)=G*J;
%   MS(5,5,2)=G*A;
%   MS(6,6,2)=G*A;
 
b=c/2;
a=t/2;

K1=0.2464+0.002222*(a/b+b/a);
K2=0.1349+0.1279*(a/b)-0.01284*(a/b)^2;
K3=0.1349+0.1279*(b/a)-0.01284*(b/a)^2;

eta=((a-b)/(a+b))^2;

Am=pi()*a*b;
U=pi()*(a+b)*(1+3*eta/(10+(4-3*eta)^.5));
A=pi()*wt*(a+b)*(1+K1*eta);
Ix=pi()/4*wt*a^2*(a+3*b)*(1+K2*eta)+pi()/16*wt^3*(3*a+b)*(1+K3*eta);
Iy=pi()/4*wt*b^2*(b+3*a)*(1+K3*eta)+pi()/16*wt^3*(3*b+a)*(1+K2*eta);
J=4*Am^2*wt/U;

MS=zeros(7,7,2);

MS(1,1,1)=rho*A;
MS(2,2,1)=rho*A;
MS(3,3,1)=rho*A;
MS(4,4,1)=rho*Ix;
MS(5,5,1)=rho*Iy;
MS(6,6,1)=rho*J;

MS(1,1,2)=E*A;
MS(2,2,2)=E*Ix;
MS(3,3,2)=E*Iy;
MS(4,4,2)=G*J;
MS(5,5,2)=G*A;
MS(6,6,2)=G*A;

end