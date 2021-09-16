function y = calcGbarxy(A,H)
%Gbarxy This function returns the average laminate shear
% modulus. Its input are two arguments:
% A - 3 x 3 stiffness matrix for balanced symmetric
% laminates.
% H - thickness of laminate
%
% See Chapter 9 of 
%     George Z. Voyiadjis and Peter I. Kattan. Mechanics of Composite 
%       Materials with MATLAB. Springer-Verlag Berlin Heidelberg, 2005.
% for more information
a = inv(A);
y = 1/(H*a(3,3));
end