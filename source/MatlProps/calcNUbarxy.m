function y = calcNUbarxy(A,H)
%NUbarxy This function returns the average laminate
% Poisson’s ratio NUxy. Its input are two arguments:
% A - 3 x 3 stiffness matrix for balanced symmetric
% laminates.
% H - thickness of laminate
a = inv(A);
y = -a(1,2)/a(1,1);
%
% See Chapter 9 of 
%     George Z. Voyiadjis and Peter I. Kattan. Mechanics of Composite 
%       Materials with MATLAB. Springer-Verlag Berlin Heidelberg, 2005.
% for more information

end