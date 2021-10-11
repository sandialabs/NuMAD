function y = calcNUbaryx(A,H)
%NUbaryx This function returns the average laminate
% Poisson’s ratio NUyx. Its input are two arguments:
% A - 3 x 3 stiffness matrix for balanced symmetric
% laminates.
% H - thickness of laminate
%
% See Chapter 9 of 
%     George Z. Voyiadjis and Peter I. Kattan. Mechanics of Composite 
%       Materials with MATLAB. Springer-Verlag Berlin Heidelberg, 2005.
% for more information

a = inv(A);
y = -a(1,2)/a(2,2);
end