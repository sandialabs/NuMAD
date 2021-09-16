function y = calcEbarx(A,H)
%Ebarx This function returns the average laminate modulus
% in the x-direction. Its input are two arguments:
% A - 3 x 3 stiffness matrix for balanced symmetric
% laminates.
% H - thickness of laminate
a = inv(A);
y = 1/(H*a(1,1));
%
% See Chapter 9 of 
%     George Z. Voyiadjis and Peter I. Kattan. Mechanics of Composite 
%       Materials with MATLAB. Springer-Verlag Berlin Heidelberg, 2005.
% for more information

end