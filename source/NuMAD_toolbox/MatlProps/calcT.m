function y = calcT(theta)
%T This function returns the transformation matrix T
% given the orientation angle "theta".
% There is only one argument representing "theta"
% The size of the matrix is 3 x 3.
% The angle "theta" must be given in degrees.
%
% See Chapter 5 of 
%     George Z. Voyiadjis and Peter I. Kattan. Mechanics of Composite 
%       Materials with MATLAB. Springer-Verlag Berlin Heidelberg, 2005.
% for more information

m = cosd(theta);
n = sind(theta);
y = [m*m n*n 2*m*n ; n*n m*m -2*m*n ; -m*n m*n m*m-n*n];
end